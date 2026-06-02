#!/usr/bin/env python3
"""
Build customer-facing WPS order folders from an EPI2ME export.

This workflow:
1. discovers per-barcode EPI2ME files
2. maps barcodes to customer metadata
3. realigns raw FASTQ reads, or raw unmapped BAM reads, to the final consensus FASTA
4. generates per-base CSVs, QC plots, and synthetic AB1 files
5. writes a customer-ready package grouped by order number
6. renders a 2-page PDF report with an Alta-style layout
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import zipfile
from collections import defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Iterable
import xml.etree.ElementTree as ET

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.path import Path as MplPath
from matplotlib.patches import Patch, PathPatch, Rectangle
import numpy as np
import pysam

from align_bam_pipeline import run_fastq_pipeline, run_pipeline
from fastq_to_ab1 import phred_from_ascii, synthesize_chromatogram, write_real_ab1
from generate_report import (
    DEFAULT_MULTIMER_DENOMINATOR,
    MULTIMER_DENOMINATOR_ALL_ELIGIBLE_READS,
    MULTIMER_DENOMINATOR_CHOICES,
    MULTIMER_TOLERANCE_FRACTION,
    READ_LENGTH_DISTRIBUTION_MIN_DISPLAY_BP,
    count_fasta_records,
    generate_report_data,
    read_first_fasta_record,
    y_axis_top_with_headroom,
)


PACKAGE_SUBDIRS = {
    "ab1": "CHROMATOGRAM_FILES_ab1",
    "fasta": "FASTA_FILES",
    "gbk": "GENBANK_FILES",
    "per_base": "PER_BASE_BREAKDOWN",
    "qc": "QC REPORTS",
}

HEADER_ALIASES = {
    "barcode": {"barcode", "barcodeid", "barcodenumber", "barcode_name"},
    "sample_name": {"samplename", "sample", "sampleid", "plasmidname", "name"},
    "sample_id": {"id", "sampleidentifier", "samplecode"},
    "order_number": {"ordernumber", "order", "orderno", "ordernum"},
    "serial_number": {"sn", "serialnumber", "samplenumber", "sample_no"},
    "order_date": {"orderdate", "dateordered"},
    "report_date": {"reportdate"},
    "run_date": {"rundate", "sequencingdate"},
    "wps_date": {"wpsdate"},
    "received_date": {"receiveddate"},
    "customer_name": {"customer", "customername"},
    "concentration": {"concngul", "concentration", "concentrationngul"},
    "size": {"size", "plasmidsize"},
    "note": {"note", "notes", "comments"},
}

THEME = {
    "title": "#10223A",
    "heading": "#003B73",
    "rule": "#005DAA",
    "teal": "#0D8686",
    "table_border": "#052866",
    "table_grid": "#052866",
    "purple": "#51459A",
    "green": "#44D8A4",
    "cyan": "#12A9D4",
    "text": "#23313f",
    "muted": "#5b6b77",
    "line": "#d5dde5",
}

# With qscore now taken directly from BAM base qualities, use a modest ONT-style cutoff.
LOW_CONFIDENCE_QSCORE = 12
SIZE_MISMATCH_TOLERANCE_FRACTION = 0.10
SIZE_MISMATCH_TOLERANCE_BP = 100
INPUT_ROOT = Path("/mnt/c/WPS data")
DEFAULT_OUTPUT_SUBDIR = "output"
DEFAULT_LOGO_PATH = Path(__file__).resolve().with_name("Alta Biotech Logo.jpg")
DEFAULT_ECOLI_REFERENCE_FASTA = Path(__file__).resolve().with_name("E. Coli Genome.fna")
HOST_DNA_MIN_ALIGNED_BP = 1300
HOST_DNA_MIN_ALIGNED_PCT = 91.0
RAW_FASTQ_MIN_RECORDS = 100
RAW_FASTQ_MIN_BASES = 100_000
RAW_FASTQ_SCAN_RECORD_LIMIT = 100_000


def slugify(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("_") or "sample"


def normalize_header(header: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", header.lower())


def normalize_barcode(value: str | int | None) -> str | None:
    """Normalize barcode values like 3, 03, 3.0, barcode3, or barcode03."""
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    explicit = re.fullmatch(r"barcode\s*#?\s*0*(\d+)", text, flags=re.IGNORECASE)
    if explicit:
        return f"barcode{int(explicit.group(1)):02d}"
    numeric = re.fullmatch(r"0*(\d+)(?:\.0+)?", text)
    if numeric:
        return f"barcode{int(numeric.group(1)):02d}"
    return None


def normalize_excel_number_text(value: str | int | float | None) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if not text:
        return ""
    try:
        number = Decimal(text)
    except InvalidOperation:
        return text
    if number == number.to_integral_value():
        return str(int(number))
    return text


def normalize_metadata_date(value: str | int | float | None) -> str:
    text = normalize_excel_number_text(value)
    if not text:
        return ""
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", text):
        return text
    try:
        serial = Decimal(text)
    except InvalidOperation:
        return text
    if serial != serial.to_integral_value():
        return text
    days = int(serial)
    if not 20000 <= days <= 60000:
        return text
    # Excel's 1900 date system, including its historical leap-year offset.
    return (dt.date(1899, 12, 30) + dt.timedelta(days=days)).isoformat()


def normalize_order_number(value: str | int | float | None) -> str:
    text = normalize_excel_number_text(value)
    digits = re.fullmatch(r"\d+", text)
    if digits:
        return text
    return text.strip()


def parse_expected_size_bp(value: str | int | float | None) -> int | None:
    text = normalize_excel_number_text(value)
    if not text or text.strip().lower() in {"unk", "unknown", "na", "n/a", "none"}:
        return None
    cleaned = text.strip().lower().replace(",", "")
    match = re.search(r"(\d+(?:\.\d+)?)\s*(kb|kilobase|kilobases|k)?\b", cleaned)
    if not match:
        return None
    number = float(match.group(1))
    unit = match.group(2)
    if unit in {"kb", "kilobase", "kilobases", "k"}:
        number *= 1000
    return int(round(number))


def lengths_match(observed_bp: int, expected_bp: int) -> bool:
    allowed = max(SIZE_MISMATCH_TOLERANCE_BP, expected_bp * SIZE_MISMATCH_TOLERANCE_FRACTION)
    return abs(observed_bp - expected_bp) <= allowed


def build_wps_sample_name(meta: dict[str, str], row_number: int | None) -> str | None:
    sample_name = meta.get("sample_name")
    sample_id = meta.get("sample_id")
    if not sample_name:
        return None
    if sample_id:
        prefix = f"{row_number:03d}_" if row_number is not None else ""
        return f"{prefix}{sample_id}_{sample_name}"
    return sample_name


def load_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    path = "xl/sharedStrings.xml"
    if path not in zf.namelist():
        return []
    root = ET.fromstring(zf.read(path))
    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    strings = []
    for si in root.findall("x:si", ns):
        parts = []
        for node in si.iter():
            if node.tag.endswith("}t") and node.text:
                parts.append(node.text)
        strings.append("".join(parts))
    return strings


def column_index_from_ref(cell_ref: str) -> int:
    letters = re.match(r"[A-Z]+", cell_ref)
    if not letters:
        return 0
    total = 0
    for ch in letters.group(0):
        total = total * 26 + (ord(ch) - ord("A") + 1)
    return total - 1


def read_xlsx_rows(path: Path) -> list[dict[str, str]]:
    with zipfile.ZipFile(path) as zf:
        sheet_name = "xl/worksheets/sheet1.xml"
        if sheet_name not in zf.namelist():
            raise ValueError(f"Could not find {sheet_name} in {path}")
        shared = load_shared_strings(zf)
        root = ET.fromstring(zf.read(sheet_name))

    ns = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    raw_rows: list[list[str]] = []
    for row in root.findall(".//x:sheetData/x:row", ns):
        cells: dict[int, str] = {}
        max_col = -1
        for cell in row.findall("x:c", ns):
            ref = cell.attrib.get("r", "A1")
            col_idx = column_index_from_ref(ref)
            max_col = max(max_col, col_idx)
            cell_type = cell.attrib.get("t")
            value = ""
            if cell_type == "inlineStr":
                node = cell.find("x:is/x:t", ns)
                value = node.text if node is not None and node.text is not None else ""
            else:
                node = cell.find("x:v", ns)
                raw = node.text if node is not None and node.text is not None else ""
                if cell_type == "s" and raw:
                    value = shared[int(raw)]
                else:
                    value = raw
            cells[col_idx] = value
        if max_col < 0:
            continue
        raw_rows.append([cells.get(idx, "") for idx in range(max_col + 1)])

    header = None
    records: list[dict[str, str]] = []
    for row in raw_rows:
        if header is None:
            if any(cell.strip() for cell in row):
                header = row
            continue
        if not any(cell.strip() for cell in row):
            continue
        row = row + [""] * max(0, len(header) - len(row))
        records.append({header[idx]: row[idx] for idx in range(len(header))})
    return records


def read_table_rows(path: Path) -> list[dict[str, str]]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))
    if suffix == ".tsv":
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    if suffix == ".xlsx":
        return read_xlsx_rows(path)
    raise ValueError(f"Unsupported metadata format: {path}")


def canonicalize_metadata_row(row: dict[str, str]) -> dict[str, str]:
    normalized = {normalize_header(key): (value or "").strip() for key, value in row.items()}
    result: dict[str, str] = {}
    for canonical, aliases in HEADER_ALIASES.items():
        for alias in aliases:
            if alias in normalized and normalized[alias]:
                result[canonical] = normalized[alias]
                break
    if "barcode" in result:
        barcode = normalize_barcode(result["barcode"])
        if barcode is None:
            raise ValueError(f"Could not parse barcode value: {result['barcode']!r}")
        result["barcode"] = barcode
    if "order_number" in result:
        result["order_number"] = normalize_order_number(result["order_number"])
    if "sample_id" in result:
        result["sample_id"] = normalize_excel_number_text(result["sample_id"])
    for date_key in ("order_date", "report_date", "run_date", "received_date", "wps_date"):
        if date_key in result:
            result[date_key] = normalize_metadata_date(result[date_key])
    if "wps_date" in result:
        for date_key in ("order_date", "report_date", "run_date"):
            result.setdefault(date_key, result["wps_date"])
    return result


def metadata_name_score(path: Path) -> int:
    name = normalize_header(path.stem)
    has_wps = "wps" in name
    has_working = "working" in name
    has_sheet = "sheet" in name
    if has_wps and has_working and has_sheet:
        return 3
    if has_working and has_sheet:
        return 2
    if has_wps and has_sheet:
        return 1
    return 0


def resolve_metadata_path(path: Path) -> Path:
    if path.is_dir():
        spreadsheet_candidates = sorted(
            item
            for item in path.iterdir()
            if item.is_file()
            and item.suffix.lower() in {".xlsx", ".csv", ".tsv"}
            and not item.name.startswith("~$")
        )
        candidates = [item for item in spreadsheet_candidates if metadata_name_score(item) > 0]
        if not candidates:
            names = ", ".join(item.name for item in spreadsheet_candidates)
            detail = f" Found: {names}" if names else ""
            raise ValueError(
                f"No WPS Working Sheet metadata .xlsx, .csv, or .tsv file found in {path}.{detail}"
            )
        best_score = max(metadata_name_score(item) for item in candidates)
        best_candidates = [item for item in candidates if metadata_name_score(item) == best_score]
        if len(best_candidates) > 1:
            names = ", ".join(item.name for item in best_candidates)
            raise ValueError(f"Multiple WPS Working Sheet metadata files found in {path}: {names}")
        return best_candidates[0].resolve()
    return path.resolve()


def load_metadata_lookup(path: Path) -> dict[str, dict[str, str]]:
    metadata_path = resolve_metadata_path(path)
    lookup = {}
    barcode_order_rows: dict[tuple[str, str], int] = {}
    barcode_to_order: dict[str, str] = {}
    for row_number, row in enumerate(read_table_rows(metadata_path), start=1):
        meta = canonicalize_metadata_row(row)
        wps_sample_name = build_wps_sample_name(meta, row_number)
        if wps_sample_name:
            meta["sample_name"] = wps_sample_name
        barcode = meta.get("barcode")
        if barcode:
            order_number = meta.get("order_number", "")
            barcode_order = (order_number, barcode)
            if barcode_order in barcode_order_rows:
                raise ValueError(
                    f"Duplicate metadata order/barcode pair Order #{order_number}, {barcode}: "
                    f"rows {barcode_order_rows[barcode_order]} and {row_number}"
                )
            existing_order = barcode_to_order.get(barcode)
            if existing_order is not None and existing_order != order_number:
                raise ValueError(
                    f"Metadata barcode {barcode} appears under multiple orders: "
                    f"Order #{existing_order} and Order #{order_number}. "
                    "A barcode can only map to one sample/order in a single run."
                )
            barcode_order_rows[barcode_order] = row_number
            barcode_to_order[barcode] = order_number
            lookup[barcode] = meta
    return lookup


def group_packaged_by_order(packaged: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    grouped: dict[str, dict[str, object]] = {}
    for item in packaged:
        order_number = item["order_number"]
        order = grouped.setdefault(
            order_number,
            {
                "order_dir": item["order_dir"],
                "sample_count": 0,
                "samples": [],
                "reports": [],
            },
        )
        order["sample_count"] = int(order["sample_count"]) + 1
        order["samples"].append(item["sample_name"])
        order["reports"].append(item["pdf"])
    return grouped


def contig_label_from_filename(path: Path) -> str | None:
    match = re.search(r"contig[\s_-]*0*(\d+)", path.stem, flags=re.IGNORECASE)
    if not match:
        return None
    return f"contig{int(match.group(1)):03d}"


def record_id_for_barcode(barcode: str, contig_label: str | None = None) -> str:
    return f"{barcode}__{contig_label}" if contig_label else barcode


def sample_stem_for_record(barcode: str, metadata: dict[str, str], contig_label: str | None = None) -> str:
    base = slugify(metadata.get("sample_name") or barcode)
    return f"{base}_{contig_label}" if contig_label else base


def add_discovered_file(
    records: dict[str, dict[str, object]],
    discovery_errors: list[dict[str, str]],
    barcode: str,
    key: str,
    path: Path,
    contig_label: str | None = None,
) -> None:
    record_id = record_id_for_barcode(barcode, contig_label)
    records[record_id]["barcode"] = barcode
    if contig_label:
        records[record_id]["contig_label"] = contig_label
    existing = records[record_id].get(key)
    if existing is not None:
        discovery_errors.append(
            {
                "barcode": barcode,
                "record_id": record_id,
                "reason": f"multiple {key} files found for {record_id}: {existing}, {path}",
            }
        )
    else:
        records[record_id][key] = path


def contig_labels_for_barcode(records: dict[str, dict[str, object]], barcode: str) -> list[str]:
    labels = {
        record.get("contig_label")
        for record in records.values()
        if record.get("barcode") == barcode and isinstance(record.get("contig_label"), str)
    }
    return sorted(label for label in labels if isinstance(label, str))


def add_grouped_discovered_files(
    records: dict[str, dict[str, object]],
    discovery_errors: list[dict[str, str]],
    barcode: str,
    key: str,
    paths: list[Path],
) -> None:
    by_label: dict[str | None, list[Path]] = defaultdict(list)
    for path in sorted(paths):
        by_label[contig_label_from_filename(path)].append(path)

    for contig_label, labelled_paths in sorted(
        ((label, group) for label, group in by_label.items() if label is not None),
        key=lambda item: item[0] or "",
    ):
        for path in labelled_paths:
            add_discovered_file(records, discovery_errors, barcode, key, path, contig_label)

    unlabelled = by_label.get(None, [])
    if not unlabelled:
        return
    if len(unlabelled) == 1:
        add_discovered_file(records, discovery_errors, barcode, key, unlabelled[0])
        return

    existing_labels = contig_labels_for_barcode(records, barcode)
    if key == "fasta" and not existing_labels:
        print(
            f"detected multiple FASTA files for {barcode}; treating this as a likely mixed-contig sample "
            f"with {len(unlabelled)} contigs"
        )
        for index, path in enumerate(unlabelled, start=1):
            add_discovered_file(records, discovery_errors, barcode, key, path, f"contig{index:03d}")
        return

    if existing_labels and len(existing_labels) == len(unlabelled):
        print(
            f"detected multiple {key} files for {barcode}; pairing them with inferred mixed-contig labels "
            f"{', '.join(existing_labels)}"
        )
        for contig_label, path in zip(existing_labels, unlabelled):
            add_discovered_file(records, discovery_errors, barcode, key, path, contig_label)
        return

    first, *rest = unlabelled
    add_discovered_file(records, discovery_errors, barcode, key, first)
    for path in rest:
        add_discovered_file(records, discovery_errors, barcode, key, path)


def barcode_from_filename(path: Path) -> str | None:
    matches = {
        normalize_barcode(match)
        for match in re.findall(r"barcode\s*#?\s*0*\d+|barcode\d+", path.name, flags=re.IGNORECASE)
    }
    matches = {barcode for barcode in matches if barcode is not None}
    if len(matches) == 1:
        return next(iter(matches))
    return None


def discover_files_for_key(
    records: dict[str, dict[str, object]],
    discovery_errors: list[dict[str, str]],
    key: str,
    directory: Path,
    matcher,
) -> None:
    if not directory.exists():
        raise ValueError(f"{key} directory does not exist: {directory}")
    if not directory.is_dir():
        raise ValueError(f"{key} path is not a directory: {directory}")
    grouped: dict[str, list[Path]] = defaultdict(list)
    for path in directory.iterdir():
        if not path.is_file():
            continue
        barcode = matcher(path)
        if barcode is None:
            continue
        grouped[barcode].append(path)
    for barcode, paths in sorted(grouped.items()):
        add_grouped_discovered_files(records, discovery_errors, barcode, key, paths)


def match_barcode_file(path: Path, extensions: set[str], required_terms: set[str] | None = None) -> str | None:
    if path.suffix.lower() not in extensions:
        return None
    name = normalize_header(path.stem)
    if required_terms and not all(term in name for term in required_terms):
        return None
    return barcode_from_filename(path)


def is_fastq_path(path: Path) -> bool:
    name = path.name.lower()
    return name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))


def infer_fastq_barcode(path: Path) -> str | None:
    barcode = barcode_from_filename(path)
    if barcode is not None:
        return barcode
    for parent in path.parents:
        if re.fullmatch(r"barcode\d+", parent.name, re.IGNORECASE):
            return normalize_barcode(parent.name)
    return None


def scan_fastq_until(path: Path, min_records: int, min_bases: int) -> tuple[int, int]:
    count = 0
    bases = 0
    with open_fastq_text(path) as handle:
        while count < RAW_FASTQ_SCAN_RECORD_LIMIT and (count < min_records or bases < min_bases):
            header = handle.readline().rstrip("\n\r")
            if not header:
                break
            seq = handle.readline().rstrip("\n\r")
            plus = handle.readline().rstrip("\n\r")
            qual = handle.readline().rstrip("\n\r")
            if not header.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
                raise ValueError(f"Malformed FASTQ record in {path}")
            count += 1
            bases += len(seq)
    return count, bases


def raw_fastq_rejection_reason(path: Path) -> str | None:
    try:
        record_count, base_count = scan_fastq_until(path, RAW_FASTQ_MIN_RECORDS, RAW_FASTQ_MIN_BASES)
    except ValueError as exc:
        return str(exc)
    if record_count < RAW_FASTQ_MIN_RECORDS:
        return f"{path} has only {record_count} FASTQ records; expected at least {RAW_FASTQ_MIN_RECORDS}"
    if base_count < RAW_FASTQ_MIN_BASES:
        return f"{path} has only {base_count:,} read bases; expected at least {RAW_FASTQ_MIN_BASES:,}"
    return None


def should_skip_discovered_input(path: Path, resolved_excludes: Iterable[Path]) -> bool:
    resolved_path = path.resolve()
    if any(path_is_inside(resolved_path, excluded) for excluded in resolved_excludes):
        return True
    return any(parent.name == "_work" or parent.name.startswith("WPS Data_Order #") for parent in path.parents)


def infer_bam_barcode(path: Path) -> tuple[str | None, str | None]:
    filename_barcodes = sorted(
        {
            normalize_barcode(match)
            for match in re.findall(r"barcode\d+", path.name, flags=re.IGNORECASE)
        }
    )
    filename_barcodes = [barcode for barcode in filename_barcodes if barcode is not None]
    folder_barcodes = [
        normalize_barcode(parent.name)
        for parent in path.parents
        if re.fullmatch(r"barcode\d+", parent.name, re.IGNORECASE)
    ]
    folder_barcode = folder_barcodes[0] if folder_barcodes else None

    if len(filename_barcodes) > 1:
        if folder_barcode and folder_barcode in filename_barcodes:
            return folder_barcode, None
        return None, f"BAM filename contains multiple barcode tokens {filename_barcodes}: {path}"

    filename_barcode = filename_barcodes[0] if filename_barcodes else None
    if filename_barcode and folder_barcode and filename_barcode != folder_barcode:
        return None, f"BAM filename barcode {filename_barcode} does not match folder barcode {folder_barcode}: {path}"
    return folder_barcode or filename_barcode, None


def discover_bam_files(
    records: dict[str, dict[str, object]],
    discovery_errors: list[dict[str, str]],
    directory: Path,
    exclude_dirs: Iterable[Path] | None = None,
) -> None:
    if not directory.exists():
        raise ValueError(f"bam directory does not exist: {directory}")
    if not directory.is_dir():
        raise ValueError(f"bam path is not a directory: {directory}")

    resolved_excludes = [
        path.resolve()
        for path in (exclude_dirs or [])
        if path.resolve() != directory.resolve()
    ]

    grouped: dict[str, list[Path]] = defaultdict(list)
    for path in directory.rglob("*.bam"):
        if not path.is_file():
            continue
        if should_skip_discovered_input(path, resolved_excludes):
            continue
        barcode, error = infer_bam_barcode(path)
        if error:
            discovery_errors.append({"barcode": "unknown", "reason": error})
            continue
        if barcode is None:
            continue
        grouped[barcode].append(path)
    for barcode, paths in sorted(grouped.items()):
        add_grouped_discovered_files(records, discovery_errors, barcode, "bam", paths)


def discover_raw_fastq_files(
    records: dict[str, dict[str, object]],
    directory: Path,
    exclude_dirs: Iterable[Path] | None = None,
) -> None:
    if not directory.exists():
        raise ValueError(f"raw FASTQ directory does not exist: {directory}")
    if not directory.is_dir():
        raise ValueError(f"raw FASTQ path is not a directory: {directory}")

    resolved_excludes = [
        path.resolve()
        for path in (exclude_dirs or [])
        if path.resolve() != directory.resolve()
    ]

    grouped: dict[str, list[Path]] = defaultdict(list)
    for path in directory.rglob("*"):
        if not path.is_file() or not is_fastq_path(path):
            continue
        if should_skip_discovered_input(path, resolved_excludes):
            continue
        barcode = infer_fastq_barcode(path)
        if barcode is None:
            continue
        grouped[barcode].append(path)

    for barcode, paths in sorted(grouped.items()):
        scanned = {}
        rejected = {}
        for path in sorted(paths):
            try:
                record_count, base_count = scan_fastq_until(path, RAW_FASTQ_MIN_RECORDS, RAW_FASTQ_MIN_BASES)
            except ValueError as exc:
                rejected[path] = str(exc)
                continue
            scanned[path] = (record_count, base_count)
            if record_count < RAW_FASTQ_MIN_RECORDS:
                rejected[path] = (
                    f"{path} has only {record_count} FASTQ records; expected at least {RAW_FASTQ_MIN_RECORDS}"
                )
            elif base_count < RAW_FASTQ_MIN_BASES:
                rejected[path] = (
                    f"{path} has only {base_count:,} read bases; expected at least {RAW_FASTQ_MIN_BASES:,}"
                )
        eligible_paths = [path for path in sorted(paths) if path not in rejected]
        if not eligible_paths:
            record = records[barcode]
            record["barcode"] = barcode
            details = "; ".join(rejected[path] for path in sorted(rejected))
            record["raw_fastq_error"] = (
                f"no appropriately sized raw FASTQ found for {barcode}; {details}"
                if details
                else f"no appropriately sized raw FASTQ found for {barcode}"
            )
            continue

        selected = max(
            eligible_paths,
            key=lambda path: (
                scanned[path][1],
                scanned[path][0],
                "final" not in normalize_header(path.stem),
                path.stat().st_size,
            ),
        )
        record = records[barcode]
        record["barcode"] = barcode
        record["raw_fastq"] = selected
        if len(paths) > 1:
            ignored = [str(path) for path in sorted(paths) if path != selected]
            record.setdefault("input_warnings", []).append(
                "multiple raw FASTQ files detected; using largest file "
                f"{selected} and ignoring {', '.join(ignored)}"
            )


def expand_shared_barcode_files(records: dict[str, dict[str, object]]) -> dict[str, dict[str, object]]:
    by_barcode: dict[str, list[str]] = defaultdict(list)
    for record_id, record in records.items():
        barcode = record.get("barcode")
        if isinstance(barcode, str):
            by_barcode[barcode].append(record_id)

    expanded = dict(records)
    for barcode, record_ids in by_barcode.items():
        contig_record_ids = [rid for rid in record_ids if expanded[rid].get("contig_label")]
        if not contig_record_ids:
            continue
        shared_record = expanded.get(barcode)
        if not shared_record:
            continue
        for key in ("bam", "fastq", "raw_fastq", "maf"):
            shared_value = shared_record.get(key)
            if shared_value is None:
                continue
            for record_id in contig_record_ids:
                expanded[record_id].setdefault(key, shared_value)
        shared_warnings = shared_record.get("input_warnings")
        if isinstance(shared_warnings, list):
            for record_id in contig_record_ids:
                expanded[record_id].setdefault("input_warnings", []).extend(shared_warnings)
        if not any(key in shared_record for key in ("fasta", "gbk")):
            expanded.pop(barcode, None)
    return expanded


def primary_record_id_for_mixed_contigs(records: dict[str, dict[str, object]], record_ids: list[str]) -> str:
    contig001 = [
        record_id
        for record_id in record_ids
        if records[record_id].get("contig_label") == "contig001"
    ]
    if contig001:
        return sorted(contig001)[0]

    def fasta_length(record_id: str) -> int:
        fasta_path = records[record_id].get("fasta")
        if not isinstance(fasta_path, Path):
            return 0
        try:
            return int(read_first_fasta_record(fasta_path)["length_bp"])
        except Exception:
            return 0

    return max(sorted(record_ids), key=fasta_length)


def collapse_mixed_contigs_to_primary(records: dict[str, dict[str, object]]) -> dict[str, dict[str, object]]:
    by_barcode: dict[str, list[str]] = defaultdict(list)
    for record_id, record in records.items():
        barcode = record.get("barcode")
        if isinstance(barcode, str):
            by_barcode[barcode].append(record_id)

    collapsed: dict[str, dict[str, object]] = {}
    for barcode, record_ids in sorted(by_barcode.items()):
        contig_record_ids = [
            record_id
            for record_id in record_ids
            if records[record_id].get("contig_label") and records[record_id].get("fasta")
        ]
        if len(contig_record_ids) <= 1:
            for record_id in record_ids:
                collapsed[record_id] = records[record_id]
            continue

        primary_record_id = primary_record_id_for_mixed_contigs(records, contig_record_ids)
        primary = dict(records[primary_record_id])
        shared_record = records.get(barcode, {})
        for key in ("bam", "fastq", "raw_fastq", "maf", "raw_fastq_error"):
            if key in shared_record:
                primary.setdefault(key, shared_record[key])
        shared_warnings = shared_record.get("input_warnings")
        if isinstance(shared_warnings, list):
            primary.setdefault("input_warnings", []).extend(shared_warnings)
        ignored_labels = [
            str(records[record_id].get("contig_label"))
            for record_id in sorted(contig_record_ids)
            if record_id != primary_record_id
        ]
        primary["mixed_contig_count"] = len(contig_record_ids)
        primary["primary_contig_label"] = primary.get("contig_label")
        primary["ignored_contig_labels"] = ignored_labels
        primary.pop("contig_label", None)
        print(
            f"detected {len(contig_record_ids)} contigs for {barcode}; "
            f"using {primary['primary_contig_label']} as the primary contig for report generation"
        )
        collapsed[barcode] = primary

    return collapsed


def discover_input_records(
    fasta_dir: Path,
    genbank_dir: Path,
    bam_dir: Path,
    fastq_dir: Path | None = None,
    maf_dir: Path | None = None,
    exclude_dirs: Iterable[Path] | None = None,
) -> tuple[dict[str, dict[str, object]], list[dict[str, str]]]:
    records: dict[str, dict[str, object]] = defaultdict(dict)
    discovery_errors: list[dict[str, str]] = []
    discover_files_for_key(
        records,
        discovery_errors,
        "fasta",
        fasta_dir,
        lambda path: match_barcode_file(path, {".fasta", ".fa"}, {"final"}),
    )
    discover_files_for_key(
        records,
        discovery_errors,
        "gbk",
        genbank_dir,
        lambda path: match_barcode_file(path, {".gbk"}, {"annotations"}),
    )
    discover_bam_files(records, discovery_errors, bam_dir, exclude_dirs=exclude_dirs)
    if fastq_dir is not None:
        discover_files_for_key(
            records,
            discovery_errors,
            "fastq",
            fastq_dir,
            lambda path: match_barcode_file(path, {".fastq", ".fq"}, {"final"}),
        )
        discover_raw_fastq_files(records, fastq_dir, exclude_dirs=exclude_dirs)
    if maf_dir is not None:
        discover_files_for_key(
            records,
            discovery_errors,
            "maf",
            maf_dir,
            lambda path: match_barcode_file(path, {".maf"}, {"assembly"}),
        )
    return collapse_mixed_contigs_to_primary(expand_shared_barcode_files(dict(records))), discovery_errors


def sample_stem_for_barcode(barcode: str, metadata: dict[str, str]) -> str:
    return slugify(metadata.get("sample_name") or barcode)


def find_sample_stem_collisions(
    records: dict[str, dict[str, object]],
    metadata_lookup: dict[str, dict[str, str]],
    requested: set[str] | None,
    invalid_record_ids: set[str],
) -> dict[str, str]:
    stems: dict[str, list[str]] = defaultdict(list)
    for record_id, record in records.items():
        barcode = record.get("barcode")
        if not isinstance(barcode, str):
            continue
        if requested is not None and barcode not in requested:
            continue
        if record_id in invalid_record_ids:
            continue
        if metadata_lookup and barcode not in metadata_lookup:
            continue
        contig_label = record.get("contig_label")
        stem = sample_stem_for_record(
            barcode,
            metadata_lookup.get(barcode, {}),
            contig_label if isinstance(contig_label, str) else None,
        )
        stems[stem].append(record_id)
    collisions = {}
    for stem, record_ids in stems.items():
        if len(record_ids) > 1:
            joined = ", ".join(sorted(record_ids))
            for record_id in record_ids:
                collisions[record_id] = f"sample name collision after slugify: {stem!r} used by {joined}"
    return collisions


def write_renamed_fasta(src_fasta: Path, dst_fasta: Path, sequence_name: str) -> dict[str, str | int]:
    record = read_first_fasta_record(src_fasta)
    sequence = "".join(ch if ord(ch) < 128 else "N" for ch in record["sequence"])
    dst_fasta.parent.mkdir(parents=True, exist_ok=True)
    with dst_fasta.open("w", encoding="ascii") as handle:
        handle.write(f">{sequence_name}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start : start + 80] + "\n")
    return {"name": sequence_name, "sequence": sequence, "length_bp": len(sequence)}


def rewrite_genbank_locus(src_gbk: Path, dst_gbk: Path, locus_name: str) -> None:
    dst_gbk.parent.mkdir(parents=True, exist_ok=True)
    lines = src_gbk.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    safe_locus = re.sub(r"[^A-Za-z0-9_.-]+", "_", locus_name).strip("_")[:16] or "record"
    out_lines = []
    for idx, line in enumerate(lines):
        if idx == 0 and line.startswith("LOCUS"):
            parts = line.split()
            if len(parts) >= 3:
                bp = parts[2]
                remainder = ""
                if "bp" in line:
                    remainder = line[line.index("bp") + 2 :]
                line = f"LOCUS       {safe_locus:<16}{bp:>11} bp{remainder}"
            else:
                line = f"LOCUS       {safe_locus}"
        out_lines.append(line)
    dst_gbk.write_text("\n".join(out_lines) + "\n", encoding="utf-8")


def parse_fastq_record(path: Path) -> tuple[str, str, str]:
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        header = handle.readline().rstrip("\n\r")
        seq = handle.readline().rstrip("\n\r")
        handle.readline()
        qual = handle.readline().rstrip("\n\r")
    if not header.startswith("@") or len(seq) != len(qual):
        raise ValueError(f"Could not parse FASTQ record from {path}")
    return header[1:], seq.upper(), qual


def open_fastq_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", errors="replace")
    return path.open("r", encoding="ascii", errors="replace")


def iter_fastq_read_lengths(path: Path):
    with open_fastq_text(path) as handle:
        while True:
            header = handle.readline().rstrip("\n\r")
            if not header:
                break
            seq = handle.readline().rstrip("\n\r")
            plus = handle.readline().rstrip("\n\r")
            qual = handle.readline().rstrip("\n\r")
            if not header.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
                raise ValueError(f"Malformed FASTQ record in {path}")
            yield header[1:].split()[0], len(seq)


def generate_ab1_files(
    fasta_path: Path,
    fastq_path: Path | None,
    output_dir: Path,
    output_stem: str,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    record = read_first_fasta_record(fasta_path)
    if fastq_path and fastq_path.exists():
        _, fastq_seq, fastq_qual = parse_fastq_record(fastq_path)
        if fastq_seq == record["sequence"]:
            q_scores = phred_from_ascii(fastq_qual, phred_offset=33)
        else:
            q_scores = None
    else:
        q_scores = None

    if q_scores is None:
        q_scores = [30] * len(record["sequence"])

    traces, peak_locs, seq, q_scores = synthesize_chromatogram(
        record["sequence"],
        q_scores,
        samples_per_base=12,
    )

    for stale_name in (
        f"{output_stem}_trace.ab1",
        f"{output_stem}_trace1.ab1",
        f"{output_stem}_trace2.ab1",
    ):
        stale_path = output_dir / stale_name
        if stale_path.exists():
            stale_path.unlink()

    out_ab1 = output_dir / f"{output_stem}_trace.ab1"
    write_real_ab1(out_ab1, traces, peak_locs, seq, q_scores)
    return [out_ab1]


def validate_expected_fastq(record: dict[str, Path]) -> list[str]:
    if "fastq" not in record:
        return ["missing expected FASTQ; AB1 generated from FASTA with default Q30 quality scores"]
    return []


def read_lengths_from_bam(bam_path: Path) -> list[int]:
    lengths = []
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            lengths.append(read.query_length or 0)
    return lengths


def read_lengths_by_name_from_bam(bam_path: Path) -> dict[str, int]:
    lengths = {}
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            read_name = read.query_name or ""
            if read_name in lengths:
                raise ValueError(f"Duplicate primary read name in raw BAM: {read_name!r}")
            lengths[read_name] = read.query_length or 0
    return lengths


def read_lengths_by_name_from_fastq(fastq_path: Path) -> dict[str, int]:
    lengths = {}
    for read_name, read_length in iter_fastq_read_lengths(fastq_path):
        if read_name in lengths:
            raise ValueError(f"Duplicate FASTQ read name in {fastq_path}: {read_name!r}")
        lengths[read_name] = read_length
    return lengths


def plot_pdf_coverage_map(per_base_csv: Path, low_conf_csv: Path, out_path: Path) -> Path:
    positions = []
    depths = []
    depth_by_pos = {}
    with per_base_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            pos = int(row["pos"])
            depth = int(row["depth"])
            positions.append(pos)
            depths.append(depth)
            depth_by_pos[pos] = depth

    low_positions = []
    with low_conf_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            low_positions.append(int(row["pos"]))
    low_depths = [depth_by_pos[pos] for pos in low_positions if pos in depth_by_pos]

    fig, ax = plt.subplots(figsize=(7.6, 3.2))
    if positions:
        ax.plot(positions, depths, color="#4f9da6", linewidth=1.3)
        if low_positions:
            ax.scatter(low_positions, low_depths, marker="x", color="#e67e22", s=16, linewidths=0.8)
        ax.set_xlim(left=0, right=max(positions))
        ax.set_ylim(bottom=0, top=y_axis_top_with_headroom(max(depths, default=0)))
    else:
        ax.text(0.5, 0.5, "No coverage data available", ha="center", va="center", transform=ax.transAxes)
        ax.set_xlim(left=0, right=1)
        ax.set_ylim(bottom=0, top=1)
    ax.set_xlabel("Base Position")
    ax.set_ylabel("Depth")
    ax.grid(alpha=0.12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def plot_read_length_vs_bases(
    raw_reads_path: Path,
    aligned_bam: Path,
    contig_length: int,
    out_path: Path,
    raw_reads_format: str = "bam",
) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mapped_names = set()
    with pysam.AlignmentFile(aligned_bam, "rb") as bam:
        seen_aligned_names = set()
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            read_name = read.query_name or ""
            if read_name in seen_aligned_names:
                raise ValueError(f"Duplicate primary read name in aligned BAM: {read_name!r}")
            seen_aligned_names.add(read_name)
            if read.is_unmapped:
                continue
            mapped_names.add(read_name)

    if raw_reads_format == "bam":
        raw_lengths_by_name = read_lengths_by_name_from_bam(raw_reads_path)
    elif raw_reads_format == "fastq":
        raw_lengths_by_name = read_lengths_by_name_from_fastq(raw_reads_path)
    else:
        raise ValueError(f"Unsupported raw read format: {raw_reads_format}")

    mapped_lengths = []
    other_lengths = []
    for read_name, qlen in raw_lengths_by_name.items():
        if qlen <= READ_LENGTH_DISTRIBUTION_MIN_DISPLAY_BP:
            continue
        if read_name in mapped_names:
            mapped_lengths.append(qlen)
        else:
            other_lengths.append(qlen)

    all_lengths = mapped_lengths + other_lengths
    fig, ax = plt.subplots(figsize=(7.6, 3.7))
    if all_lengths:
        max_len = max(all_lengths)
        bin_size = max(250, int(math.ceil(max_len / 24 / 50.0)) * 50)
        start = int(math.floor(min(all_lengths) / bin_size) * bin_size)
        stop = int(math.ceil(max_len / bin_size) * bin_size) + bin_size
        bins = np.arange(start, stop + bin_size, bin_size)
        centers = bins[:-1] + bin_size / 2
        mapped_bases_kb, _ = np.histogram(
            mapped_lengths,
            bins=bins,
            weights=np.array(mapped_lengths, dtype=float) / 1000.0,
        )
        other_bases_kb, _ = np.histogram(
            other_lengths,
            bins=bins,
            weights=np.array(other_lengths, dtype=float) / 1000.0,
        )

        band_left = contig_length * (1.0 - MULTIMER_TOLERANCE_FRACTION)
        band_right = contig_length * (1.0 + MULTIMER_TOLERANCE_FRACTION)
        ymax = max((mapped_bases_kb + other_bases_kb).max(), 1)
        ytop = y_axis_top_with_headroom(ymax)
        label_y = ymax + (ytop - ymax) * 0.5
        ax.axvspan(band_left, band_right, color="#d9d9d9", alpha=0.6, lw=0)
        ax.text(
            (band_left + band_right) / 2,
            label_y,
            "Monomer",
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
            zorder=4,
        )

        width = bin_size * 0.78
        ax.bar(centers, mapped_bases_kb, width=width, color=THEME["purple"], label="Mapped reads")
        ax.bar(centers, other_bases_kb, width=width, bottom=mapped_bases_kb, color=THEME["green"], label="Unmapped reads")
        ax.set_xlim(left=max(0, start - bin_size * 0.25), right=stop)
        ax.set_ylim(bottom=0, top=ytop)
    else:
        ax.text(0.5, 0.5, "No read-length data available", ha="center", va="center", transform=ax.transAxes)
        ax.set_xlim(left=0, right=1)
        ax.set_ylim(bottom=0, top=1)
    ax.set_xlabel("Read Length (bp)")
    ax.set_ylabel("Total Bases (kb)")
    ax.legend(loc="upper right", frameon=True, facecolor="white")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def find_existing_path(paths: Iterable[Path | None]) -> Path | None:
    for path in paths:
        if path is not None and path.exists():
            return path
    return None


def validate_length_consistency(metadata: dict[str, str], fasta_length_bp: int, report_summary: dict) -> list[str]:
    warnings = []
    expected_size_bp = parse_expected_size_bp(metadata.get("size"))
    if expected_size_bp is not None and not lengths_match(fasta_length_bp, expected_size_bp):
        warnings.append(
            "metadata Size does not match FASTA contig length; using FASTA length: "
            f"metadata={expected_size_bp:,} bp, FASTA={fasta_length_bp:,} bp"
        )

    genbank_summary = report_summary.get("genbank_summary") or {}
    genbank_length_bp = genbank_summary.get("length_bp")
    if genbank_length_bp is not None and int(genbank_length_bp) != int(fasta_length_bp):
        raise ValueError(
            "GenBank LOCUS length does not match FASTA contig length: "
            f"GenBank={int(genbank_length_bp):,} bp, FASTA={fasta_length_bp:,} bp"
        )
    return warnings


def add_logo(fig, logos: list[Path]) -> None:
    for logo in logos[:1]:
        if not logo.exists():
            continue
        ax = fig.add_axes([0.055, 0.918, 0.145, 0.075])
        ax.set_axis_off()
        ax.imshow(mpimg.imread(logo))


def draw_report_title(fig, logos: list[Path]) -> None:
    add_logo(fig, logos)
    header_ax = fig.add_axes([0.05, 0.91, 0.90, 0.018])
    header_ax.set_axis_off()
    header_ax.add_patch(Rectangle((0, 0.45), 1, 0.12, color=THEME["rule"], lw=0))

    ax = fig.add_axes([0.06, 0.855, 0.88, 0.048])
    ax.set_axis_off()
    ax.text(
        0.5,
        0.5,
        "Whole Plasmid Sequencing Report",
        ha="center",
        va="center",
        fontsize=20.5,
        fontweight="bold",
        color=THEME["title"],
    )


def draw_section_heading(fig, y: float, title: str) -> None:
    ax = fig.add_axes([0.06, y, 0.88, 0.04])
    ax.set_axis_off()
    ax.text(0.5, 0.5, title, ha="center", va="center", fontsize=13.5, fontweight="bold", color=THEME["heading"])


def rounded_table_path(fig, bbox, radius_y=0.34) -> MplPath:
    width_in = bbox[2] * fig.get_figwidth()
    height_in = bbox[3] * fig.get_figheight()
    radius_y = min(radius_y, 0.49)
    radius_x = min(radius_y * height_in / width_in, 0.49) if width_in else radius_y
    kappa = 0.5522847498
    vertices = [
        (radius_x, 0.0),
        (1.0 - radius_x, 0.0),
        (1.0 - radius_x + kappa * radius_x, 0.0),
        (1.0, radius_y - kappa * radius_y),
        (1.0, radius_y),
        (1.0, 1.0 - radius_y),
        (1.0, 1.0 - radius_y + kappa * radius_y),
        (1.0 - radius_x + kappa * radius_x, 1.0),
        (1.0 - radius_x, 1.0),
        (radius_x, 1.0),
        (radius_x - kappa * radius_x, 1.0),
        (0.0, 1.0 - radius_y + kappa * radius_y),
        (0.0, 1.0 - radius_y),
        (0.0, radius_y),
        (0.0, radius_y - kappa * radius_y),
        (radius_x - kappa * radius_x, 0.0),
        (radius_x, 0.0),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
    ]
    return MplPath(vertices, codes)


def draw_table(fig, bbox, headers, values, col_widths=None) -> None:
    ax = fig.add_axes(bbox)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    if col_widths is None:
        col_widths = [1 / len(headers)] * len(headers)
    width_total = sum(col_widths) or 1
    col_widths = [width / width_total for width in col_widths]

    rounded_border = PathPatch(
        rounded_table_path(fig, bbox),
        facecolor="white",
        edgecolor=THEME["table_border"],
        linewidth=2.2,
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.add_patch(rounded_border)

    grid_lines = []
    header_y = 0.52
    grid_lines.extend(ax.plot([0.0, 1.0], [header_y, header_y], color=THEME["table_grid"], linewidth=0.75))
    x = 0.0
    centers = []
    for width in col_widths:
        centers.append(x + width / 2)
        x += width
        if x < 0.999:
            grid_lines.extend(ax.plot([x, x], [0.0, 1.0], color=THEME["table_grid"], linewidth=0.75))
    for line in grid_lines:
        line.set_clip_path(rounded_border)

    def value_font_size(text: str, col_width: float) -> float:
        char_capacity = max(1.0, col_width * 88.0)
        if len(text) <= char_capacity:
            return 9.2
        return max(6.2, 9.2 * char_capacity / len(text))

    for x_center, width, header, value in zip(centers, col_widths, headers, values):
        ax.text(
            x_center,
            0.74,
            header,
            ha="center",
            va="center",
            fontsize=9.2,
            fontweight="bold",
            color="#333333",
            wrap=True,
        )
        ax.text(
            x_center,
            0.25,
            value,
            ha="center",
            va="center",
            fontsize=value_font_size(str(value), width),
            color="#333333",
            wrap=True,
        )


def draw_image(fig, bounds, image_path: Path) -> None:
    ax = fig.add_axes(bounds)
    ax.set_axis_off()
    ax.imshow(mpimg.imread(image_path))


def draw_footer(fig) -> None:
    footer_ax = fig.add_axes([0.05, 0.026, 0.90, 0.045])
    footer_ax.set_axis_off()
    footer_ax.add_patch(Rectangle((0, 0.86), 1, 0.035, color=THEME["rule"], lw=0))
    footer_ax.text(
        0.5,
        0.38,
        "Alta Biotech, LLC \u2022 2115 N Scranton St, 3040B, Aurora CO 80045 \u2022 720-640-9400 \u2022 Support@altabiotech.com \u2022 www.altabiotech.com",
        ha="center",
        va="center",
        fontsize=8.0,
        color="black",
    )


def report_date_value(metadata: dict[str, str]) -> str:
    return (
        metadata.get("report_date")
        or metadata.get("run_date")
        or metadata.get("order_date")
        or metadata.get("received_date")
        or ""
    )


def ecoli_contamination_pct(report_summary: dict) -> float | None:
    contamination = report_summary.get("contamination") or {}
    if contamination.get("ecoli_genomic_contamination_pct") is not None:
        return contamination["ecoli_genomic_contamination_pct"]
    return 0.0


def compute_host_dna_from_alignment(host_bam: Path) -> dict[str, object]:
    total_read_count = 0
    total_read_bases = 0
    host_read_count = 0
    host_read_bases = 0
    host_aligned_bases = 0
    primary_names = set()

    with pysam.AlignmentFile(host_bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            read_name = read.query_name or ""
            if read_name in primary_names:
                raise ValueError(f"Duplicate primary read name in host-aligned BAM: {read_name!r}")
            primary_names.add(read_name)
            read_length = read.query_length or 0
            total_read_count += 1
            total_read_bases += read_length
            if read_length <= 0 or read.is_unmapped:
                continue
            aligned_bp = read.query_alignment_length or 0
            aligned_pct = aligned_bp / read_length * 100.0
            if aligned_bp > HOST_DNA_MIN_ALIGNED_BP and aligned_pct > HOST_DNA_MIN_ALIGNED_PCT:
                host_read_count += 1
                host_read_bases += read_length
                host_aligned_bases += aligned_bp

    host_read_pct = host_read_count / total_read_count * 100.0 if total_read_count else 0.0
    host_read_base_pct = host_read_bases / total_read_bases * 100.0 if total_read_bases else 0.0
    host_aligned_base_pct = host_aligned_bases / total_read_bases * 100.0 if total_read_bases else 0.0
    return {
        "method": "ecoli_reference_alignment",
        "reference_fasta": str(DEFAULT_ECOLI_REFERENCE_FASTA),
        "host_aligned_bp_threshold": HOST_DNA_MIN_ALIGNED_BP,
        "host_aligned_pct_threshold": HOST_DNA_MIN_ALIGNED_PCT,
        "host_dna_pct": round(host_read_pct, 3),
        "host_read_pct": round(host_read_pct, 3),
        "host_read_base_pct": round(host_read_base_pct, 3),
        "host_aligned_base_pct": round(host_aligned_base_pct, 3),
        "host_read_count": host_read_count,
        "host_read_bases": host_read_bases,
        "host_aligned_bases": host_aligned_bases,
        "total_read_count": total_read_count,
        "total_read_bases": total_read_bases,
    }


def format_percent(value: float | None) -> str:
    if value is None:
        return "N/A"
    return f"{value:.2f}%"


def format_yes_no(value: bool | None) -> str:
    if value is None:
        return "N/A"
    return "Yes" if value else "No"


def validate_percent_value(name: str, value: object, allow_none: bool = False) -> float | None:
    if value is None and allow_none:
        return None
    if not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} is not a finite percentage: {value!r}")
    value = float(value)
    if not 0.0 <= value <= 100.0:
        raise ValueError(f"{name} is outside 0..100%: {value}")
    return value


def multimer_pdf_table(assembly: dict) -> tuple[list[str], list[float | None], list[float]]:
    required = [
        ("monomer_pct", "Monomer"),
        ("dimer_pct", "Dimer"),
        ("trimer_pct", "Trimer"),
        ("tetramer_pct", "Tetramer"),
    ]
    if assembly.get("multimer_denominator") == MULTIMER_DENOMINATOR_ALL_ELIGIBLE_READS:
        required.append(("unclassified_multimer_read_pct", "Unclassified"))

    missing = [key for key, _label in required if key not in assembly]
    if missing:
        raise ValueError(f"Missing multimer fields: {', '.join(missing)}")

    if not assembly.get("multimer_calculated"):
        return [label for _key, label in required], [None for _key, _label in required], [1 / len(required)] * len(required)

    values = [validate_percent_value(label, assembly[key]) for key, label in required]
    total = sum(value for value in values if value is not None)
    if not 99.0 <= total <= 101.0:
        raise ValueError(f"Multimer percentages should sum to ~100%, got {total:.2f}%")
    return [label for _key, label in required], values, [1 / len(required)] * len(required)


def render_pdf_report(
    report_summary: dict,
    metadata: dict[str, str],
    output_pdf: Path,
    sample_stem: str,
    coverage_png: Path,
    read_length_bases_png: Path,
    feature_map_png: Path | None,
    logos: list[Path],
) -> None:
    del feature_map_png
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    contig = report_summary["contig"]
    assembly = report_summary["assembly_status"]
    coverage = report_summary["coverage"]
    contamination = ecoli_contamination_pct(report_summary)
    multimer_headers, multimer_values, multimer_widths = multimer_pdf_table(assembly)
    with PdfPages(output_pdf) as pdf:
        fig = plt.figure(figsize=(8.5, 11))
        fig.patch.set_facecolor("white")
        draw_report_title(fig, logos)

        draw_section_heading(fig, 0.785, "Sample & Order Information")
        draw_table(
            fig,
            [0.035, 0.69, 0.93, 0.085],
            ["Sample Name", "Lot Number", "Order Number", "Report Date"],
            [
                metadata.get("sample_name") or sample_stem,
                metadata.get("sample_id", ""),
                metadata.get("order_number", "UNKNOWN"),
                report_date_value(metadata),
            ],
            col_widths=[0.52, 0.12, 0.18, 0.18],
        )

        draw_section_heading(fig, 0.615, "Assembly Summary")
        draw_table(
            fig,
            [0.04, 0.52, 0.92, 0.08],
            ["Contig Length (bp)", "Bases Mapped", "Reads Mapped", "Host DNA %", "Is Circular?"],
            [
                f"{contig['length_bp']:,}",
                f"{assembly.get('bases_mapped', 0):,} ({assembly.get('bases_mapped_pct', 0):.2f}%)",
                f"{assembly.get('reads_mapped', 0):,} ({assembly.get('reads_mapped_pct', 0):.2f}%)",
                format_percent(contamination),
                format_yes_no(contig.get("is_circular")),
            ],
        )

        draw_section_heading(fig, 0.445, "Nanopore Performance")
        draw_table(
            fig,
            [0.035, 0.35, 0.93, 0.08],
            ["Mean Read Depth", "Min/Max Depth", "Coverage", "Low Confidence Bases", "Single Contig?"],
            [
                f"{round(coverage.get('mean_depth', 0)):,}",
                f"{coverage.get('min_depth', 0):,} / {coverage.get('max_depth', 0):,}",
                f"{round(coverage.get('mean_depth', 0)):,}x",
                f"{coverage.get('low_confidence_count', 0):,}",
                format_yes_no(assembly.get("single_contig")),
            ],
            col_widths=[0.20, 0.20, 0.20, 0.23, 0.17],
        )

        draw_section_heading(fig, 0.275, "Multimer Analysis")
        draw_table(
            fig,
            [0.035, 0.18, 0.93, 0.08],
            multimer_headers,
            [format_percent(value) for value in multimer_values],
            col_widths=multimer_widths,
        )
        draw_footer(fig)

        pdf.savefig(fig)
        plt.close(fig)

        fig = plt.figure(figsize=(8.5, 11))
        fig.patch.set_facecolor("white")
        draw_report_title(fig, logos)

        draw_section_heading(fig, 0.785, "Coverage & Distribution Maps")

        cov_head = fig.add_axes([0.055, 0.705, 0.84, 0.04])
        cov_head.set_axis_off()
        cov_head.text(0.0, 0.5, "Coverage Map", ha="left", va="center", fontsize=12.5, fontweight="bold", color="#335F91")
        draw_image(fig, [0.18, 0.465, 0.64, 0.23], coverage_png)
        cov_note = fig.add_axes([0.18, 0.43, 0.64, 0.03])
        cov_note.set_axis_off()
        cov_note.text(
            0.5,
            0.5,
            'low confidence positions are marked with orange "X"',
            ha="center",
            va="center",
            fontsize=10,
            color="#666666",
            style="italic",
        )

        dist_head = fig.add_axes([0.055, 0.345, 0.84, 0.04])
        dist_head.set_axis_off()
        dist_head.text(0.0, 0.5, "Read Length Distribution", ha="left", va="center", fontsize=12.5, fontweight="bold", color="#335F91")
        draw_image(fig, [0.18, 0.075, 0.64, 0.27], read_length_bases_png)
        draw_footer(fig)

        pdf.savefig(fig)
        plt.close(fig)


def collect_logos(paths: Iterable[str]) -> list[Path]:
    logos = []
    for item in paths:
        path = Path(item)
        if path.exists():
            logos.append(path)
    if not logos and DEFAULT_LOGO_PATH.exists():
        logos.append(DEFAULT_LOGO_PATH)
    return logos


def resolve_input_dir_from_folder_name(folder_name: str) -> Path:
    text = folder_name.strip()
    folder = Path(text)
    if not text or folder.is_absolute() or len(folder.parts) != 1 or "\\" in text:
        raise ValueError(
            f"--folder-name must be one folder name under {INPUT_ROOT}, not a full path"
        )
    if folder.parts[0] in {".", ".."}:
        raise ValueError(f"--folder-name must name a real folder under {INPUT_ROOT}")
    return (INPUT_ROOT / folder.parts[0]).resolve()


def path_is_inside(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    return True


def remove_output_path(path: Path, output_root: Path) -> None:
    if not path_is_inside(path, output_root) or not path.exists():
        return
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


def remove_empty_package_dirs(order_dir: Path) -> None:
    for subdir in PACKAGE_SUBDIRS.values():
        path = order_dir / subdir
        if path.exists() and path.is_dir():
            try:
                path.rmdir()
            except OSError:
                pass
    try:
        order_dir.rmdir()
    except OSError:
        pass


def cleanup_previous_sample_output(output_root: Path, barcode: str, sample_stem: str) -> None:
    work_root = output_root / "_work"
    touched_order_dirs: set[Path] = set()
    if work_root.exists():
        for summary_path in work_root.glob("*/package_summary.json"):
            try:
                summary = json.loads(summary_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            if summary.get("barcode") != barcode and summary.get("sample_name") != sample_stem:
                continue
            paths = summary.get("paths") or {}
            order_dir = paths.get("order_dir")
            if order_dir:
                touched_order_dirs.add(Path(order_dir))
            for value in paths.values():
                if isinstance(value, list):
                    for item in value:
                        remove_output_path(Path(item), output_root)
                elif value and Path(value).name != "":
                    path = Path(value)
                    if path.is_file() or path.suffix:
                        remove_output_path(path, output_root)
            remove_output_path(summary_path.parent, output_root)

    current_work_dir = work_root / sample_stem
    current_summary = current_work_dir / "package_summary.json"
    if current_summary.exists():
        try:
            summary = json.loads(current_summary.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            summary = {}
        existing_barcode = summary.get("barcode")
        if existing_barcode and existing_barcode != barcode:
            raise ValueError(
                f"existing work directory {current_work_dir} belongs to {existing_barcode}, not {barcode}"
            )
    remove_output_path(current_work_dir, output_root)

    for order_dir in touched_order_dirs:
        remove_empty_package_dirs(order_dir)


def package_sample(
    record_id: str,
    record: dict[str, object],
    metadata: dict[str, str],
    output_root: Path,
    logos: list[Path],
    threads: int = 1,
    sort_memory: str = "768M",
    keep_intermediates: bool = False,
    allow_aligned_input: bool = False,
    use_bam: bool = False,
    multimer_denominator: str = DEFAULT_MULTIMER_DENOMINATOR,
) -> dict[str, object]:
    barcode = record.get("barcode")
    if not isinstance(barcode, str):
        raise ValueError(f"record {record_id} is missing barcode")
    contig_label_value = record.get("contig_label")
    contig_label = contig_label_value if isinstance(contig_label_value, str) else None
    sample_name = metadata.get("sample_name") or barcode
    sample_stem = sample_stem_for_record(barcode, metadata, contig_label)
    output_stem = sample_stem if contig_label else f"{sample_stem}_contig"
    sequence_name = output_stem
    order_number = metadata.get("order_number")
    if not order_number:
        raise ValueError("metadata is missing order_number; refusing to package under WPS Data_Order #UNKNOWN")
    fasta_path = record["fasta"]
    gbk_path = record["gbk"]
    bam_path = record.get("bam")
    raw_fastq_path = record.get("raw_fastq")
    fastq_path = record.get("fastq")
    maf_path = record.get("maf")
    if not isinstance(fasta_path, Path) or not isinstance(gbk_path, Path):
        raise ValueError(f"record {record_id} has invalid required file paths")
    if bam_path is not None and not isinstance(bam_path, Path):
        raise ValueError(f"record {record_id} has invalid BAM path")
    if raw_fastq_path is not None and not isinstance(raw_fastq_path, Path):
        raise ValueError(f"record {record_id} has invalid raw FASTQ path")
    if use_bam and bam_path is None:
        raise ValueError(f"record {record_id} needs a BAM file when --use-bam is set")
    if not use_bam and raw_fastq_path is None:
        detail = record.get("raw_fastq_error") or "no appropriately sized raw FASTQ was detected"
        raise ValueError(f"record {record_id} needs a raw FASTQ unless --use-bam is set: {detail}")
    if fastq_path is not None and not isinstance(fastq_path, Path):
        raise ValueError(f"record {record_id} has invalid FASTQ path")
    if maf_path is not None and not isinstance(maf_path, Path):
        raise ValueError(f"record {record_id} has invalid MAF path")

    fasta_record_count = count_fasta_records(fasta_path)
    if fasta_record_count != 1:
        raise ValueError(f"expected exactly one FASTA record for {barcode}, found {fasta_record_count}")

    order_dir = output_root / f"WPS Data_Order #{order_number}"
    cleanup_previous_sample_output(output_root, barcode, sample_stem)
    package_dirs = {name: order_dir / subdir for name, subdir in PACKAGE_SUBDIRS.items()}
    for path in package_dirs.values():
        path.mkdir(parents=True, exist_ok=True)

    work_dir = output_root / "_work" / sample_stem
    work_dir.mkdir(parents=True, exist_ok=True)

    fasta_out = package_dirs["fasta"] / f"{output_stem}.fa"
    gbk_out = package_dirs["gbk"] / f"{output_stem}.gbk"
    renamed = write_renamed_fasta(fasta_path, fasta_out, sequence_name)
    rewrite_genbank_locus(gbk_path, gbk_out, sequence_name)

    alignment_dir = work_dir / "alignment"
    if raw_fastq_path is not None and not use_bam:
        alignment_result = run_fastq_pipeline(
            raw_fastq_path,
            fasta_out,
            alignment_dir,
            minimap2_preset="map-ont",
            threads=threads,
            sort_memory=sort_memory,
            keep_intermediates=keep_intermediates,
        )
        raw_reads_path = raw_fastq_path
        raw_reads_format = "fastq"
    else:
        assert isinstance(bam_path, Path)
        alignment_result = run_pipeline(
            bam_path,
            fasta_out,
            alignment_dir,
            minimap2_preset="map-ont",
            threads=threads,
            sort_memory=sort_memory,
            keep_intermediates=keep_intermediates,
            allow_aligned_input=allow_aligned_input,
        )
        raw_reads_path = bam_path
        raw_reads_format = "bam"
    aligned_bam = Path(alignment_result["sorted_bam"])
    fasta_index = Path(f"{fasta_out}.fai")
    if fasta_index.exists():
        fasta_index.unlink()

    warnings = validate_expected_fastq(record)
    input_warnings = record.get("input_warnings")
    if isinstance(input_warnings, list):
        warnings.extend(str(warning) for warning in input_warnings)
    if raw_fastq_path is not None and not use_bam:
        warnings.append(f"using raw FASTQ for alignment and read metrics: {raw_fastq_path}")
    if raw_fastq_path is not None and use_bam:
        warnings.append(f"raw FASTQ detected but --use-bam requested; using BAM for alignment and read metrics: {bam_path}")
    if record.get("mixed_contig_count"):
        warnings.append(
            "multiple contigs detected; report generated from primary contig "
            f"{record.get('primary_contig_label', 'unknown')}"
        )
    host_contamination_details = None
    host_contamination_pct = None
    host_aligned_bam = None
    if DEFAULT_ECOLI_REFERENCE_FASTA.exists():
        if raw_fastq_path is not None and not use_bam:
            host_alignment_result = run_fastq_pipeline(
                raw_fastq_path,
                DEFAULT_ECOLI_REFERENCE_FASTA,
                work_dir / "host_alignment",
                minimap2_preset="map-ont",
                threads=threads,
                sort_memory=sort_memory,
                keep_intermediates=keep_intermediates,
            )
        else:
            assert isinstance(bam_path, Path)
            host_alignment_result = run_pipeline(
                bam_path,
                DEFAULT_ECOLI_REFERENCE_FASTA,
                work_dir / "host_alignment",
                minimap2_preset="map-ont",
                threads=threads,
                sort_memory=sort_memory,
                keep_intermediates=keep_intermediates,
                allow_aligned_input=allow_aligned_input,
            )
        host_bam = Path(host_alignment_result["sorted_bam"])
        host_aligned_bam = str(host_bam)
        host_contamination_details = compute_host_dna_from_alignment(host_bam)
        host_contamination_pct = float(host_contamination_details["host_dna_pct"])
    else:
        warnings.append(f"missing E. coli host genome reference: {DEFAULT_ECOLI_REFERENCE_FASTA}")

    report_dir = work_dir / "report"
    report_summary = generate_report_data(
        aligned_bam=aligned_bam,
        contig_fasta=fasta_out,
        out_dir=report_dir,
        reference_fasta=fasta_out,
        maf_path=maf_path,
        gbk_path=gbk_out,
        sample_name=sample_stem,
        low_confidence_qscore=LOW_CONFIDENCE_QSCORE,
        ecoli_contamination_pct=host_contamination_pct,
        ecoli_contamination_details=host_contamination_details,
        multimer_denominator=multimer_denominator,
    )
    warnings.extend(validate_length_consistency(metadata, renamed["length_bp"], report_summary))

    per_base_src = Path(report_summary["outputs"]["per_base_details_csv"])
    low_conf_src = Path(report_summary["outputs"]["low_confidence_bases_csv"])
    coverage_png = plot_pdf_coverage_map(per_base_src, low_conf_src, work_dir / "coverage_map_pdf.png")
    feature_map_value = report_summary["outputs"].get("feature_map_png")
    feature_map_png = find_existing_path([Path(feature_map_value)]) if feature_map_value else None

    per_base_dst = package_dirs["per_base"] / f"{output_stem}_per_base_details.csv"
    low_conf_dst = package_dirs["per_base"] / f"{output_stem}_low_confidence_bases.csv"
    shutil.copyfile(per_base_src, per_base_dst)
    shutil.copyfile(low_conf_src, low_conf_dst)

    ab1_paths = generate_ab1_files(fasta_out, fastq_path, package_dirs["ab1"], output_stem)

    bases_plot = plot_read_length_vs_bases(
        raw_reads_path,
        aligned_bam,
        renamed["length_bp"],
        work_dir / "read_length_vs_bases.png",
        raw_reads_format=raw_reads_format,
    )

    pdf_out = package_dirs["qc"] / f"{sample_stem}_report.pdf"
    render_pdf_report(
        report_summary=report_summary,
        metadata={**metadata, "barcode": barcode},
        output_pdf=pdf_out,
        sample_stem=sample_stem,
        coverage_png=coverage_png,
        read_length_bases_png=bases_plot,
        feature_map_png=feature_map_png,
        logos=logos,
    )

    summary_out = work_dir / "package_summary.json"
    summary_out.write_text(
        json.dumps(
            {
                "barcode": barcode,
                "record_id": record_id,
                "contig_label": contig_label,
                "mixed_contig_count": record.get("mixed_contig_count"),
                "primary_contig_label": record.get("primary_contig_label"),
                "ignored_contig_labels": record.get("ignored_contig_labels"),
                "sample_name": sample_stem,
                "order_number": order_number,
                "multimer_denominator": multimer_denominator,
                "warnings": warnings,
                "paths": {
                    "order_dir": str(order_dir),
                    "pdf": str(pdf_out),
                    "alignment_input": str(raw_reads_path),
                    "alignment_input_type": raw_reads_format,
                    "fasta": str(fasta_out),
                    "gbk": str(gbk_out),
                    "ab1": [str(path) for path in ab1_paths],
                    "per_base_details": str(per_base_dst),
                    "low_confidence": str(low_conf_dst),
                    "aligned_bam": str(aligned_bam),
                    "host_aligned_bam": host_aligned_bam,
                },
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    return {
        "barcode": barcode,
        "record_id": record_id,
        "contig_label": contig_label,
        "mixed_contig_count": record.get("mixed_contig_count"),
        "primary_contig_label": record.get("primary_contig_label"),
        "ignored_contig_labels": record.get("ignored_contig_labels"),
        "sample_name": sample_stem,
        "order_number": order_number,
        "order_dir": str(order_dir),
        "pdf": str(pdf_out),
        "warnings": warnings,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--folder-name",
        required=True,
        help=(
            f"Name of the folder under {INPUT_ROOT} containing all run inputs: "
            "barcodeXX.final.fasta/fa, barcodeXX.annotations.gbk, raw/unmapped BAMs, "
            "barcodeXX.final.fastq/fq, optional MAF files, and exactly one WPS Working Sheet metadata CSV, TSV, or XLSX file. "
            "Mixed-contig samples may use barcodeXX.contig001.final.fasta/fa and matching contig labels on related files. "
            "Missing consensus FASTQ files warn but do not stop packaging. "
            "BAMs may be directly inside it or inside barcodeXX subfolders. "
            "Raw FASTQ/FASTQ.GZ files with barcode tokens are required for alignment by default; pass --use-bam to force BAM input."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            f"Output directory for WPS order folders. "
            f"Default: <folder-name>/{DEFAULT_OUTPUT_SUBDIR}"
        ),
    )
    parser.add_argument(
        "--barcodes",
        nargs="*",
        default=None,
        help="Optional barcode filter, for example: barcode01 barcode02 or 1 2",
    )
    parser.add_argument(
        "--logo",
        action="append",
        default=[],
        help="Optional logo image(s) to place in the PDF header. Pass up to two times.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads for minimap2 and samtools alignment steps. Default: 1.",
    )
    parser.add_argument(
        "--sort-memory",
        default="768M",
        help="Memory per samtools sort thread, for example 768M or 2G. Default: 768M.",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep alignment intermediates reads.fastq, aligned.sam, and aligned.unsorted.bam.",
    )
    parser.add_argument(
        "--allow-aligned-input",
        action="store_true",
        help="Allow BAMs that already contain mapped primary reads.",
    )
    parser.add_argument(
        "--use-bam",
        action="store_true",
        help=(
            "Use BAM files for alignment, host DNA, read-length distribution, and multimer metrics even when "
            "larger raw FASTQ/FASTQ.GZ files are detected. By default, an appropriately sized barcoded FASTQ is required."
        ),
    )
    parser.add_argument(
        "--multimer-denominator",
        choices=MULTIMER_DENOMINATOR_CHOICES,
        default=DEFAULT_MULTIMER_DENOMINATOR,
        help=(
            "Denominator for monomer/dimer/trimer/tetramer percentages. "
            "classified-reads reports percentages only among reads classified as 1x-4x; "
            "all-eligible-reads includes unclassified eligible mapped reads in the denominator."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    folder_name = args.folder_name.strip()
    input_dir = resolve_input_dir_from_folder_name(folder_name)
    if not input_dir.exists():
        raise ValueError(f"input directory does not exist: {input_dir}")
    if not input_dir.is_dir():
        raise ValueError(f"input path is not a directory: {input_dir}")

    output_dir = Path(args.output_dir).resolve() if args.output_dir else input_dir / DEFAULT_OUTPUT_SUBDIR
    metadata_path = resolve_metadata_path(input_dir)
    metadata_lookup = load_metadata_lookup(metadata_path)
    records, discovery_errors = discover_input_records(
        fasta_dir=input_dir,
        genbank_dir=input_dir,
        bam_dir=input_dir,
        fastq_dir=input_dir,
        maf_dir=input_dir,
        exclude_dirs=[output_dir],
    )

    requested = None
    if args.barcodes:
        requested = {normalize_barcode(item) for item in args.barcodes}

    logos = collect_logos(args.logo)
    packaged = []
    skipped = list(discovery_errors)
    invalid_record_ids = {item.get("record_id", item["barcode"]) for item in discovery_errors}
    invalid_barcodes = {
        item["barcode"]
        for item in discovery_errors
        if item.get("record_id", item["barcode"]) == item["barcode"]
    }
    considered_record_ids = set(records)
    if requested is not None:
        considered_record_ids = {
            record_id
            for record_id in considered_record_ids
            if isinstance(records[record_id].get("barcode"), str) and records[record_id]["barcode"] in requested
        }

    sample_stem_collisions = find_sample_stem_collisions(records, metadata_lookup, requested, invalid_record_ids)
    for record_id, reason in sorted(sample_stem_collisions.items()):
        record = records.get(record_id, {})
        barcode = record.get("barcode", record_id)
        skipped.append({"barcode": str(barcode), "record_id": record_id, "reason": reason})
    invalid_record_ids.update(sample_stem_collisions)

    record_barcodes = {record["barcode"] for record in records.values() if isinstance(record.get("barcode"), str)}
    for record_id in sorted(considered_record_ids):
        record = records[record_id]
        barcode = record.get("barcode")
        if not isinstance(barcode, str):
            continue
        if record_id in invalid_record_ids or barcode in invalid_barcodes:
            continue
        if metadata_lookup and barcode not in metadata_lookup:
            skipped.append({"barcode": barcode, "record_id": record_id, "reason": "barcode not present in metadata"})

    for barcode in sorted(set(metadata_lookup) - record_barcodes):
        if requested is not None and barcode not in requested:
            continue
        skipped.append({"barcode": barcode, "reason": "metadata row has no matching EPI2ME files"})

    for record_id in sorted(records):
        record = records[record_id]
        barcode = record.get("barcode")
        if not isinstance(barcode, str):
            skipped.append({"barcode": record_id, "record_id": record_id, "reason": "record is missing barcode"})
            continue
        if requested is not None and barcode not in requested:
            continue
        if record_id in invalid_record_ids or barcode in invalid_barcodes:
            continue
        if metadata_lookup and barcode not in metadata_lookup:
            continue
        missing = [key for key in ("fasta", "gbk") if key not in record]
        if args.use_bam:
            if "bam" not in record:
                missing.append("bam")
        elif "raw_fastq" not in record:
            detail = record.get("raw_fastq_error") or "no appropriately sized raw FASTQ was detected"
            raise ValueError(f"{barcode}: {detail}. Pass --use-bam to force BAM input.")
        if missing:
            skipped.append({"barcode": barcode, "record_id": record_id, "reason": f"missing required files: {', '.join(missing)}"})
            continue
        metadata = metadata_lookup.get(barcode, {})
        try:
            packaged.append(
                package_sample(
                    record_id,
                    record,
                    metadata,
                    output_dir,
                    logos,
                    threads=args.threads,
                    sort_memory=args.sort_memory,
                    keep_intermediates=args.keep_intermediates,
                    allow_aligned_input=args.allow_aligned_input,
                    use_bam=args.use_bam,
                    multimer_denominator=args.multimer_denominator,
                )
            )
        except subprocess.CalledProcessError as exc:
            skipped.append({"barcode": barcode, "record_id": record_id, "reason": f"command failed ({exc.returncode}): {' '.join(map(str, exc.cmd))}"})
        except Exception as exc:
            skipped.append({"barcode": barcode, "record_id": record_id, "reason": str(exc)})

    grouped_orders = group_packaged_by_order(packaged)

    summary = {
        "folder_name": folder_name,
        "input_root": str(INPUT_ROOT),
        "input_dir": str(input_dir),
        "metadata": str(metadata_path) if metadata_path else None,
        "multimer_denominator": args.multimer_denominator,
        "output_dir": str(output_dir),
        "packaged": packaged,
        "orders": grouped_orders,
        "skipped": skipped,
    }
    summary_path = output_dir / "run_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    print(f"run_summary_json: {summary_path}")
    print(f"packaged_count: {len(packaged)}")
    print(f"order_count: {len(grouped_orders)}")
    print(f"skipped_count: {len(skipped)}")
    for order_number, order in sorted(grouped_orders.items()):
        print(f"order: {order_number} -> {order['order_dir']} ({order['sample_count']} sample(s))")
    for item in packaged:
        print(f"packaged: {item.get('record_id', item['barcode'])} -> {item['order_dir']}")
        for warning in item.get("warnings", []):
            print(f"warning: {item.get('record_id', item['barcode'])} ({warning})")
    for item in skipped:
        print(f"skipped: {item.get('record_id', item['barcode'])} ({item['reason']})")


if __name__ == "__main__":
    main()
