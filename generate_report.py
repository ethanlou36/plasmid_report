#!/usr/bin/env python3
"""
Generate a report-style summary from aligned plasmid sequencing outputs.

Inputs:
- aligned BAM
- contig FASTA
- optional reference FASTA
- optional MAF
- optional GenBank

Outputs:
- nested summary dictionary (returned from generate_report_data)
- JSON written to out_dir/report_summary.json
- figures written to out_dir
"""

import argparse
import csv
import json
import math
import os
import re
import statistics
import tempfile
from collections import Counter
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "altabiotech_mplcache"))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam
from Bio import SeqIO
from Bio.SeqIO.InsdcIO import _insdc_location_string

from bam_to_per_base_data import summarize_bam_to_table

MULTIMER_TOLERANCE_FRACTION = 0.15
MIN_MULTIMER_ALIGNMENT_FRACTION = 0.0
MIN_MULTIMER_MAPQ = 0
MULTIMER_ELIGIBILITY_RULE = (
    "primary mapped reads classified by full read length; no MAPQ or aligned-fraction cutoff"
)
READ_LENGTH_DISTRIBUTION_MIN_DISPLAY_BP = 1000
PLOT_Y_AXIS_HEADROOM_FRACTION = 0.10
NON_MULTIMER_PEAK_MIN_BASE_FRACTION = 0.08
NON_MULTIMER_PEAK_MIN_READ_COUNT = 3
MULTIMER_DENOMINATOR_CLASSIFIED_READS = "classified-reads"
MULTIMER_DENOMINATOR_ALL_ELIGIBLE_READS = "all-eligible-reads"
MULTIMER_DENOMINATOR_CHOICES = (
    MULTIMER_DENOMINATOR_CLASSIFIED_READS,
    MULTIMER_DENOMINATOR_ALL_ELIGIBLE_READS,
)
DEFAULT_MULTIMER_DENOMINATOR = MULTIMER_DENOMINATOR_CLASSIFIED_READS


def y_axis_top_with_headroom(max_value, headroom_fraction=PLOT_Y_AXIS_HEADROOM_FRACTION):
    if max_value <= 0:
        return 1
    return max_value * (1.0 + headroom_fraction)


def read_first_fasta_record(path):
    name = None
    chunks = []
    with open(path, "r", encoding="utf-8-sig", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].split()[0]
            else:
                chunks.append(line)
    if name is None:
        raise ValueError(f"No FASTA record found in {path}")
    sequence = "".join(chunks).upper()
    return {"name": name, "sequence": sequence, "length_bp": len(sequence)}


def count_fasta_records(path):
    count = 0
    with open(path, "r", encoding="utf-8-sig", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def parse_plasmidasaurus_summary(path):
    text = Path(path).read_text(encoding="ascii", errors="replace")
    moles = {}
    mass = {}

    header_match = re.search(r"^\s+(.*)$", text, flags=re.MULTILINE)
    if header_match:
        header = header_match.group(1)
        multiples = [int(match) for match in re.findall(r"(\d+)-mer", header)]
        moles_match = re.search(r"^moles\s+([0-9.\s]+)$", text, flags=re.MULTILINE)
        mass_match = re.search(r"^mass\s+([0-9.\s]+)$", text, flags=re.MULTILINE)
        if moles_match:
            values = [float(value) for value in moles_match.group(1).split()]
            for multiple, value in zip(multiples, values):
                moles[f"{multiple}-mer"] = value
        if mass_match:
            values = [float(value) for value in mass_match.group(1).split()]
            for multiple, value in zip(multiples, values):
                mass[f"{multiple}-mer"] = value

    contamination_match = re.search(
        r"E\. coli genomic contamination:\s*([0-9.]+)%",
        text,
        flags=re.IGNORECASE,
    )
    contamination_pct = float(contamination_match.group(1)) if contamination_match else None
    return {
        "multimer_by_moles_pct": moles or None,
        "multimer_by_mass_pct": mass or None,
        "ecoli_genomic_contamination_pct": contamination_pct,
    }


def _summary_qualifier_value(values):
    """Keep the prior scalar qualifier shape while retaining repeated values."""
    if not isinstance(values, (list, tuple)):
        return values
    if not values or values == [""]:
        return True
    if len(values) == 1:
        return values[0]
    return list(values)


def _summary_location_segments(location):
    """Return Biopython locations as one-based, end-inclusive segments."""
    if location is None:
        return []

    segments = []
    for part in location.parts:
        try:
            start = int(part.start) + 1
            end = int(part.end)
        except TypeError:
            start = None
            end = None
        segments.append(
            {
                "start": start,
                "end": end,
                "strand": part.strand,
            }
        )
    return segments


def parse_genbank_summary(path):
    with open(path, "r", encoding="utf-8-sig", errors="replace") as handle:
        record = SeqIO.read(handle, "genbank")

    features = []
    for feature in record.features:
        qualifiers = {
            key: _summary_qualifier_value(values)
            for key, values in feature.qualifiers.items()
        }
        location = (
            _insdc_location_string(feature.location, len(record))
            if feature.location is not None
            else ""
        )
        features.append(
            {
                "type": feature.type,
                "location": location,
                "qualifiers": qualifiers,
                "segments": _summary_location_segments(feature.location),
            }
        )

    feature_counts = Counter(feature["type"] for feature in features)
    labels = []
    for feature in features:
        label = feature["qualifiers"].get("label")
        if isinstance(label, list):
            labels.extend(value for value in label if value)
        elif label:
            labels.append(label)

    return {
        "locus_name": record.name,
        "length_bp": len(record),
        "is_circular": str(record.annotations.get("topology", "")).lower() == "circular",
        "feature_count": len(features),
        "feature_type_counts": dict(sorted(feature_counts.items())),
        "labels": labels,
        "features": features,
    }


def parse_maf_summary(path, contig_length):
    blocks = []
    current_score = None
    current_s_lines = []

    with open(path, "r", encoding="ascii", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                if len(current_s_lines) >= 2:
                    blocks.append({"score": current_score, "s_lines": current_s_lines[:2]})
                current_score = None
                current_s_lines = []
                continue
            if line.startswith("a "):
                match = re.search(r"score=([^\s]+)", line)
                current_score = float(match.group(1)) if match else None
            elif line.startswith("s "):
                current_s_lines.append(line.split())

    if len(current_s_lines) >= 2:
        blocks.append({"score": current_score, "s_lines": current_s_lines[:2]})

    summaries = []
    for block in blocks:
        first = block["s_lines"][0]
        second = block["s_lines"][1]
        seq1 = first[6]
        seq2 = second[6]
        matches = 0
        mismatches = 0
        gaps = 0
        for base1, base2 in zip(seq1, seq2):
            if base1 == "-" or base2 == "-":
                gaps += 1
            elif base1.upper() == base2.upper():
                matches += 1
            else:
                mismatches += 1

        aligned_columns = len(seq1)
        start1 = int(first[2])
        start2 = int(second[2])
        size1 = int(first[3])
        size2 = int(second[3])
        is_full_self = start1 == 0 and start2 == 0 and size1 == contig_length and size2 == contig_length
        summaries.append(
            {
                "score": block["score"],
                "start_1": start1,
                "start_2": start2,
                "aligned_columns": aligned_columns,
                "aligned_bases_1": size1,
                "aligned_bases_2": size2,
                "matches": matches,
                "mismatches": mismatches,
                "gaps": gaps,
                "identity_pct": (matches / aligned_columns * 100.0) if aligned_columns else 0.0,
                "is_full_self": is_full_self,
            }
        )

    repeat_blocks = [summary for summary in summaries if not summary["is_full_self"]]
    largest_repeat = max(repeat_blocks, key=lambda item: item["aligned_bases_1"], default=None)
    return {
        "block_count": len(summaries),
        "repeat_block_count": len(repeat_blocks),
        "largest_repeat_block_bp": largest_repeat["aligned_bases_1"] if largest_repeat else 0,
        "largest_repeat_identity_pct": round(largest_repeat["identity_pct"], 3) if largest_repeat else 0.0,
        "largest_repeat_span_fraction": (
            largest_repeat["aligned_bases_1"] / contig_length if largest_repeat and contig_length else 0.0
        ),
        "top_repeat_blocks": repeat_blocks[:10],
    }


def read_per_base_rows(csv_path):
    with open(csv_path, "r", encoding="ascii", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def compute_n50(lengths):
    if not lengths:
        return 0
    total = sum(lengths)
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= total / 2:
            return length
    return 0


def classify_multimer(read_length, contig_length, tolerance_fraction=MULTIMER_TOLERANCE_FRACTION, max_multiple=4):
    if contig_length <= 0:
        return None
    ratio = read_length / contig_length
    candidates = []
    for multiple in range(1, max_multiple + 1):
        delta = abs(ratio - multiple)
        if delta <= tolerance_fraction or math.isclose(
            delta,
            tolerance_fraction,
            rel_tol=1e-12,
            abs_tol=1e-12,
        ):
            candidates.append((delta, multiple))
    if not candidates:
        return None
    candidates.sort()
    return candidates[0][1]


def read_length_distribution_bins(read_lengths):
    if not read_lengths:
        return None
    max_len = max(read_lengths)
    bin_size = max(250, int(math.ceil(max_len / 24 / 50.0)) * 50)
    start = int(math.floor(min(read_lengths) / bin_size) * bin_size)
    stop = int(math.ceil(max_len / bin_size) * bin_size) + bin_size
    return bin_size, list(range(start, stop + bin_size, bin_size))


def read_length_peak_is_multimer(length_bp, contig_length):
    return classify_multimer(length_bp, contig_length=contig_length) is not None


def detect_non_multimer_read_length_peaks(read_lengths, contig_length):
    displayed_lengths = [
        length
        for length in read_lengths
        if length > READ_LENGTH_DISTRIBUTION_MIN_DISPLAY_BP
    ]
    bins_info = read_length_distribution_bins(displayed_lengths)
    if contig_length <= 0 or bins_info is None:
        return {
            "single_contig_by_read_lengths": True,
            "non_multimer_peak_count": 0,
            "non_multimer_peak_read_count": 0,
            "non_multimer_peak_base_count": 0,
            "non_multimer_peak_base_pct": 0.0,
            "non_multimer_peak_intervals": [],
        }

    _bin_size, bins = bins_info
    bin_count = len(bins) - 1
    base_totals = [0] * bin_count
    read_counts = [0] * bin_count
    total_bases = sum(displayed_lengths)
    for length in displayed_lengths:
        index = min(max(int((length - bins[0]) // (bins[1] - bins[0])), 0), bin_count - 1)
        base_totals[index] += length
        read_counts[index] += 1

    intervals = []
    for index, base_total in enumerate(base_totals):
        if read_counts[index] < NON_MULTIMER_PEAK_MIN_READ_COUNT:
            continue
        base_fraction = base_total / total_bases if total_bases else 0.0
        if base_fraction < NON_MULTIMER_PEAK_MIN_BASE_FRACTION:
            continue
        previous_bases = base_totals[index - 1] if index > 0 else 0
        next_bases = base_totals[index + 1] if index + 1 < bin_count else 0
        if base_total < previous_bases or base_total < next_bases:
            continue
        center = (bins[index] + bins[index + 1]) / 2.0
        if read_length_peak_is_multimer(center, contig_length):
            continue
        intervals.append(
            {
                "start_bp": bins[index],
                "end_bp": bins[index + 1],
                "center_bp": round(center, 1),
                "read_count": read_counts[index],
                "base_count": base_total,
                "base_pct": round(base_fraction * 100.0, 3),
            }
        )

    peak_read_count = sum(item["read_count"] for item in intervals)
    peak_base_count = sum(item["base_count"] for item in intervals)
    return {
        "single_contig_by_read_lengths": len(intervals) == 0,
        "non_multimer_peak_count": len(intervals),
        "non_multimer_peak_read_count": peak_read_count,
        "non_multimer_peak_base_count": peak_base_count,
        "non_multimer_peak_base_pct": round((peak_base_count / total_bases * 100.0), 3) if total_bases else 0.0,
        "non_multimer_peak_intervals": intervals,
    }


def length_in_intervals(length, intervals):
    return any(interval["start_bp"] <= length < interval["end_bp"] for interval in intervals)


def multimer_breakdown(
    read_lengths,
    contig_length,
    tolerance_fraction=MULTIMER_TOLERANCE_FRACTION,
    max_multiple=4,
):
    counts = {multiple: 0 for multiple in range(1, max_multiple + 1)}
    masses = {multiple: 0 for multiple in range(1, max_multiple + 1)}

    for read_length in read_lengths:
        multiple = classify_multimer(
            read_length,
            contig_length=contig_length,
            tolerance_fraction=tolerance_fraction,
            max_multiple=max_multiple,
        )
        if multiple is None:
            continue
        counts[multiple] += 1
        masses[multiple] += read_length

    total_read_count = len(read_lengths)
    total_base_count = sum(read_lengths)
    classified_read_count = sum(counts.values())
    classified_base_count = sum(masses.values())
    unclassified_read_count = total_read_count - classified_read_count
    unclassified_base_count = total_base_count - classified_base_count
    moles_pct = {
        f"{multiple}-mer": round((counts[multiple] / total_read_count * 100.0), 3)
        if total_read_count
        else None
        for multiple in counts
    }
    classified_moles_pct = {
        f"{multiple}-mer": round((counts[multiple] / classified_read_count * 100.0), 3)
        if classified_read_count
        else None
        for multiple in counts
    }
    mass_pct = {
        f"{multiple}-mer": round((masses[multiple] / total_base_count * 100.0), 3)
        if total_base_count
        else None
        for multiple in masses
    }
    classified_mass_pct = {
        f"{multiple}-mer": round((masses[multiple] / classified_base_count * 100.0), 3)
        if classified_base_count
        else None
        for multiple in masses
    }
    return {
        "counts": {f"{multiple}-mer": counts[multiple] for multiple in counts},
        "bases": {f"{multiple}-mer": masses[multiple] for multiple in masses},
        "moles_pct": moles_pct,
        "mass_pct": mass_pct,
        "classified_moles_pct": classified_moles_pct,
        "classified_mass_pct": classified_mass_pct,
        "unclassified_read_pct": round((unclassified_read_count / total_read_count * 100.0), 3)
        if total_read_count
        else None,
        "unclassified_base_pct": round((unclassified_base_count / total_base_count * 100.0), 3)
        if total_base_count
        else None,
        "eligible_read_count": total_read_count,
        "eligible_base_count": total_base_count,
        "classified_read_count": classified_read_count,
        "classified_base_count": classified_base_count,
        "unclassified_read_count": unclassified_read_count,
        "unclassified_base_count": unclassified_base_count,
        "calculated": total_read_count > 0,
    }


def is_multimer_eligible_alignment(read):
    return not read.is_unmapped


def selected_multimer_percentages(
    multimer: dict,
    multimer_denominator: str,
) -> tuple[dict[str, float | None], dict[str, float | None], bool]:
    if multimer_denominator == MULTIMER_DENOMINATOR_CLASSIFIED_READS:
        return multimer["classified_moles_pct"], multimer["classified_mass_pct"], multimer["classified_read_count"] > 0
    if multimer_denominator == MULTIMER_DENOMINATOR_ALL_ELIGIBLE_READS:
        return multimer["moles_pct"], multimer["mass_pct"], multimer["eligible_read_count"] > 0
    raise ValueError(
        f"Unsupported multimer denominator {multimer_denominator!r}; "
        f"expected one of {', '.join(MULTIMER_DENOMINATOR_CHOICES)}"
    )


def bam_summary(bam_path, contig_length, multimer_denominator=DEFAULT_MULTIMER_DENOMINATOR):
    total_records = 0
    total_read_bases = 0
    primary_read_lengths = []
    mapped_primary_read_lengths = []
    multimer_eligible_read_lengths = []
    mapped_primary_count = 0
    mapped_bases = 0
    primary_names = set()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        contig = bam.references[0] if bam.nreferences == 1 else None
        for read in bam.fetch(until_eof=True):
            total_records += 1
            if read.is_secondary or read.is_supplementary:
                continue
            read_name = read.query_name or ""
            if read_name in primary_names:
                raise ValueError(f"Duplicate primary read name in aligned BAM: {read_name!r}")
            primary_names.add(read_name)
            qlen = read.query_length or 0
            total_read_bases += qlen
            primary_read_lengths.append(qlen)
            if read.is_unmapped:
                continue
            mapped_primary_count += 1
            mapped_primary_read_lengths.append(qlen)
            mapped_bases += read.query_alignment_length or 0
            if is_multimer_eligible_alignment(read):
                multimer_eligible_read_lengths.append(qlen)

    read_length_contig_detection = detect_non_multimer_read_length_peaks(primary_read_lengths, contig_length)
    non_multimer_peak_intervals = read_length_contig_detection["non_multimer_peak_intervals"]
    adjusted_multimer_eligible_read_lengths = [
        length
        for length in multimer_eligible_read_lengths
        if not length_in_intervals(length, non_multimer_peak_intervals)
    ]
    multimer = multimer_breakdown(adjusted_multimer_eligible_read_lengths, contig_length)
    selected_moles_pct, selected_mass_pct, multimer_calculated = selected_multimer_percentages(
        multimer,
        multimer_denominator,
    )
    return {
        "total_records": total_records,
        "total_bases": total_read_bases,
        "primary_reads": len(primary_read_lengths),
        "mapped_primary_reads": mapped_primary_count,
        "mapped_read_pct": (mapped_primary_count / len(primary_read_lengths) * 100.0) if primary_read_lengths else 0.0,
        "mapped_bases": mapped_bases,
        "mapped_base_pct": (mapped_bases / total_read_bases * 100.0) if total_read_bases else 0.0,
        "mean_read_length": round(statistics.fmean(primary_read_lengths), 3) if primary_read_lengths else 0.0,
        "median_read_length": statistics.median(primary_read_lengths) if primary_read_lengths else 0,
        "read_length_n50": compute_n50(primary_read_lengths),
        "mapped_mean_read_length": round(statistics.fmean(mapped_primary_read_lengths), 3)
        if mapped_primary_read_lengths
        else 0.0,
        "monomer_pct": selected_mass_pct["1-mer"],
        "dimer_pct": selected_mass_pct["2-mer"],
        "trimer_pct": selected_mass_pct["3-mer"],
        "tetramer_pct": selected_mass_pct["4-mer"],
        "multimer_by_moles_pct": selected_moles_pct,
        "multimer_by_mass_pct": selected_mass_pct,
        "multimer_by_all_eligible_reads_pct": multimer["moles_pct"],
        "multimer_by_classified_reads_pct": multimer["classified_moles_pct"],
        "multimer_calculated": multimer_calculated,
        "multimer_denominator": multimer_denominator,
        "multimer_eligible_read_count": multimer["eligible_read_count"],
        "multimer_eligible_base_count": multimer["eligible_base_count"],
        "multimer_excluded_non_contig_peak_read_count": len(multimer_eligible_read_lengths)
        - len(adjusted_multimer_eligible_read_lengths),
        "multimer_excluded_non_contig_peak_base_count": sum(multimer_eligible_read_lengths)
        - sum(adjusted_multimer_eligible_read_lengths),
        "unclassified_multimer_read_pct": multimer["unclassified_read_pct"],
        "unclassified_multimer_base_pct": multimer["unclassified_base_pct"],
        "classified_multimer_read_count": multimer["classified_read_count"],
        "classified_multimer_base_count": multimer["classified_base_count"],
        "unclassified_multimer_read_count": multimer["unclassified_read_count"],
        "unclassified_multimer_base_count": multimer["unclassified_base_count"],
        "multimer_min_alignment_fraction": MIN_MULTIMER_ALIGNMENT_FRACTION,
        "multimer_min_mapq": MIN_MULTIMER_MAPQ,
        "multimer_eligibility_rule": MULTIMER_ELIGIBILITY_RULE,
        "read_length_contig_detection": read_length_contig_detection,
        "primary_read_lengths": primary_read_lengths,
    }


def coverage_summary(per_base_rows, low_conf_rows, contig_length):
    depths = [int(row["depth"]) for row in per_base_rows]
    low_positions = [int(row["pos"]) for row in low_conf_rows]
    return {
        "mean_depth": round(statistics.fmean(depths), 3) if depths else 0.0,
        "median_depth": statistics.median(depths) if depths else 0,
        "min_depth": min(depths) if depths else 0,
        "max_depth": max(depths) if depths else 0,
        "covered_bases": sum(1 for depth in depths if depth > 0),
        "coverage_breadth_pct": (sum(1 for depth in depths if depth > 0) / contig_length * 100.0)
        if contig_length
        else 0.0,
        "low_confidence_count": len(low_positions),
        "low_confidence_positions_preview": low_positions[:20],
    }


def plot_coverage_map(per_base_rows, low_conf_rows, out_path, title):
    positions = [int(row["pos"]) for row in per_base_rows]
    depths = [int(row["depth"]) for row in per_base_rows]
    depth_by_pos = {int(row["pos"]): int(row["depth"]) for row in per_base_rows}
    low_positions = [int(row["pos"]) for row in low_conf_rows]
    low_depths = [depth_by_pos[pos] for pos in low_positions if pos in depth_by_pos]

    fig, ax = plt.subplots(figsize=(11, 4.8))
    ax.plot(positions, depths, color="#2f6c9e", linewidth=1.1)
    if low_positions:
        ax.scatter(low_positions, low_depths, marker="x", color="#e67e22", s=30, linewidths=1.1)
    ax.set_xlim(left=0, right=max(positions))
    ax.set_ylim(bottom=0, top=y_axis_top_with_headroom(max(depths, default=0)))
    ax.margins(x=0)
    ax.set_xlabel("Base Position")
    ax.set_ylabel("Depth")
    ax.set_title(title)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_read_length_histogram(read_lengths, out_path, title):
    fig, ax = plt.subplots(figsize=(8.5, 5))
    read_lengths = [length for length in read_lengths if length > READ_LENGTH_DISTRIBUTION_MIN_DISPLAY_BP]
    bins = min(80, max(15, int(math.sqrt(len(read_lengths))))) if read_lengths else 20
    ax.hist(read_lengths, bins=bins, color="#6d8f72", edgecolor="white")
    ax.set_xlabel("Read Length (bp)")
    ax.set_ylabel("Read Count")
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def location_segments(location):
    return [(int(start), int(end)) for start, end in re.findall(r"(\d+)\.\.(\d+)", location)]


def feature_location_segments(feature):
    """Prefer parsed GenBank segments, with a fallback for legacy summaries."""
    segments = feature.get("segments")
    if segments is not None:
        return [
            (segment["start"], segment["end"])
            for segment in segments
            if segment.get("start") is not None and segment.get("end") is not None
        ]
    return location_segments(feature.get("location", ""))


def plot_feature_map(gbk_summary, contig_length, out_path, title):
    features = gbk_summary["features"]
    if not features:
        return None

    color_by_type = {
        "CDS": "#4c78a8",
        "gene": "#f58518",
        "promoter": "#54a24b",
        "rep_origin": "#e45756",
        "terminator": "#72b7b2",
        "misc_feature": "#b279a2",
        "ncRNA": "#ff9da6",
        "intron": "#9d755d",
    }

    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.hlines(0, 0, contig_length, color="black", linewidth=1)

    row_count = max(1, min(8, len(features)))
    for idx, feature in enumerate(features):
        y = 0.2 + (idx % row_count) * 0.18
        feature_type = feature["type"]
        color = color_by_type.get(feature_type, "#7f7f7f")
        segments = feature_location_segments(feature)
        for start, end in segments:
            ax.broken_barh([(start - 1, end - start + 1)], (y, 0.12), facecolors=color)
        label = feature["qualifiers"].get("label", feature_type)
        if isinstance(label, list):
            label = ", ".join(str(value) for value in label)
        if segments:
            ax.text(segments[0][0], y + 0.14, label, fontsize=7, va="bottom")

    ax.set_xlim(0, contig_length)
    ax.set_ylim(-0.05, 0.2 + row_count * 0.18 + 0.18)
    ax.set_xlabel("Base Position")
    ax.set_yticks([])
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.15)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def generate_report_data(
    aligned_bam,
    contig_fasta,
    out_dir,
    reference_fasta=None,
    maf_path=None,
    gbk_path=None,
    sample_name=None,
    low_confidence_qscore=12,
    plasmidasaurus_summary_txt=None,
    ecoli_contamination_pct=None,
    ecoli_contamination_details=None,
    multimer_denominator=DEFAULT_MULTIMER_DENOMINATOR,
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    contig = read_first_fasta_record(contig_fasta)
    reference = read_first_fasta_record(reference_fasta) if reference_fasta else None
    gbk_summary = parse_genbank_summary(gbk_path) if gbk_path else None
    maf_summary = parse_maf_summary(maf_path, contig["length_bp"]) if maf_path else None
    bam_stats = bam_summary(aligned_bam, contig["length_bp"], multimer_denominator=multimer_denominator)
    vendor_summary = parse_plasmidasaurus_summary(plasmidasaurus_summary_txt) if plasmidasaurus_summary_txt else None
    if ecoli_contamination_pct is None and vendor_summary is not None:
        ecoli_contamination_pct = vendor_summary["ecoli_genomic_contamination_pct"]
    if ecoli_contamination_pct is None:
        ecoli_contamination_pct = 0.0

    per_base_csv = out_dir / "per_base_details.csv"
    low_conf_csv = out_dir / "low_confidence_bases.csv"
    summarize_bam_to_table(
        bam_path=aligned_bam,
        output_csv=per_base_csv,
        reference_path=reference_fasta,
        include_zero_depth=bool(reference_fasta),
        low_confidence_out=low_conf_csv,
        low_confidence_qscore=low_confidence_qscore,
    )

    per_base_rows = read_per_base_rows(per_base_csv)
    low_conf_rows = read_per_base_rows(low_conf_csv)
    coverage_stats = coverage_summary(per_base_rows, low_conf_rows, contig["length_bp"])

    coverage_png = out_dir / "coverage_map.png"
    read_len_png = out_dir / "read_length_distribution.png"
    feature_map_png = out_dir / "feature_map.png"

    plot_coverage_map(
        per_base_rows,
        low_conf_rows,
        coverage_png,
        title=f"{contig['name']} Coverage Map",
    )
    plot_read_length_histogram(
        bam_stats["primary_read_lengths"],
        read_len_png,
        title="Read Length Distribution",
    )
    feature_map_written = None
    if gbk_summary is not None:
        feature_map_written = plot_feature_map(
            gbk_summary,
            contig["length_bp"],
            feature_map_png,
            title="Annotation Map",
        )

    report = {
        "sample_name": sample_name or contig["name"],
        "contig": {
            "name": contig["name"],
            "length_bp": contig["length_bp"],
            "is_circular": gbk_summary["is_circular"] if gbk_summary is not None else None,
            "fasta_record_count": count_fasta_records(contig_fasta),
        },
        "reference": (
            {
                "name": reference["name"],
                "length_bp": reference["length_bp"],
            }
            if reference is not None
            else None
        ),
        "sequencing_information": {
            "total_dna_reads": bam_stats["primary_reads"],
            "total_dna_bases": bam_stats["total_bases"],
            "mean_read_length": bam_stats["mean_read_length"],
            "median_read_length": bam_stats["median_read_length"],
            "read_length_n50": bam_stats["read_length_n50"],
        },
        "assembly_status": {
            "contig": contig["name"],
            "length_bp": contig["length_bp"],
            "reads_mapped": bam_stats["mapped_primary_reads"],
            "reads_mapped_pct": round(bam_stats["mapped_read_pct"], 3),
            "bases_mapped": bam_stats["mapped_bases"],
            "bases_mapped_pct": round(bam_stats["mapped_base_pct"], 3),
            "coverage_x": round(coverage_stats["mean_depth"], 3),
            "median_coverage_x": coverage_stats["median_depth"],
            "is_circular": gbk_summary["is_circular"] if gbk_summary is not None else None,
            "monomer_pct": bam_stats["monomer_pct"],
            "dimer_pct": bam_stats["dimer_pct"],
            "trimer_pct": bam_stats["trimer_pct"],
            "tetramer_pct": bam_stats["tetramer_pct"],
            "multimer_by_moles_pct": bam_stats["multimer_by_moles_pct"],
            "multimer_by_mass_pct": bam_stats["multimer_by_mass_pct"],
            "multimer_by_all_eligible_reads_pct": bam_stats["multimer_by_all_eligible_reads_pct"],
            "multimer_by_classified_reads_pct": bam_stats["multimer_by_classified_reads_pct"],
            "multimer_calculated": bam_stats["multimer_calculated"],
            "multimer_denominator": bam_stats["multimer_denominator"],
            "multimer_eligible_read_count": bam_stats["multimer_eligible_read_count"],
            "multimer_eligible_base_count": bam_stats["multimer_eligible_base_count"],
            "multimer_excluded_non_contig_peak_read_count": bam_stats["multimer_excluded_non_contig_peak_read_count"],
            "multimer_excluded_non_contig_peak_base_count": bam_stats["multimer_excluded_non_contig_peak_base_count"],
            "unclassified_multimer_read_pct": bam_stats["unclassified_multimer_read_pct"],
            "unclassified_multimer_base_pct": bam_stats["unclassified_multimer_base_pct"],
            "classified_multimer_read_count": bam_stats["classified_multimer_read_count"],
            "classified_multimer_base_count": bam_stats["classified_multimer_base_count"],
            "unclassified_multimer_read_count": bam_stats["unclassified_multimer_read_count"],
            "unclassified_multimer_base_count": bam_stats["unclassified_multimer_base_count"],
            "multimer_min_alignment_fraction": bam_stats["multimer_min_alignment_fraction"],
            "multimer_min_mapq": bam_stats["multimer_min_mapq"],
            "multimer_eligibility_rule": bam_stats["multimer_eligibility_rule"],
            "read_length_contig_detection": bam_stats["read_length_contig_detection"],
            "single_contig": count_fasta_records(contig_fasta) == 1
            and bam_stats["read_length_contig_detection"]["single_contig_by_read_lengths"],
        },
        "coverage": coverage_stats,
        "contamination": {
            "ecoli_genomic_contamination_pct": ecoli_contamination_pct,
            "ecoli_genomic_contamination_details": ecoli_contamination_details,
        },
        "maf_summary": maf_summary,
        "genbank_summary": (
            {
                "locus_name": gbk_summary["locus_name"],
                "length_bp": gbk_summary["length_bp"],
                "feature_count": gbk_summary["feature_count"],
                "feature_type_counts": gbk_summary["feature_type_counts"],
                "labels": gbk_summary["labels"],
            }
            if gbk_summary is not None
            else None
        ),
        "plasmidasaurus_summary": vendor_summary,
        "outputs": {
            "per_base_details_csv": str(per_base_csv),
            "low_confidence_bases_csv": str(low_conf_csv),
            "coverage_map_png": str(coverage_png),
            "read_length_distribution_png": str(read_len_png),
            "feature_map_png": str(feature_map_written) if feature_map_written else None,
        },
    }

    summary_json = out_dir / "report_summary.json"
    with open(summary_json, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
    report["outputs"]["report_summary_json"] = str(summary_json)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aligned-bam", required=True, help="Aligned BAM input")
    parser.add_argument("--contig-fasta", required=True, help="Consensus/contig FASTA")
    parser.add_argument("--reference-fasta", default=None, help="Optional reference FASTA")
    parser.add_argument("--maf", default=None, help="Optional MAF alignment file")
    parser.add_argument("--gbk", default=None, help="Optional annotated GenBank file")
    parser.add_argument(
        "--plasmidasaurus-summary-txt",
        default=None,
        help="Optional Plasmidasaurus summary TXT used to import E. coli contamination and vendor summary values",
    )
    parser.add_argument("--out-dir", required=True, help="Output directory for JSON and figures")
    parser.add_argument("--sample-name", default=None, help="Optional sample name override")
    parser.add_argument(
        "--ecoli-contamination-pct",
        type=float,
        default=None,
        help="Optional E. coli genomic contamination percentage to include in the report",
    )
    parser.add_argument(
        "--low-confidence-qscore",
        type=int,
        default=12,
        help="Mean BAM base-quality threshold for marking low-confidence positions",
    )
    parser.add_argument(
        "--multimer-denominator",
        choices=MULTIMER_DENOMINATOR_CHOICES,
        default=DEFAULT_MULTIMER_DENOMINATOR,
        help=(
            "Denominator for monomer/dimer/trimer/tetramer percentages. "
            "classified-reads reports percentages only among reads classified as 1x-4x; "
            "all-eligible-reads includes unclassified eligible mapped reads in the denominator. "
            "Reported multimer percentages are base-weighted."
        ),
    )
    args = parser.parse_args()

    report = generate_report_data(
        aligned_bam=args.aligned_bam,
        contig_fasta=args.contig_fasta,
        reference_fasta=args.reference_fasta,
        maf_path=args.maf,
        gbk_path=args.gbk,
        out_dir=args.out_dir,
        sample_name=args.sample_name,
        low_confidence_qscore=args.low_confidence_qscore,
        plasmidasaurus_summary_txt=args.plasmidasaurus_summary_txt,
        ecoli_contamination_pct=args.ecoli_contamination_pct,
        multimer_denominator=args.multimer_denominator,
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
