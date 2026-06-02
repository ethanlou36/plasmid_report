#!/usr/bin/env python3
"""
Align an unaligned BAM to a reference and produce SAM/BAM outputs.

Pipeline:
1. Verify the input BAM is not already aligned
2. Convert unaligned BAM to FASTQ with pysam
3. Index the reference FASTA with samtools
4. Align reads to the reference with minimap2
5. Convert SAM to BAM with samtools
6. Sort the BAM with samtools
7. Index the sorted BAM with samtools
8. Remove bulky intermediates unless requested
"""

import argparse
import gzip
import hashlib
from pathlib import Path
import shutil
import subprocess

import pysam


def require_tool(name):
    path = shutil.which(name)
    if path is None:
        raise RuntimeError(f"{name} is not installed or not on PATH")
    return path


def run_command(args, stdout=None):
    subprocess.run(args, check=True, stdout=stdout)


def sanitize_fastq_read_name(read_name):
    text = str(read_name or "read")
    cleaned = "".join(ch if 33 <= ord(ch) <= 126 and not ch.isspace() else "_" for ch in text)
    cleaned = cleaned.strip("@") or "read"
    if cleaned == text and len(cleaned) <= 240:
        return cleaned
    digest = hashlib.sha1(text.encode("utf-8", errors="replace")).hexdigest()[:10]
    return f"{cleaned[:229]}_{digest}"


def inspect_bam_input(bam_path):
    primary_reads = 0
    mapped_primary_reads = 0
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            primary_reads += 1
            if not read.is_unmapped:
                mapped_primary_reads += 1
    return {
        "primary_reads": primary_reads,
        "mapped_primary_reads": mapped_primary_reads,
    }


def open_fastq_text(path):
    path = Path(path)
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="ascii", errors="replace")
    return path.open("r", encoding="ascii", errors="replace")


def inspect_fastq_input(fastq_path):
    records = 0
    bases = 0
    with open_fastq_text(fastq_path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().rstrip("\n\r")
            plus = handle.readline()
            qual = handle.readline().rstrip("\n\r")
            if not plus or not header.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
                raise ValueError(f"Malformed FASTQ record in {fastq_path}")
            records += 1
            bases += len(seq)
    if records == 0:
        raise ValueError(f"No FASTQ records found in {fastq_path}")
    return {
        "primary_reads": records,
        "total_bases": bases,
    }


def bam_to_fastq(_samtools, bam_path, out_fastq):
    """
    Convert a BAM to FASTQ with pysam.

    `samtools fastq` produced empty output for the unmapped MinKNOW BAMs in this
    workflow even though the records carry sequence and quality. Writing FASTQ
    directly from pysam is more reliable here.
    """
    written = 0
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam, open(
        out_fastq, "w", encoding="ascii"
    ) as handle:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary:
                continue
            seq = read.query_sequence
            if not seq:
                continue
            quals = read.query_qualities
            qual_text = (
                "".join(chr(min(max(q, 0), 93) + 33) for q in quals)
                if quals is not None
                else "I" * len(seq)
            )
            handle.write(f"@{sanitize_fastq_read_name(read.query_name)}\n{seq}\n+\n{qual_text}\n")
            written += 1

    if written == 0:
        raise ValueError(f"No FASTQ records were written from {bam_path}")


def index_reference_fasta(samtools, reference_fasta):
    run_command([samtools, "faidx", str(reference_fasta)])


def count_alignment_records(samtools, bam_path, mapped_only=False, threads=1):
    command = [samtools, "view", "-@", str(threads), "-c"]
    if mapped_only:
        command.extend(["-F", "4"])
    command.append(str(bam_path))
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return int(result.stdout.strip())


def remove_intermediates(paths):
    removed = []
    for path in paths:
        path = Path(path)
        if path.exists():
            path.unlink()
            removed.append(str(path))
    return removed


def align_fastq_to_reference(
    reads_fastq,
    reference_path,
    out_dir,
    minimap2_preset="map-ont",
    threads=1,
    sort_memory="768M",
    keep_intermediates=False,
    removable_intermediates=(),
):
    minimap2 = require_tool("minimap2")
    samtools = require_tool("samtools")

    if threads < 1:
        raise ValueError("threads must be >= 1")
    if not str(sort_memory).strip():
        raise ValueError("sort_memory must be a non-empty samtools memory value, for example 768M or 2G")

    reference_path = Path(reference_path)
    if reference_path.suffix.lower() not in {".fa", ".fasta", ".fna"}:
        raise ValueError("Reference must be a FASTA file (.fa, .fasta, or .fna)")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sam_path = out_dir / "aligned.sam"
    unsorted_bam = out_dir / "aligned.unsorted.bam"
    sorted_bam = out_dir / "aligned.sorted.bam"

    index_reference_fasta(samtools, reference_path)

    with open(sam_path, "wb") as sam_handle:
        run_command(
            [
                minimap2,
                "-t",
                str(threads),
                "-ax",
                minimap2_preset,
                str(reference_path),
                str(reads_fastq),
            ],
            stdout=sam_handle,
        )

    run_command(
        [
            samtools,
            "view",
            "-@",
            str(threads),
            "-bS",
            "-o",
            str(unsorted_bam),
            str(sam_path),
        ]
    )
    run_command(
        [
            samtools,
            "sort",
            "-@",
            str(threads),
            "-m",
            str(sort_memory),
            "-o",
            str(sorted_bam),
            str(unsorted_bam),
        ]
    )
    run_command([samtools, "index", "-@", str(threads), str(sorted_bam)])

    aligned_reads = count_alignment_records(samtools, sorted_bam, mapped_only=True, threads=threads)
    sam_reads = count_alignment_records(samtools, sorted_bam, threads=threads)

    intermediates = [*removable_intermediates, sam_path, unsorted_bam]
    removed_intermediates = [] if keep_intermediates else remove_intermediates(intermediates)

    return {
        "reads_fastq": str(reads_fastq) if reads_fastq.exists() else None,
        "reference_fasta": str(reference_path),
        "reference_fai": str(reference_path) + ".fai",
        "sam": str(sam_path) if sam_path.exists() else None,
        "unsorted_bam": str(unsorted_bam) if unsorted_bam.exists() else None,
        "sorted_bam": str(sorted_bam),
        "sorted_bai": str(sorted_bam) + ".bai",
        "aligned_reads": aligned_reads,
        "sam_reads": sam_reads,
        "aligner": "minimap2",
        "threads": threads,
        "sort_memory": str(sort_memory),
        "intermediates_kept": keep_intermediates,
        "intermediates_removed": removed_intermediates,
    }


def run_pipeline(
    bam_path,
    reference_path,
    out_dir,
    minimap2_preset="map-ont",
    threads=1,
    sort_memory="768M",
    keep_intermediates=False,
    allow_aligned_input=False,
):
    bam_stats = inspect_bam_input(bam_path)
    if bam_stats["mapped_primary_reads"] and not allow_aligned_input:
        raise ValueError(
            f"Input BAM appears to contain {bam_stats['mapped_primary_reads']} already-mapped primary reads. "
            "This pipeline expects raw/unaligned MinKNOW BAM input. Pass allow_aligned_input=True "
            "or --allow-aligned-input only if you intentionally want to realign these reads."
        )

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    reads_fastq = out_dir / "reads.fastq"
    bam_to_fastq(None, bam_path, reads_fastq)
    result = align_fastq_to_reference(
        reads_fastq,
        reference_path,
        out_dir,
        minimap2_preset=minimap2_preset,
        threads=threads,
        sort_memory=sort_memory,
        keep_intermediates=keep_intermediates,
        removable_intermediates=[reads_fastq],
    )
    result.update(
        {
            "input_type": "bam",
            "input_primary_reads": bam_stats["primary_reads"],
            "input_mapped_primary_reads": bam_stats["mapped_primary_reads"],
        }
    )
    return result


def run_fastq_pipeline(
    fastq_path,
    reference_path,
    out_dir,
    minimap2_preset="map-ont",
    threads=1,
    sort_memory="768M",
    keep_intermediates=False,
):
    fastq_stats = inspect_fastq_input(fastq_path)
    result = align_fastq_to_reference(
        Path(fastq_path),
        reference_path,
        out_dir,
        minimap2_preset=minimap2_preset,
        threads=threads,
        sort_memory=sort_memory,
        keep_intermediates=keep_intermediates,
    )
    result.update(
        {
            "input_type": "fastq",
            "input_primary_reads": fastq_stats["primary_reads"],
            "input_total_bases": fastq_stats["total_bases"],
            "input_mapped_primary_reads": 0,
        }
    )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, help="Input unaligned BAM")
    parser.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA file (.fa, .fasta, or .fna)",
    )
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument(
        "--minimap2-preset",
        default="map-ont",
        help="Preset passed to minimap2.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads for minimap2 and samtools subcommands. Default: 1.",
    )
    parser.add_argument(
        "--sort-memory",
        default="768M",
        help="Memory per samtools sort thread, for example 768M or 2G. Default: 768M.",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep reads.fastq, aligned.sam, and aligned.unsorted.bam for debugging.",
    )
    parser.add_argument(
        "--allow-aligned-input",
        action="store_true",
        help="Allow BAMs that already contain mapped primary reads.",
    )
    args = parser.parse_args()

    try:
        result = run_pipeline(
            bam_path=args.bam,
            reference_path=args.reference,
            out_dir=args.out_dir,
            minimap2_preset=args.minimap2_preset,
            threads=args.threads,
            sort_memory=args.sort_memory,
            keep_intermediates=args.keep_intermediates,
            allow_aligned_input=args.allow_aligned_input,
        )
    except ValueError as exc:
        parser.error(str(exc))

    for key, value in result.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
