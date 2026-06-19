#!/usr/bin/env python3
"""
Summarize an aligned BAM into per-base CSV output.

Output columns:
pos,base,depth,match_count,vaf,G,A,T,C,ins,del,qscore,confidence

This follows the example report semantics:
- covered positions only by default
- base is the majority observed nucleotide at each position
- depth uses primary, non-supplementary alignments and includes deletion-supporting reads
- vaf = match_count / depth
- qscore = mean BAM base quality for reads supporting the called base
"""

import argparse
import csv
import statistics
import warnings
from pathlib import Path

import pysam

PER_BASE_FIELDNAMES = [
    "pos",
    "base",
    "depth",
    "match_count",
    "vaf",
    "G",
    "A",
    "T",
    "C",
    "ins",
    "del",
    "qscore",
    "confidence",
]


def open_maybe_gzip(path):
    path = str(path)
    if path.endswith(".gz"):
        import gzip

        return gzip.open(path, "rt", encoding="utf-8-sig", errors="replace")
    return open(path, "r", encoding="utf-8-sig", errors="replace")


def parse_fasta(path):
    name = None
    chunks = []
    with open_maybe_gzip(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def parse_fastq(path):
    with open_maybe_gzip(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not (seq and plus and qual):
                raise ValueError("Truncated FASTQ file")
            header = header.rstrip("\n\r")
            seq = seq.rstrip("\n\r").upper()
            plus = plus.rstrip("\n\r")
            if not header.startswith("@"):
                raise ValueError(f"Bad FASTQ header: {header!r}")
            if not plus.startswith("+"):
                raise ValueError(f"Bad FASTQ separator: {plus!r}")
            yield header[1:].split()[0], seq


def parse_genbank(path):
    name = None
    collecting = False
    chunks = []
    with open(path, "r", encoding="utf-8-sig", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("LOCUS"):
                parts = line.split()
                name = parts[1] if len(parts) > 1 else "record"
            elif line.startswith("ORIGIN"):
                collecting = True
            elif line.startswith("//"):
                if name is not None and chunks:
                    yield name, "".join(chunks).upper()
                name = None
                collecting = False
                chunks = []
            elif collecting:
                seq = "".join(ch for ch in line if ch.isalpha())
                if seq:
                    chunks.append(seq)
    if name is not None and chunks:
        yield name, "".join(chunks).upper()


def load_reference_sequences(reference_path):
    if reference_path is None:
        return {}
    path = str(reference_path).lower()
    if path.endswith((".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")):
        return {name: seq for name, seq in parse_fasta(reference_path)}
    if path.endswith((".fq", ".fastq", ".fq.gz", ".fastq.gz")):
        return {name: seq for name, seq in parse_fastq(reference_path)}
    if path.endswith((".gb", ".gbk", ".genbank")):
        return {name: seq for name, seq in parse_genbank(reference_path)}
    raise ValueError(f"Unsupported reference format: {reference_path}")


def choose_consensus_base(base_counts, reference_base=None):
    max_count = max(base_counts.values())
    if max_count <= 0:
        return reference_base if reference_base in base_counts else "N"
    tied = [base for base, count in base_counts.items() if count == max_count]
    if reference_base in tied:
        return reference_base
    return sorted(tied)[0]


def compute_qscore(base_qualities):
    if not base_qualities:
        return 0
    return int(round(statistics.fmean(base_qualities)))


def confidence_label(
    depth,
    qscore,
    low_confidence_qscore,
    ambiguous_count=0,
    refskip_count=0,
    filtered_low_quality_count=0,
):
    flags = []
    if depth == 0:
        flags.append("ZERO_DEPTH")
    elif qscore < low_confidence_qscore:
        flags.append("LOW_QSCORE")
    else:
        flags.append("PASS")
    if ambiguous_count:
        flags.append(f"AMBIGUOUS={ambiguous_count}")
    if refskip_count:
        flags.append(f"REFSKIP={refskip_count}")
    if filtered_low_quality_count:
        flags.append(f"FILTERED_LOW_Q={filtered_low_quality_count}")
    return ";".join(flags)


def zero_depth_row(reference_pos, reference_seq=None):
    ref_base = reference_seq[reference_pos] if reference_seq is not None and reference_pos < len(reference_seq) else "N"
    return {
        "pos": reference_pos + 1,
        "base": ref_base,
        "depth": 0,
        "match_count": 0,
        "vaf": "0.000000",
        "G": 0,
        "A": 0,
        "T": 0,
        "C": 0,
        "ins": 0,
        "del": 0,
        "qscore": 0,
        "confidence": "ZERO_DEPTH",
    }


def build_row(
    reference_pos,
    pileup_column,
    reference_base=None,
    min_base_quality=0,
    low_confidence_qscore=12,
):
    base_counts = {base: 0 for base in "GATC"}
    base_qualities = {base: [] for base in "GATC"}
    depth = 0
    ins_count = 0
    del_count = 0
    ambiguous_count = 0
    refskip_count = 0
    filtered_low_quality_count = 0

    for pileup_read in pileup_column.pileups:
        read = pileup_read.alignment
        if read.is_secondary or read.is_supplementary:
            continue

        if pileup_read.is_refskip:
            refskip_count += 1
            continue

        if pileup_read.is_del:
            del_count += 1
            depth += 1
            continue

        query_pos = pileup_read.query_position
        if query_pos is None:
            continue

        query_sequence = read.query_sequence
        if not query_sequence or query_pos >= len(query_sequence):
            continue

        query_qualities = read.query_qualities
        quality = int(query_qualities[query_pos]) if query_qualities is not None else 0
        if quality < min_base_quality:
            filtered_low_quality_count += 1
            continue

        depth += 1
        query_base = query_sequence[query_pos].upper()
        if query_base in base_counts:
            base_counts[query_base] += 1
            base_qualities[query_base].append(quality)
        else:
            ambiguous_count += 1

        # `ins` is a read/event count anchored at this reference position,
        # not a count of inserted bases.
        if pileup_read.indel > 0:
            ins_count += 1

    consensus_base = choose_consensus_base(base_counts, reference_base=reference_base)
    match_count = base_counts.get(consensus_base, 0)
    vaf = (match_count / depth) if depth else 0.0
    qscore = compute_qscore(base_qualities.get(consensus_base, []))
    confidence = confidence_label(
        depth=depth,
        qscore=qscore,
        low_confidence_qscore=low_confidence_qscore,
        ambiguous_count=ambiguous_count,
        refskip_count=refskip_count,
        filtered_low_quality_count=filtered_low_quality_count,
    )

    return {
        "pos": reference_pos + 1,
        "base": consensus_base,
        "depth": depth,
        "match_count": match_count,
        "vaf": f"{vaf:.6f}",
        "G": base_counts["G"],
        "A": base_counts["A"],
        "T": base_counts["T"],
        "C": base_counts["C"],
        "ins": ins_count,
        "del": del_count,
        "qscore": qscore,
        "confidence": confidence,
    }


def empty_projected_position():
    return {
        "base_counts": {base: 0 for base in "GATC"},
        "base_qualities": {base: [] for base in "GATC"},
        "depth": 0,
        "ins_count": 0,
        "del_count": 0,
        "ambiguous_count": 0,
        "refskip_count": 0,
        "filtered_low_quality_count": 0,
    }


def projected_position_to_row(
    reference_pos,
    position,
    reference_base=None,
    low_confidence_qscore=12,
):
    base_counts = position["base_counts"]
    consensus_base = choose_consensus_base(base_counts, reference_base=reference_base)
    match_count = base_counts.get(consensus_base, 0)
    depth = position["depth"]
    vaf = (match_count / depth) if depth else 0.0
    qscore = compute_qscore(position["base_qualities"].get(consensus_base, []))
    confidence = confidence_label(
        depth=depth,
        qscore=qscore,
        low_confidence_qscore=low_confidence_qscore,
        ambiguous_count=position["ambiguous_count"],
        refskip_count=position["refskip_count"],
        filtered_low_quality_count=position["filtered_low_quality_count"],
    )
    return {
        "pos": reference_pos + 1,
        "base": consensus_base,
        "depth": depth,
        "match_count": match_count,
        "vaf": f"{vaf:.6f}",
        "G": base_counts["G"],
        "A": base_counts["A"],
        "T": base_counts["T"],
        "C": base_counts["C"],
        "ins": position["ins_count"],
        "del": position["del_count"],
        "qscore": qscore,
        "confidence": confidence,
    }


def iter_circular_projected_rows(
    bam,
    contig,
    reference_seq,
    min_base_quality=0,
    low_confidence_qscore=12,
):
    reference_length = len(reference_seq)
    if reference_length <= 0:
        raise ValueError("reference sequence must be non-empty")

    positions = [empty_projected_position() for _ in range(reference_length)]

    for read in bam.fetch(contig):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        query_sequence = read.query_sequence
        if not query_sequence:
            continue
        query_qualities = read.query_qualities
        last_projected_pos = None
        in_insertion = False

        for query_pos, ref_pos, cigar_op in read.get_aligned_pairs(matches_only=False, with_cigar=True):
            if ref_pos is None:
                if query_pos is not None and cigar_op == pysam.CINS and last_projected_pos is not None:
                    if not in_insertion:
                        positions[last_projected_pos]["ins_count"] += 1
                    in_insertion = True
                continue

            projected_pos = ref_pos % reference_length
            last_projected_pos = projected_pos
            in_insertion = False
            position = positions[projected_pos]

            if cigar_op == pysam.CREF_SKIP:
                position["refskip_count"] += 1
                continue

            if query_pos is None:
                position["del_count"] += 1
                position["depth"] += 1
                continue

            if query_pos >= len(query_sequence):
                continue

            quality = int(query_qualities[query_pos]) if query_qualities is not None else 0
            if quality < min_base_quality:
                position["filtered_low_quality_count"] += 1
                continue

            position["depth"] += 1
            query_base = query_sequence[query_pos].upper()
            if query_base in position["base_counts"]:
                position["base_counts"][query_base] += 1
                position["base_qualities"][query_base].append(quality)
            else:
                position["ambiguous_count"] += 1

    for reference_pos, position in enumerate(positions):
        yield projected_position_to_row(
            reference_pos=reference_pos,
            position=position,
            reference_base=reference_seq[reference_pos],
            low_confidence_qscore=low_confidence_qscore,
        )


def iter_rows(
    bam,
    contig,
    start,
    end,
    reference_seq=None,
    include_zero_depth=False,
    min_base_quality=0,
    low_confidence_qscore=12,
):
    next_zero_depth_pos = start
    for pileup_column in bam.pileup(
        contig=contig,
        start=start,
        end=end,
        stepper="all",
        truncate=True,
        min_base_quality=0,
    ):
        reference_pos = pileup_column.reference_pos
        if include_zero_depth:
            while next_zero_depth_pos < reference_pos:
                yield zero_depth_row(next_zero_depth_pos, reference_seq)
                next_zero_depth_pos += 1

        ref_base = None
        if reference_seq is not None and reference_pos < len(reference_seq):
            ref_base = reference_seq[reference_pos]
        yield build_row(
            reference_pos=reference_pos,
            pileup_column=pileup_column,
            reference_base=ref_base,
            min_base_quality=min_base_quality,
            low_confidence_qscore=low_confidence_qscore,
        )

        if include_zero_depth:
            next_zero_depth_pos = reference_pos + 1

    if include_zero_depth and end is None:
        if reference_seq is not None:
            end = len(reference_seq)
        else:
            end = bam.get_reference_length(contig)
    if include_zero_depth:
        while next_zero_depth_pos < end:
            yield zero_depth_row(next_zero_depth_pos, reference_seq)
            next_zero_depth_pos += 1


def resolve_contig(bam, contig):
    if contig:
        if contig not in bam.references:
            raise ValueError(f"Contig {contig!r} not found in BAM header")
        return contig
    if bam.nreferences == 1:
        return bam.references[0]
    raise ValueError("BAM has multiple references; pass --contig explicitly")


def summarize_bam_to_table(
    bam_path,
    output_csv,
    reference_path=None,
    reference_fastq_path=None,
    contig=None,
    start=None,
    end=None,
    include_zero_depth=False,
    low_confidence_out=None,
    low_confidence_qscore=12,
    min_base_quality=0,
):
    reference_path = reference_path or reference_fastq_path
    reference_lookup = load_reference_sequences(reference_path)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if bam.nreferences == 0:
            raise ValueError("BAM has no reference sequences (@SQ); align the reads first")

        contig = resolve_contig(bam, contig)
        reference_seq = reference_lookup.get(contig)
        if reference_lookup and reference_seq is None and len(reference_lookup) == 1:
            reference_name = next(iter(reference_lookup))
            warnings.warn(
                f"Reference contig {reference_name!r} does not match BAM contig {contig!r}; "
                "using the only reference sequence for tie-breaking and zero-depth rows.",
                RuntimeWarning,
                stacklevel=2,
            )
            reference_seq = next(iter(reference_lookup.values()))
        elif reference_lookup and reference_seq is None:
            warnings.warn(
                f"BAM contig {contig!r} was not found in the reference file; "
                "continuing without reference tie-breaking.",
                RuntimeWarning,
                stacklevel=2,
            )

        if start is None:
            start = 0
        if end is not None and end <= start:
            raise ValueError("end must be greater than start")

        low_writer = None
        low_handle = None
        if low_confidence_out is not None:
            low_handle = open(low_confidence_out, "w", newline="", encoding="ascii")
            low_writer = csv.DictWriter(low_handle, fieldnames=PER_BASE_FIELDNAMES)
            low_writer.writeheader()

        try:
            with open(output_csv, "w", newline="", encoding="ascii") as handle:
                writer = csv.DictWriter(handle, fieldnames=PER_BASE_FIELDNAMES)
                writer.writeheader()
                for row in iter_rows(
                    bam=bam,
                    contig=contig,
                    start=start,
                    end=end,
                    reference_seq=reference_seq,
                    include_zero_depth=include_zero_depth,
                    min_base_quality=min_base_quality,
                    low_confidence_qscore=low_confidence_qscore,
                ):
                    writer.writerow(row)
                    if low_writer is not None and int(row["qscore"]) < low_confidence_qscore:
                        low_writer.writerow(row)
        finally:
            if low_handle is not None:
                low_handle.close()


def summarize_circular_projected_bam_to_table(
    bam_path,
    output_csv,
    reference_path,
    contig=None,
    low_confidence_out=None,
    low_confidence_qscore=12,
    min_base_quality=0,
):
    reference_lookup = load_reference_sequences(reference_path)
    if not reference_lookup:
        raise ValueError(f"No reference sequences found in {reference_path}")
    if len(reference_lookup) != 1:
        raise ValueError("Circular projection currently requires a single-reference FASTA")
    reference_seq = next(iter(reference_lookup.values()))

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if bam.nreferences == 0:
            raise ValueError("BAM has no reference sequences (@SQ); align the reads first")
        contig = resolve_contig(bam, contig)

        low_writer = None
        low_handle = None
        if low_confidence_out is not None:
            low_handle = open(low_confidence_out, "w", newline="", encoding="ascii")
            low_writer = csv.DictWriter(low_handle, fieldnames=PER_BASE_FIELDNAMES)
            low_writer.writeheader()

        try:
            with open(output_csv, "w", newline="", encoding="ascii") as handle:
                writer = csv.DictWriter(handle, fieldnames=PER_BASE_FIELDNAMES)
                writer.writeheader()
                for row in iter_circular_projected_rows(
                    bam=bam,
                    contig=contig,
                    reference_seq=reference_seq,
                    min_base_quality=min_base_quality,
                    low_confidence_qscore=low_confidence_qscore,
                ):
                    writer.writerow(row)
                    if low_writer is not None and int(row["qscore"]) < low_confidence_qscore:
                        low_writer.writerow(row)
        finally:
            if low_handle is not None:
                low_handle.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, help="Aligned input BAM")
    parser.add_argument("--out", required=True, help="Output CSV path")
    parser.add_argument(
        "--reference",
        default=None,
        help="Optional reference (.fasta, .fastq, or .gbk) used for tie-breaking and zero-depth rows",
    )
    parser.add_argument("--contig", default=None, help="Contig to summarize")
    parser.add_argument("--start", type=int, default=None, help="0-based start position")
    parser.add_argument("--end", type=int, default=None, help="0-based end position (exclusive)")
    parser.add_argument(
        "--include-zero-depth",
        action="store_true",
        help="Emit zero-depth positions within the selected interval",
    )
    parser.add_argument(
        "--low-confidence-out",
        default=None,
        help="Optional CSV path for rows below --low-confidence-qscore",
    )
    parser.add_argument(
        "--low-confidence-qscore",
        type=int,
        default=12,
        help="Rows below this mean BAM base-quality threshold are copied to --low-confidence-out",
    )
    parser.add_argument(
        "--min-base-quality",
        type=int,
        default=0,
        help="Minimum base quality threshold for pileup counting",
    )
    parser.add_argument(
        "--circular-project",
        action="store_true",
        help=(
            "Interpret the aligned BAM as reads mapped to a doubled circular reference and "
            "project per-base counts back onto --reference coordinates"
        ),
    )
    args = parser.parse_args()

    if args.circular_project:
        if args.reference is None:
            raise ValueError("--circular-project requires --reference")
        summarize_circular_projected_bam_to_table(
            bam_path=args.bam,
            output_csv=args.out,
            reference_path=args.reference,
            contig=args.contig,
            low_confidence_out=args.low_confidence_out,
            low_confidence_qscore=args.low_confidence_qscore,
            min_base_quality=args.min_base_quality,
        )
    else:
        summarize_bam_to_table(
            bam_path=args.bam,
            output_csv=args.out,
            reference_path=args.reference,
            contig=args.contig,
            start=args.start,
            end=args.end,
            include_zero_depth=args.include_zero_depth,
            low_confidence_out=args.low_confidence_out,
            low_confidence_qscore=args.low_confidence_qscore,
            min_base_quality=args.min_base_quality,
        )


if __name__ == "__main__":
    main()
