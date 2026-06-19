#!/usr/bin/env python3
import csv
import sys
import tempfile
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bam_to_per_base_data import summarize_circular_projected_bam_to_table


def write_cross_boundary_bam(path: Path) -> None:
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"LN": 16, "SN": "plasmid__circular_doubled"}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        read = pysam.AlignedSegment()
        read.query_name = "crosses_artificial_breakpoint"
        read.query_sequence = "GTAC"
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 6
        read.mapping_quality = 60
        read.cigar = [(0, 4)]
        read.query_qualities = pysam.qualitystring_to_array("????")
        bam.write(read)
    pysam.index(str(path))


def write_filtering_and_deletion_bam(path: Path) -> None:
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"LN": 8, "SN": "plasmid__circular_doubled"}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        primary = pysam.AlignedSegment()
        primary.query_name = "primary_with_deletion"
        primary.query_sequence = "ACT"
        primary.flag = 0
        primary.reference_id = 0
        primary.reference_start = 0
        primary.mapping_quality = 60
        primary.cigar = [(0, 2), (2, 1), (0, 1)]
        primary.query_qualities = pysam.qualitystring_to_array("???")
        bam.write(primary)

        secondary = pysam.AlignedSegment()
        secondary.query_name = "secondary_ignored"
        secondary.query_sequence = "GGGG"
        secondary.flag = 256
        secondary.reference_id = 0
        secondary.reference_start = 0
        secondary.mapping_quality = 60
        secondary.cigar = [(0, 4)]
        secondary.query_qualities = pysam.qualitystring_to_array("????")
        bam.write(secondary)

        supplementary = pysam.AlignedSegment()
        supplementary.query_name = "supplementary_ignored"
        supplementary.query_sequence = "TTTT"
        supplementary.flag = 2048
        supplementary.reference_id = 0
        supplementary.reference_start = 0
        supplementary.mapping_quality = 60
        supplementary.cigar = [(0, 4)]
        supplementary.query_qualities = pysam.qualitystring_to_array("????")
        bam.write(supplementary)
    pysam.index(str(path))


def test_circular_projection_wraps_doubled_reference_positions():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        reference = tmpdir / "plasmid.fa"
        reference.write_text(">plasmid\nACGTACGT\n", encoding="ascii")
        bam = tmpdir / "aligned.bam"
        write_cross_boundary_bam(bam)

        per_base_csv = tmpdir / "per_base.csv"
        low_conf_csv = tmpdir / "low_conf.csv"
        summarize_circular_projected_bam_to_table(
            bam_path=bam,
            output_csv=per_base_csv,
            reference_path=reference,
            contig="plasmid__circular_doubled",
            low_confidence_out=low_conf_csv,
            low_confidence_qscore=12,
        )

        with per_base_csv.open("r", encoding="ascii", newline="") as handle:
            rows = list(csv.DictReader(handle))

        assert [int(row["depth"]) for row in rows] == [1, 1, 0, 0, 0, 0, 1, 1]
        assert [row["base"] for row in rows] == ["A", "C", "G", "T", "A", "C", "G", "T"]
        assert rows[0]["qscore"] == "30"
        assert rows[1]["qscore"] == "30"
        assert rows[6]["qscore"] == "30"
        assert rows[7]["qscore"] == "30"


def test_circular_projection_counts_deletions_and_skips_secondary_alignments():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        reference = tmpdir / "plasmid.fa"
        reference.write_text(">plasmid\nACGT\n", encoding="ascii")
        bam = tmpdir / "aligned.bam"
        write_filtering_and_deletion_bam(bam)

        per_base_csv = tmpdir / "per_base.csv"
        summarize_circular_projected_bam_to_table(
            bam_path=bam,
            output_csv=per_base_csv,
            reference_path=reference,
            contig="plasmid__circular_doubled",
            low_confidence_qscore=12,
        )

        with per_base_csv.open("r", encoding="ascii", newline="") as handle:
            rows = list(csv.DictReader(handle))

        assert [int(row["depth"]) for row in rows] == [1, 1, 1, 1]
        assert [int(row["del"]) for row in rows] == [0, 0, 1, 0]
        assert [int(row["G"]) for row in rows] == [0, 0, 0, 0]
        assert rows[2]["base"] == "G"
        assert rows[2]["match_count"] == "0"
        assert rows[2]["qscore"] == "0"


if __name__ == "__main__":
    test_circular_projection_wraps_doubled_reference_positions()
    test_circular_projection_counts_deletions_and_skips_secondary_alignments()
