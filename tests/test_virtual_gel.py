import json
from pathlib import Path

import numpy as np
import pysam

from epi2me_to_final_package import cleanup_previous_sample_output, generate_order_virtual_gels
from virtual_gel import make_virtual_gel, resolve_virtual_gel_y_max, select_virtual_gel_bands


def write_unaligned_bam(path: Path) -> None:
    with pysam.AlignmentFile(path, "wb", header={"HD": {"VN": "1.0"}}) as bam:
        for name, length in (("read1", 4800), ("read2", 5000), ("read3", 12795)):
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "A" * length
            read.flag = 4
            read.query_qualities = pysam.qualitystring_to_array("I" * length)
            bam.write(read)


def test_make_virtual_gel_writes_png(tmp_path: Path):
    out_path = tmp_path / "virtual_gel.png"
    make_virtual_gel(
        {
            "001 - sample-a": [900, 1100, 4800, 5000, 5050, 9900],
            "002 - sample-b": [1200, 4300, 4400, 4500, 8700, 8800],
        },
        out_path,
        title="Virtual Gel - WPS Order TEST",
    )

    assert out_path.exists()
    assert out_path.stat().st_size > 0


def test_virtual_gel_y_max_uses_full_displayed_range():
    y_max = resolve_virtual_gel_y_max(
        [
            # Rare long reads should still be visible on the gel axis.
            # This guards against percentile clipping.
            np.array([4800, 5000, 5050, 12795], dtype=float)
        ],
        y_max=None,
        min_display_bp=1000,
    )

    assert y_max == 13000


def test_virtual_gel_selects_expected_multimer_bands_and_drops_smear():
    bands = select_virtual_gel_bands(
        [5000] * 100 + [7500] * 40 + [10000] * 20 + [20000] * 3,
        expected_length_bp=5000,
    )

    assert [band.multiple for band in bands] == [1, 2]
    assert all(abs(band.center_bp - expected) < 100 for band, expected in zip(bands, [5000, 10000]))


def test_generate_order_virtual_gel_uses_package_summary_fastq(tmp_path: Path):
    output_root = tmp_path / "output"
    order_dir = output_root / "WPS Data_Order #TEST"
    work_dir = output_root / "_work" / "001_sample"
    work_dir.mkdir(parents=True)

    fastq_path = tmp_path / "reads.fastq"
    fastq_path.write_text(
        "@read1\n" + "A" * 4800 + "\n+\n" + "I" * 4800 + "\n"
        "@read2\n" + "A" * 5000 + "\n+\n" + "I" * 5000 + "\n",
        encoding="ascii",
    )
    (work_dir / "package_summary.json").write_text(
        json.dumps(
            {
                "order_number": "TEST",
                "sample_name": "001_sample",
                "virtual_gel_label": "001 - sample",
                "order_sample_number": "1",
                "contig_length_bp": 5000,
                "paths": {
                    "order_dir": str(order_dir),
                    "alignment_input": str(fastq_path),
                    "alignment_input_type": "fastq",
                },
            }
        ),
        encoding="utf-8",
    )

    virtual_gels, warnings = generate_order_virtual_gels(output_root, ["TEST"])

    assert warnings == []
    assert "TEST" in virtual_gels
    assert Path(virtual_gels["TEST"]).exists()
    assert Path(virtual_gels["TEST"]).parent == order_dir


def test_generate_order_virtual_gel_uses_package_summary_bam(tmp_path: Path):
    output_root = tmp_path / "output"
    order_dir = output_root / "WPS Data_Order #TEST"
    work_dir = output_root / "_work" / "001_sample"
    work_dir.mkdir(parents=True)

    bam_path = tmp_path / "reads.bam"
    write_unaligned_bam(bam_path)
    (work_dir / "package_summary.json").write_text(
        json.dumps(
            {
                "order_number": "TEST",
                "sample_name": "001_sample",
                "virtual_gel_label": "001 - sample",
                "order_sample_number": "1",
                "contig_length_bp": 5000,
                "paths": {
                    "order_dir": str(order_dir),
                    "alignment_input": str(bam_path),
                    "alignment_input_type": "bam",
                },
            }
        ),
        encoding="utf-8",
    )

    virtual_gels, warnings = generate_order_virtual_gels(output_root, ["TEST"])

    assert warnings == []
    assert "TEST" in virtual_gels
    assert Path(virtual_gels["TEST"]).exists()
    assert Path(virtual_gels["TEST"]).parent == order_dir


def test_cleanup_previous_sample_output_removes_stale_order_virtual_gel(tmp_path: Path):
    output_root = tmp_path / "output"
    order_dir = output_root / "WPS Data_Order #TEST"
    qc_dir = order_dir / "QC REPORTS"
    qc_dir.mkdir(parents=True)
    old_stale_gel = qc_dir / "Order_TEST_virtual_gel.png"
    old_stale_gel.write_bytes(b"stale")
    stale_gel = order_dir / "Order_TEST_virtual_gel.png"
    stale_gel.write_bytes(b"stale")

    work_dir = output_root / "_work" / "001_sample"
    work_dir.mkdir(parents=True)
    (work_dir / "package_summary.json").write_text(
        json.dumps(
            {
                "barcode": "barcode01",
                "sample_name": "001_sample",
                "paths": {
                    "order_dir": str(order_dir),
                    "pdf": str(qc_dir / "001_sample_report.pdf"),
                },
            }
        ),
        encoding="utf-8",
    )

    cleanup_previous_sample_output(output_root, "barcode01", "001_sample")

    assert not stale_gel.exists()
    assert not old_stale_gel.exists()
