from pathlib import Path

from epi2me_to_final_package import raw_fastq_rejection_reason


def write_fastq(path: Path, record_count: int, read_length: int) -> None:
    sequence = "A" * read_length
    quality = "I" * read_length
    with path.open("w", encoding="ascii") as handle:
        for index in range(record_count):
            handle.write(f"@read{index}\n{sequence}\n+\n{quality}\n")


def test_fifty_reads_and_one_hundred_thousand_bases_are_eligible(tmp_path: Path):
    fastq = tmp_path / "barcode01.fastq"
    write_fastq(fastq, record_count=50, read_length=2_000)

    assert raw_fastq_rejection_reason(fastq) is None


def test_fewer_than_fifty_reads_are_rejected_even_with_enough_bases(tmp_path: Path):
    fastq = tmp_path / "barcode01.fastq"
    write_fastq(fastq, record_count=49, read_length=2_100)

    reason = raw_fastq_rejection_reason(fastq)

    assert "only 49 FASTQ records" in reason
    assert "expected at least 50" in reason


def test_fifty_reads_are_rejected_when_total_bases_are_too_low(tmp_path: Path):
    fastq = tmp_path / "barcode01.fastq"
    write_fastq(fastq, record_count=50, read_length=1_000)

    reason = raw_fastq_rejection_reason(fastq)

    assert "only 50,000 read bases" in reason
    assert "expected at least 100,000" in reason
