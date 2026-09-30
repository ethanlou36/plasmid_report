from pathlib import Path

from epi2me_to_final_package import raw_fastq_rejection_reason


def write_fastq(path: Path, record_count: int, read_length: int) -> None:
    sequence = "A" * read_length
    quality = "I" * read_length
    with path.open("w", encoding="ascii") as handle:
        for index in range(record_count):
            handle.write(f"@read{index}\n{sequence}\n+\n{quality}\n")


def test_twenty_reads_and_sixty_thousand_bases_are_eligible(tmp_path: Path):
    fastq = tmp_path / "barcode01.fastq"
    write_fastq(fastq, record_count=20, read_length=3_000)

    assert raw_fastq_rejection_reason(fastq) is None


def test_fewer_than_twenty_reads_are_rejected_even_with_enough_bases(tmp_path: Path):
    fastq = tmp_path / "barcode01.fastq"
    write_fastq(fastq, record_count=19, read_length=3_200)

    reason = raw_fastq_rejection_reason(fastq)

    assert "only 19 FASTQ records" in reason
    assert "expected at least 20" in reason


def test_twenty_reads_are_rejected_when_total_bases_are_too_low(tmp_path: Path):
    fastq = tmp_path / "barcode01.fastq"
    write_fastq(fastq, record_count=20, read_length=2_999)

    reason = raw_fastq_rejection_reason(fastq)

    assert "only 59,980 read bases" in reason
    assert "expected at least 60,000" in reason
