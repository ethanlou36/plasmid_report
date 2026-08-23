from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import (
    BeforePosition,
    CompoundLocation,
    ExactPosition,
    SeqFeature,
    SimpleLocation,
)
from Bio.SeqRecord import SeqRecord

from plasmid_normalization import (
    DEFAULT_ORIGIN_MOTIF,
    OriginNormalizationError,
    apply_transform_to_fastq,
    canonicalize_sequence,
    normalize_genbank_record,
    normalize_plasmid_files,
    normalize_plasmid_pair,
)


MOTIF = DEFAULT_ORIGIN_MOTIF


def make_record(sequence: str, features=None, topology="circular") -> SeqRecord:
    record = SeqRecord(Seq(sequence), id="test", name="test", description="test plasmid")
    record.annotations["molecule_type"] = "DNA"
    if topology is not None:
        record.annotations["topology"] = topology
    record.features = list(features or [])
    return record


def test_canonicalize_sequence_handles_noop_rotation_and_reverse_complement():
    canonical = MOTIF + "ACGCGT"

    unchanged, unchanged_transform = canonicalize_sequence(canonical)
    rotated, rotated_transform = canonicalize_sequence("GCGT" + MOTIF + "AC")
    reverse_input = str(Seq("GG" + canonical).reverse_complement())
    reverse, reverse_transform = canonicalize_sequence(reverse_input)

    assert unchanged == canonical
    assert unchanged_transform.status == "already_at_origin"
    assert rotated == MOTIF + "ACGCGT"
    assert rotated_transform.rotation_offset_bp == 4
    assert reverse == canonical + "GG"
    assert reverse_transform.orientation == "reverse_complement"
    assert reverse_transform.rotation_offset_bp == 2


def test_canonicalize_sequence_finds_motif_across_old_origin():
    sequence = MOTIF[5:] + "ACGCGC" + MOTIF[:5]

    canonical, transform = canonicalize_sequence(sequence)

    assert canonical == MOTIF + "ACGCGC"
    assert transform.orientation == "forward"
    assert transform.rotation_offset_bp == 13


def test_missing_and_ambiguous_motifs_fail_explicitly():
    with pytest.raises(OriginNormalizationError, match="was not found"):
        canonicalize_sequence("ACGCGTACGCGT")

    with pytest.raises(OriginNormalizationError, match="multiple distinct"):
        canonicalize_sequence(MOTIF + "AC" + MOTIF + "GT")


def test_fastq_qualities_follow_forward_and_reverse_complement_transforms():
    forward_input = "GG" + MOTIF + "AC"
    forward_quality = "ABCDEFGHIJKLMNOP"
    canonical, transform = canonicalize_sequence(forward_input)
    transformed_sequence, transformed_quality = apply_transform_to_fastq(
        forward_input,
        forward_quality,
        transform,
    )
    assert transformed_sequence == canonical
    assert transformed_quality == forward_quality[2:] + forward_quality[:2]

    reverse_input = str(Seq(forward_input).reverse_complement())
    reverse_quality = "abcdefghijklmnop"
    reverse_canonical, reverse_transform = canonicalize_sequence(reverse_input)
    reverse_sequence, transformed_reverse_quality = apply_transform_to_fastq(
        reverse_input,
        reverse_quality,
        reverse_transform,
    )
    expected_oriented_quality = reverse_quality[::-1]
    expected_quality = expected_oriented_quality[2:] + expected_oriented_quality[:2]
    assert reverse_sequence == reverse_canonical == canonical
    assert transformed_reverse_quality == expected_quality


def test_genbank_rotation_preserves_crossing_feature_extraction_and_source():
    prefix = "GCGCG"
    sequence = prefix + MOTIF + "AACCGGTTAACC"
    length = len(sequence)
    features = [
        SeqFeature(SimpleLocation(0, length, strand=1), type="source"),
        SeqFeature(SimpleLocation(2, 10, strand=1), type="misc_feature", qualifiers={"label": ["plus"]}),
        SeqFeature(SimpleLocation(2, 10, strand=-1), type="CDS", qualifiers={"label": ["minus"]}),
        SeqFeature(
            CompoundLocation(
                [SimpleLocation(1, 3, strand=1), SimpleLocation(15, 18, strand=1)],
                operator="join",
            ),
            type="misc_feature",
            qualifiers={"note": ["compound"]},
        ),
    ]
    record = make_record(sequence, features)
    _, transform = canonicalize_sequence(sequence)

    normalized = normalize_genbank_record(record, transform)

    assert str(normalized.seq).startswith(MOTIF)
    assert normalized.features[0].location == SimpleLocation(0, length, strand=1)
    assert isinstance(normalized.features[1].location, CompoundLocation)
    assert isinstance(normalized.features[2].location, CompoundLocation)
    for before, after in zip(record.features, normalized.features):
        if before.type != "source":
            assert str(before.extract(record.seq)).upper() == str(after.extract(normalized.seq)).upper()
        assert before.qualifiers == after.qualifiers


def test_genbank_reverse_complement_flips_strand_and_preserves_feature_sequence():
    canonical = MOTIF + "AACCGGTT"
    input_sequence = str(Seq("GG" + canonical).reverse_complement())
    feature = SeqFeature(
        SimpleLocation(3, 10, strand=1),
        type="misc_feature",
        qualifiers={"label": ["oriented feature"]},
    )
    record = make_record(input_sequence, [feature])
    record.letter_annotations["phred_quality"] = list(range(len(input_sequence)))
    _, transform = canonicalize_sequence(input_sequence)

    normalized = normalize_genbank_record(record, transform)

    assert transform.status == "reverse_complemented_and_rotated"
    assert str(normalized.seq) == canonical + "GG"
    assert normalized.features[0].location.strand == -1
    oriented_quality = list(reversed(record.letter_annotations["phred_quality"]))
    assert normalized.letter_annotations["phred_quality"] == (
        oriented_quality[2:] + oriented_quality[:2]
    )
    assert str(feature.extract(record.seq)).upper() == str(
        normalized.features[0].extract(normalized.seq)
    ).upper()


def test_attached_coordinate_oracle_rep_origin_8043_becomes_base_one():
    length = 10822
    sequence = "A" * 8042 + MOTIF + "C" * (length - 8042 - len(MOTIF))
    feature = SeqFeature(
        SimpleLocation(8042, 8631, strand=1),
        type="rep_origin",
        qualifiers={"label": ["ori"]},
    )
    record = make_record(sequence, [feature])
    _, transform = canonicalize_sequence(sequence)

    normalized = normalize_genbank_record(record, transform)

    assert transform.rotation_offset_bp == 8042
    assert normalized.features[0].location == SimpleLocation(0, 589, strand=1)


def test_rotation_rejects_fuzzy_locations_and_non_circular_records():
    sequence = "GG" + MOTIF + "ACGT"
    fuzzy = SeqFeature(
        SimpleLocation(BeforePosition(0), ExactPosition(5), strand=1),
        type="misc_feature",
    )
    _, transform = canonicalize_sequence(sequence)

    with pytest.raises(OriginNormalizationError, match="fuzzy or remote"):
        normalize_genbank_record(make_record(sequence, [fuzzy]), transform)
    with pytest.raises(OriginNormalizationError, match="not marked circular"):
        normalize_genbank_record(make_record(sequence, topology="linear"), transform)


def test_pair_normalization_accepts_different_rotations_and_rejects_disagreement():
    canonical = MOTIF + "ACGCGTTA"
    fasta_sequence = canonical[3:] + canonical[:3]
    genbank_sequence = canonical[7:] + canonical[:7]
    record = make_record(
        genbank_sequence,
        [SeqFeature(SimpleLocation(0, len(genbank_sequence)), type="source")],
    )

    normalized = normalize_plasmid_pair(fasta_sequence, record, sequence_name="sample")

    assert normalized.sequence == canonical
    assert str(normalized.genbank_record.seq) == canonical
    assert normalized.fasta_transform.rotation_offset_bp == len(canonical) - 3
    assert normalized.genbank_transform.rotation_offset_bp == len(canonical) - 7

    mismatched_canonical = MOTIF + "ACGCGTTC"
    mismatched = make_record(mismatched_canonical[7:] + mismatched_canonical[:7])
    with pytest.raises(OriginNormalizationError, match="do not match"):
        normalize_plasmid_pair(fasta_sequence, mismatched)


def test_normalize_plasmid_files_writes_matching_valid_outputs(tmp_path: Path):
    canonical = MOTIF + "ACGCGTTA"
    input_sequence = "GG" + canonical
    fasta_path = tmp_path / "input.fasta"
    gbk_path = tmp_path / "input.gbk"
    SeqIO.write(SeqRecord(Seq(input_sequence), id="old", description=""), fasta_path, "fasta")
    record = make_record(
        input_sequence,
        [
            SeqFeature(SimpleLocation(0, len(input_sequence)), type="source"),
            SeqFeature(
                SimpleLocation(1, 5, strand=1),
                type="misc_feature",
                qualifiers={"label": ["unicode α"]},
            ),
        ],
    )
    SeqIO.write(record, gbk_path, "genbank")
    fasta_out = tmp_path / "package" / "sample.fa"
    gbk_out = tmp_path / "package" / "sample.gbk"

    result = normalize_plasmid_files(
        fasta_path,
        gbk_path,
        fasta_out,
        gbk_out,
        "001_sample_contig",
    )

    written_fasta = SeqIO.read(fasta_out, "fasta")
    written_genbank = SeqIO.read(gbk_out, "genbank")
    assert str(written_fasta.seq) == str(written_genbank.seq) == canonical + "GG"
    assert written_genbank.name == "001_sample_conti"
    assert written_genbank.features[1].qualifiers["label"] == ["unicode α"]
    assert result["length_bp"] == len(input_sequence)
    assert result["origin_normalization"]["fasta"]["rotation_offset_bp"] == 2
    assert len(result["origin_normalization"]["canonical_sequence_sha256"]) == 64
