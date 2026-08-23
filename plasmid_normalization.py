"""Canonicalize circular plasmid FASTA and GenBank records to one origin.

The sequence transformation and the GenBank feature transformation live here so
that callers cannot accidentally rotate the two deliverables independently.
Coordinates are handled with Biopython's structured location objects; GenBank
location text is never rewritten with regular expressions.
"""

from __future__ import annotations

import copy
import hashlib
import os
import re
import tempfile
import warnings
from dataclasses import dataclass
from pathlib import Path

from Bio import BiopythonParserWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, ExactPosition, SimpleLocation
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqRecord import SeqRecord


DEFAULT_ORIGIN_MOTIF = "TTGAGATCCTTT"


class OriginNormalizationError(ValueError):
    """Raised when a plasmid cannot be normalized without guessing."""


@dataclass(frozen=True)
class CircularTransform:
    """A reverse-complement choice followed by a zero-based circular cut."""

    orientation: str
    rotation_offset_bp: int

    @property
    def reverse_complemented(self) -> bool:
        return self.orientation == "reverse_complement"

    @property
    def status(self) -> str:
        if not self.reverse_complemented and self.rotation_offset_bp == 0:
            return "already_at_origin"
        if not self.reverse_complemented:
            return "rotated"
        if self.rotation_offset_bp == 0:
            return "reverse_complemented"
        return "reverse_complemented_and_rotated"

    def as_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "orientation": self.orientation,
            "reverse_complemented": self.reverse_complemented,
            "rotation_offset_bp": self.rotation_offset_bp,
            "bases_moved_from_front": self.rotation_offset_bp,
            "motif_start_in_oriented_sequence_1based": self.rotation_offset_bp + 1,
        }


@dataclass
class NormalizedPlasmidPair:
    sequence: str
    fasta_transform: CircularTransform
    genbank_transform: CircularTransform
    genbank_record: SeqRecord


def _validated_sequence_and_motif(sequence: str, motif: str, label: str) -> tuple[str, str]:
    sequence = str(sequence).upper()
    motif = str(motif).upper()
    if not sequence:
        raise OriginNormalizationError(f"{label} is empty")
    if not motif:
        raise OriginNormalizationError("origin motif is empty")
    if set(motif) - set("ACGT"):
        raise OriginNormalizationError(
            f"origin motif must contain only A, C, G, and T; got {motif!r}"
        )
    if len(motif) > len(sequence):
        raise OriginNormalizationError(
            f"origin motif {motif!r} is longer than {label} ({len(sequence)} bp)"
        )
    return sequence, motif


def _circular_match_offsets(sequence: str, motif: str) -> list[int]:
    extended = sequence + sequence[: len(motif) - 1]
    return [i for i in range(len(sequence)) if extended[i : i + len(motif)] == motif]


def _rotate(values, offset: int):
    return values[offset:] + values[:offset]


def apply_transform_to_sequence(sequence: str, transform: CircularTransform) -> str:
    sequence = str(sequence).upper()
    if transform.reverse_complemented:
        sequence = str(Seq(sequence).reverse_complement())
    offset = transform.rotation_offset_bp
    if not 0 <= offset < len(sequence):
        raise OriginNormalizationError(
            f"rotation offset {offset} is outside a {len(sequence)} bp sequence"
        )
    return _rotate(sequence, offset)


def apply_transform_to_fastq(
    sequence: str,
    qualities: str,
    transform: CircularTransform,
) -> tuple[str, str]:
    """Apply a plasmid transform to a FASTQ sequence and parallel qualities."""
    if len(sequence) != len(qualities):
        raise OriginNormalizationError(
            f"FASTQ sequence and quality lengths differ ({len(sequence)} != {len(qualities)})"
        )
    normalized_sequence = apply_transform_to_sequence(sequence, transform)
    if transform.reverse_complemented:
        qualities = qualities[::-1]
    return normalized_sequence, _rotate(qualities, transform.rotation_offset_bp)


def _candidate_map(
    sequence: str,
    motif: str,
    label: str,
) -> dict[str, list[CircularTransform]]:
    sequence, motif = _validated_sequence_and_motif(sequence, motif, label)
    candidates: dict[str, list[CircularTransform]] = {}
    orientations = (
        ("forward", sequence),
        ("reverse_complement", str(Seq(sequence).reverse_complement())),
    )
    for orientation, oriented_sequence in orientations:
        for offset in _circular_match_offsets(oriented_sequence, motif):
            transform = CircularTransform(orientation, offset)
            canonical = _rotate(oriented_sequence, offset)
            candidates.setdefault(canonical, []).append(transform)
    if not candidates:
        raise OriginNormalizationError(
            f"origin motif {motif!r} was not found on either strand of {label}"
        )
    return candidates


def _preferred_transform(transforms: list[CircularTransform]) -> CircularTransform:
    return min(
        transforms,
        key=lambda item: (
            item.reverse_complemented,
            item.rotation_offset_bp != 0,
            item.rotation_offset_bp,
        ),
    )


def _ambiguous_message(
    candidates: dict[str, list[CircularTransform]],
    motif: str,
    label: str,
) -> str:
    hits = [
        f"{transform.orientation} position {transform.rotation_offset_bp + 1}"
        for transforms in candidates.values()
        for transform in transforms
    ]
    return (
        f"origin motif {motif.upper()!r} gives multiple distinct canonical starts "
        f"for {label} ({'; '.join(hits)}); refusing to choose one"
    )


def canonicalize_sequence(
    sequence: str,
    motif: str = DEFAULT_ORIGIN_MOTIF,
    *,
    label: str = "sequence",
) -> tuple[str, CircularTransform]:
    """Return one unambiguous circular sequence beginning with ``motif``."""
    candidates = _candidate_map(sequence, motif, label)
    if len(candidates) != 1:
        raise OriginNormalizationError(_ambiguous_message(candidates, motif, label))
    canonical, transforms = next(iter(candidates.items()))
    return canonical, _preferred_transform(transforms)


def _has_only_exact_local_positions(location) -> bool:
    return location is None or all(
        part.ref is None
        and part.ref_db is None
        and type(part.start) is ExactPosition
        and type(part.end) is ExactPosition
        for part in location.parts
    )


def _has_remote_reference(location) -> bool:
    return location is not None and any(
        part.ref is not None or part.ref_db is not None for part in location.parts
    )


def _new_simple_location(
    start: int,
    end: int,
    source: SimpleLocation,
) -> SimpleLocation:
    return SimpleLocation(ExactPosition(start), ExactPosition(end), strand=source.strand)


def _rotate_simple_location(
    location: SimpleLocation,
    offset: int,
    sequence_length: int,
) -> list[SimpleLocation]:
    start = int(location.start)
    end = int(location.end)
    if start == 0 and end == sequence_length:
        return [_new_simple_location(0, sequence_length, location)]
    if start == end:
        point = (start - offset) % sequence_length
        return [_new_simple_location(point, point, location)]
    if start >= offset:
        return [_new_simple_location(start - offset, end - offset, location)]
    if end <= offset:
        return [
            _new_simple_location(
                start + sequence_length - offset,
                end + sequence_length - offset,
                location,
            )
        ]

    # The feature crosses the new origin. CompoundLocation concatenates parts in
    # list order, so negative-strand parts use the opposite order to preserve
    # feature.extract() exactly.
    high = _new_simple_location(start + sequence_length - offset, sequence_length, location)
    low = _new_simple_location(0, end - offset, location)
    return [low, high] if location.strand == -1 else [high, low]


def _rotate_location(location, offset: int, sequence_length: int):
    if location is None or offset == 0:
        return copy.deepcopy(location)
    if not _has_only_exact_local_positions(location):
        raise OriginNormalizationError(
            f"cannot safely rotate unsupported fuzzy or remote GenBank location {location}"
        )

    compound = isinstance(location, CompoundLocation)
    original_parts = location.parts if compound else [location]
    parts = [
        rotated
        for part in original_parts
        for rotated in _rotate_simple_location(part, offset, sequence_length)
    ]
    return parts[0] if len(parts) == 1 else CompoundLocation(
        parts, operator=location.operator if compound else "join"
    )


def normalize_genbank_record(
    record: SeqRecord,
    transform: CircularTransform,
    *,
    sequence_name: str | None = None,
) -> SeqRecord:
    """Apply ``transform`` while preserving feature extraction semantics."""
    sequence_length = len(record)
    if sequence_length == 0:
        raise OriginNormalizationError("GenBank record is empty")

    topology = str(record.annotations.get("topology", "")).lower()
    changed = transform.reverse_complemented or transform.rotation_offset_bp != 0
    if changed and topology != "circular":
        raise OriginNormalizationError(
            "GenBank record is not marked circular; refusing to rotate or reverse-complement it"
        )
    if changed and any(_has_remote_reference(feature.location) for feature in record.features):
        raise OriginNormalizationError(
            "cannot safely transform a GenBank feature that references a remote sequence"
        )

    if transform.reverse_complemented:
        oriented = record.reverse_complement(
            id=True,
            name=True,
            description=True,
            features=True,
            annotations=True,
            letter_annotations=True,
            dbxrefs=True,
        )
    else:
        oriented = copy.deepcopy(record)

    offset = transform.rotation_offset_bp
    normalized_sequence = apply_transform_to_sequence(str(record.seq), transform)
    unsupported = next(
        (
            feature.location
            for feature in oriented.features
            if not _has_only_exact_local_positions(feature.location)
        ),
        None,
    )
    if offset and unsupported is not None:
        raise OriginNormalizationError(
            f"cannot safely rotate unsupported fuzzy or remote GenBank location {unsupported}"
        )

    normalized = copy.deepcopy(oriented)
    normalized.seq = Seq(normalized_sequence)
    normalized.letter_annotations = {
        key: copy.deepcopy(_rotate(values, offset))
        for key, values in oriented.letter_annotations.items()
    }
    normalized.features = copy.deepcopy(oriented.features)
    for feature in normalized.features:
        feature.location = _rotate_location(feature.location, offset, sequence_length)

    if sequence_name:
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", sequence_name).strip("_")[:16]
        normalized.name = safe_name or "record"
    normalized.annotations.setdefault("molecule_type", "DNA")

    for index, (before, after) in enumerate(zip(oriented.features, normalized.features)):
        if (
            before.type == "source"
            and before.location is not None
            and int(before.location.start) == 0
            and int(before.location.end) == sequence_length
        ):
            # A whole-record source still covers the entire molecule, but its
            # extracted linear string necessarily follows the newly chosen cut.
            continue
        try:
            before_sequence = str(before.extract(oriented.seq)).upper()
            after_sequence = str(after.extract(normalized.seq)).upper()
        except (KeyError, ValueError) as exc:
            raise OriginNormalizationError(
                f"could not validate GenBank feature {index + 1} ({before.type}): {exc}"
            ) from exc
        if before_sequence != after_sequence:
            raise OriginNormalizationError(
                f"GenBank feature {index + 1} ({before.type}) changed sequence during normalization"
            )
    return normalized


def normalize_plasmid_pair(
    fasta_sequence: str,
    genbank_record: SeqRecord,
    *,
    motif: str = DEFAULT_ORIGIN_MOTIF,
    sequence_name: str | None = None,
    fasta_label: str = "FASTA",
    genbank_label: str = "GenBank",
) -> NormalizedPlasmidPair:
    """Choose a shared canonical sequence and transform both input records."""
    fasta_candidates = _candidate_map(fasta_sequence, motif, fasta_label)
    genbank_candidates = _candidate_map(str(genbank_record.seq), motif, genbank_label)
    shared = fasta_candidates.keys() & genbank_candidates.keys()
    if not shared:
        raise OriginNormalizationError(
            "FASTA and GenBank sequences do not match after origin normalization"
        )
    if len(shared) != 1:
        raise OriginNormalizationError(
            f"origin motif {motif.upper()!r} gives multiple shared canonical starts for "
            f"{fasta_label} and {genbank_label}; refusing to choose one"
        )

    canonical = shared.pop()
    fasta_transform = _preferred_transform(fasta_candidates[canonical])
    genbank_transform = _preferred_transform(genbank_candidates[canonical])
    normalized_genbank = normalize_genbank_record(
        genbank_record,
        genbank_transform,
        sequence_name=sequence_name,
    )
    if str(normalized_genbank.seq).upper() != canonical:
        raise OriginNormalizationError(
            "internal error: normalized GenBank sequence does not equal the canonical FASTA"
        )
    return NormalizedPlasmidPair(
        sequence=canonical,
        fasta_transform=fasta_transform,
        genbank_transform=genbank_transform,
        genbank_record=normalized_genbank,
    )


def _read_one_record(path: Path, file_format: str) -> SeqRecord:
    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        with warnings.catch_warnings():
            # Several EPI2ME GenBank exports have nonstandard LOCUS spacing.
            # Biopython recovers the correct name/length; all important fields
            # are validated again before the normalized files are published.
            warnings.filterwarnings(
                "ignore",
                message="Attempting to parse malformed locus line",
                category=BiopythonParserWarning,
            )
            records = list(SeqIO.parse(handle, file_format))
    if len(records) != 1:
        raise OriginNormalizationError(
            f"expected exactly one {file_format} record in {path}, found {len(records)}"
        )
    return records[0]


def _temporary_output_path(destination: Path) -> Path:
    with tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent, delete=False
    ) as handle:
        return Path(handle.name)


def normalize_plasmid_files(
    src_fasta: Path,
    src_genbank: Path,
    dst_fasta: Path,
    dst_genbank: Path,
    sequence_name: str,
    *,
    motif: str = DEFAULT_ORIGIN_MOTIF,
) -> dict[str, object]:
    """Normalize one FASTA/GenBank pair and atomically publish validated files."""
    src_fasta = Path(src_fasta)
    src_genbank = Path(src_genbank)
    dst_fasta = Path(dst_fasta)
    dst_genbank = Path(dst_genbank)
    if {src_fasta.resolve(), src_genbank.resolve()} & {
        dst_fasta.resolve(), dst_genbank.resolve()
    }:
        raise OriginNormalizationError(
            "normalization destinations must not overwrite the original FASTA or GenBank inputs"
        )
    if dst_fasta.resolve() == dst_genbank.resolve():
        raise OriginNormalizationError("FASTA and GenBank destinations must be different files")
    fasta_record = _read_one_record(src_fasta, "fasta")
    genbank_record = _read_one_record(src_genbank, "genbank")
    normalized = normalize_plasmid_pair(
        str(fasta_record.seq),
        genbank_record,
        motif=motif,
        sequence_name=sequence_name,
        fasta_label=str(src_fasta),
        genbank_label=str(src_genbank),
    )

    dst_fasta.parent.mkdir(parents=True, exist_ok=True)
    dst_genbank.parent.mkdir(parents=True, exist_ok=True)
    temp_fasta = _temporary_output_path(dst_fasta)
    temp_genbank = _temporary_output_path(dst_genbank)
    try:
        output_fasta_record = SeqRecord(
            Seq(normalized.sequence),
            id=sequence_name,
            name=sequence_name,
            description="",
        )
        with temp_fasta.open("w", encoding="ascii") as handle:
            FastaWriter(handle, wrap=80).write_file([output_fasta_record])
        SeqIO.write(normalized.genbank_record, temp_genbank, "genbank")

        checked_fasta = _read_one_record(temp_fasta, "fasta")
        checked_genbank = _read_one_record(temp_genbank, "genbank")
        checked_fasta_sequence = str(checked_fasta.seq).upper()
        checked_genbank_sequence = str(checked_genbank.seq).upper()
        if checked_fasta_sequence != normalized.sequence or checked_genbank_sequence != normalized.sequence:
            raise OriginNormalizationError("written FASTA/GenBank sequences failed validation")
        if len(checked_genbank.features) != len(genbank_record.features):
            raise OriginNormalizationError("written GenBank failed feature-count validation")
        for index, (expected, checked) in enumerate(
            zip(normalized.genbank_record.features, checked_genbank.features)
        ):
            if expected.type != checked.type or expected.qualifiers != checked.qualifiers:
                raise OriginNormalizationError(
                    f"written GenBank feature {index + 1} failed metadata validation"
                )
            if expected.location != checked.location:
                raise OriginNormalizationError(
                    f"written GenBank feature {index + 1} failed location validation"
                )
        os.replace(temp_fasta, dst_fasta)
        os.replace(temp_genbank, dst_genbank)
    finally:
        for temporary_path in (temp_fasta, temp_genbank):
            if temporary_path.exists():
                temporary_path.unlink()

    audit = {
        "motif": motif.upper(),
        "canonical_sequence_sha256": hashlib.sha256(normalized.sequence.encode("ascii")).hexdigest(),
        "fasta": normalized.fasta_transform.as_dict(),
        "genbank": normalized.genbank_transform.as_dict(),
    }
    return {
        "name": sequence_name,
        "sequence": normalized.sequence,
        "length_bp": len(normalized.sequence),
        "origin_normalization": audit,
    }
