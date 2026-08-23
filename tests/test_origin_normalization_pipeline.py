import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

import epi2me_to_final_package as pipeline
from plasmid_normalization import DEFAULT_ORIGIN_MOTIF


def test_package_normalizes_reference_before_first_alignment(tmp_path: Path, monkeypatch):
    motif = DEFAULT_ORIGIN_MOTIF
    canonical = motif + "ACGCGTTA"
    input_sequence = "GG" + canonical
    fasta_path = tmp_path / "barcode01.final.fasta"
    gbk_path = tmp_path / "barcode01.annotations.gbk"
    consensus_fastq = tmp_path / "barcode01.final.fastq"
    raw_fastq = tmp_path / "barcode01.fastq"

    SeqIO.write(SeqRecord(Seq(input_sequence), id="old", description=""), fasta_path, "fasta")
    genbank = SeqRecord(Seq(input_sequence), id="old", name="old", description="test")
    genbank.annotations.update({"molecule_type": "DNA", "topology": "circular"})
    genbank.features = [
        SeqFeature(SimpleLocation(0, len(input_sequence)), type="source"),
        SeqFeature(SimpleLocation(1, 5, strand=1), type="misc_feature"),
    ]
    SeqIO.write(genbank, gbk_path, "genbank")
    fastq_text = f"@consensus\n{input_sequence}\n+\n{'I' * len(input_sequence)}\n"
    consensus_fastq.write_text(fastq_text, encoding="ascii")
    raw_fastq.write_text(fastq_text, encoding="ascii")

    observed = {}

    def fake_alignment(raw_reads, reference_fasta, output_dir, **kwargs):
        normalized = SeqIO.read(reference_fasta, "fasta")
        observed["alignment_reference"] = str(normalized.seq)
        return {"sorted_bam": str(Path(output_dir) / "aligned.sorted.bam")}

    def fake_report(**kwargs):
        report_dir = Path(kwargs["out_dir"])
        report_dir.mkdir(parents=True, exist_ok=True)
        normalized_fasta = SeqIO.read(kwargs["contig_fasta"], "fasta")
        normalized_genbank = SeqIO.read(kwargs["gbk_path"], "genbank")
        observed["report_fasta"] = str(normalized_fasta.seq)
        observed["report_genbank"] = str(normalized_genbank.seq)
        per_base = report_dir / "per_base_details.csv"
        low_confidence = report_dir / "low_confidence_bases.csv"
        per_base.write_text("position,base,depth,qscore\n", encoding="ascii")
        low_confidence.write_text("position,base,depth,qscore\n", encoding="ascii")
        return {
            "genbank_summary": {"length_bp": len(normalized_genbank)},
            "maf_summary": None,
            "outputs": {
                "per_base_details_csv": str(per_base),
                "low_confidence_bases_csv": str(low_confidence),
                "feature_map_png": None,
                "report_summary_json": str(report_dir / "report_summary.json"),
            },
        }

    def fake_plot(*args, **kwargs):
        path = Path(args[-1])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"plot")
        return path

    def fake_pdf(**kwargs):
        output = Path(kwargs["output_pdf"])
        output.write_bytes(b"%PDF-test")

    monkeypatch.setattr(pipeline, "run_fastq_pipeline", fake_alignment)
    monkeypatch.setattr(pipeline, "generate_report_data", fake_report)
    monkeypatch.setattr(pipeline, "plot_pdf_coverage_map", fake_plot)
    monkeypatch.setattr(pipeline, "plot_read_length_vs_bases", fake_plot)
    monkeypatch.setattr(pipeline, "render_pdf_report", fake_pdf)
    monkeypatch.setattr(pipeline, "DEFAULT_ECOLI_REFERENCE_FASTA", tmp_path / "missing-ecoli.fa")

    output_root = tmp_path / "output"
    result = pipeline.package_sample(
        "barcode01",
        {
            "barcode": "barcode01",
            "fasta": fasta_path,
            "gbk": gbk_path,
            "fastq": consensus_fastq,
            "raw_fastq": raw_fastq,
        },
        {
            "sample_name": "sample",
            "order_number": "TEST",
            "order_sample_number": "1",
        },
        output_root,
        logos=[],
        circular_coverage=False,
    )

    assert observed == {
        "alignment_reference": canonical + "GG",
        "report_fasta": canonical + "GG",
        "report_genbank": canonical + "GG",
    }
    assert result["origin_normalization"]["fasta"]["rotation_offset_bp"] == 2
    assert result["origin_normalization"]["consensus_fastq_quality"]["source"] == (
        "normalized_consensus_fastq"
    )

    summary_path = output_root / "_work" / result["sample_name"] / "package_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["origin_normalization"] == result["origin_normalization"]
    assert Path(summary["paths"]["fasta"]).read_text(encoding="ascii").splitlines()[1].startswith(motif)
    report_summary_path = summary_path.parent / "report" / "report_summary.json"
    report_summary = json.loads(report_summary_path.read_text(encoding="utf-8"))
    assert report_summary["origin_normalization"] == result["origin_normalization"]
