from pathlib import Path

from matplotlib.axes import Axes

from generate_report import parse_genbank_summary, plot_feature_map


GENBANK_WITH_COMPOUND_FEATURES = """\
LOCUS       TESTREC                  100 bp    DNA     circular SYN 01-JAN-2000
DEFINITION  synthetic test.
ACCESSION   TESTREC
VERSION     TESTREC.1
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            other sequences.
FEATURES             Location/Qualifiers
     source          1..100
                     /organism="synthetic construct"
     misc_feature    join(91..100,1..10)
                     /label="origin-spanning"
                     /note="first part of a long
                     note continued across lines"
     CDS             complement(join(40..50,70..80))
                     /label="reverse feature"
                     /pseudo
ORIGIN
        1 aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa
       61 aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa
//
"""


def test_parse_genbank_summary_handles_compound_locations_and_multiline_qualifiers(
    tmp_path: Path,
):
    gbk_path = tmp_path / "compound.gbk"
    gbk_path.write_text(GENBANK_WITH_COMPOUND_FEATURES, encoding="utf-8")

    summary = parse_genbank_summary(gbk_path)

    assert summary["locus_name"] == "TESTREC"
    assert summary["length_bp"] == 100
    assert summary["is_circular"] is True
    assert summary["feature_count"] == 3
    assert summary["feature_type_counts"] == {"CDS": 1, "misc_feature": 1, "source": 1}
    assert summary["labels"] == ["origin-spanning", "reverse feature"]

    origin_feature = summary["features"][1]
    assert origin_feature["location"] == "join(91..100,1..10)"
    assert origin_feature["segments"] == [
        {"start": 91, "end": 100, "strand": 1},
        {"start": 1, "end": 10, "strand": 1},
    ]
    assert origin_feature["qualifiers"]["note"] == (
        "first part of a long note continued across lines"
    )

    reverse_feature = summary["features"][2]
    assert reverse_feature["location"] == "complement(join(40..50,70..80))"
    assert reverse_feature["segments"] == [
        {"start": 70, "end": 80, "strand": -1},
        {"start": 40, "end": 50, "strand": -1},
    ]
    assert reverse_feature["qualifiers"]["pseudo"] is True


def test_plot_feature_map_uses_structured_segments(tmp_path: Path, monkeypatch):
    plotted_ranges = []
    original_broken_barh = Axes.broken_barh

    def capture_broken_barh(self, xranges, yrange, **kwargs):
        plotted_ranges.extend(xranges)
        return original_broken_barh(self, xranges, yrange, **kwargs)

    monkeypatch.setattr(Axes, "broken_barh", capture_broken_barh)
    output_path = tmp_path / "feature_map.png"
    result = plot_feature_map(
        {
            "features": [
                {
                    "type": "misc_feature",
                    "location": "not-a-parseable-location",
                    "qualifiers": {"label": "origin-spanning"},
                    "segments": [
                        {"start": 91, "end": 100, "strand": 1},
                        {"start": 1, "end": 10, "strand": 1},
                    ],
                }
            ]
        },
        contig_length=100,
        out_path=output_path,
        title="Annotation Map",
    )

    assert result == output_path
    assert output_path.exists()
    assert plotted_ranges == [(90, 10), (0, 10)]
