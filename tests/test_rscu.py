from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

from bioseeker.analysis.rscu_analysis import (
    calculate_rscu_from_counts,
    count_codons_for_sequence,
    global_codon_counts,
    global_rscu,
)
from bioseeker.utils.genetic_code import get_genetic_code


def test_genetic_code_excludes_stops():
    code = get_genetic_code(1)
    all_codons = set(str(c) for codons in code.synonymous_codons().values() for c in codons)
    stop_codons = {"TAA", "TAG", "TGA"}
    assert not any(stop in all_codons for stop in stop_codons), "Stop codons should be excluded"


def test_count_codons():
    seq_aligned = "TTT-TTC-TTC-TAA"
    counts = count_codons_for_sequence(seq_aligned)
    assert counts["TTT"] == 1
    assert counts["TTC"] == 2
    assert counts["TAA"] == 1


def test_calculate_rscu_values():
    counts = Counter({"TTT": 2, "TTC": 0})
    rscu = calculate_rscu_from_counts(counts)
    assert rscu["TTT"] == 2.0
    assert rscu["TTC"] == 0.0


def test_global_codon_counts():
    with TemporaryDirectory() as tmpdir:
        input_dir = Path(tmpdir)
        fasta_content = """>SpeciesA
TTTTTT
>SpeciesB
TTCTTC
"""
        fasta_file = input_dir / "test.afa"
        fasta_file.write_text(fasta_content)

        counts = global_codon_counts([fasta_file])
        assert "SpeciesA" in counts
        assert counts["SpeciesA"]["TTT"] == 2
        assert counts["SpeciesB"]["TTC"] == 2


def test_global_rscu_aggregation():
    with TemporaryDirectory() as tmpdir:
        input_dir = Path(tmpdir) / "input"
        input_dir.mkdir()

        fasta1 = """>SpeciesA
TTT
>SpeciesB
TTC
"""
        (input_dir / "1.afa").write_text(fasta1)

        fasta2 = """>SpeciesA
TTT
>SpeciesB
TTC
"""
        (input_dir / "2.afa").write_text(fasta2)

        df = global_rscu([input_dir / "1.afa", input_dir / "2.afa"])

        assert "SpeciesA" in df.columns
        assert "SpeciesB" in df.columns

        # RSCU(TTT) for SpeciesA should be 2.0 (since TTC is 0)
        assert df.loc["TTT", "SpeciesA"] == 2.0
        assert df.loc["TTC", "SpeciesB"] == 2.0
