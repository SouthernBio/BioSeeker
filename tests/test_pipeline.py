from pathlib import Path

import pandas as pd

from bioseeker.analysis.assembler import ConservationAggregator
from bioseeker.analysis.conservation import conservation_for_records
from bioseeker.cli import main

DATA_DIR = Path(__file__).parent / "data"


def test_aggregator_sums_multiple_results():
    r1 = conservation_for_records([("A", "ATGGCC"), ("B", "ATGGCC")])
    r2 = conservation_for_records([("A", "ATGGCA"), ("B", "ATGTCA")])
    agg = ConservationAggregator()
    agg.add(r1)
    agg.add(r2)
    c_df, b_df = agg.get_results()

    atg = c_df[c_df["Codon"] == "ATG"].iloc[0]
    assert atg["ReferenceCount"] == 2
    assert atg["ConservationCount"] == 2
    assert atg["ConservationRate"] == 1.0

    gca = c_df[c_df["Codon"] == "GCA"].iloc[0]
    assert gca["ReferenceCount"] == 1
    assert gca["ConservationCount"] == 0

    pair = b_df[b_df["CodonPair"] == "ATGGCC"].iloc[0]
    assert pair["ConservationCount"] == 1
    assert agg.n_alignments == 2
    assert agg.n_reference_codons == 4


def test_analyze_end_to_end(tmp_path):
    rc = main(
        ["-q", "analyze", "-i", str(DATA_DIR), "-o", str(tmp_path), "--orfs", "0", "--no-plot"]
    )
    assert rc == 0
    stats = pd.read_csv(tmp_path / "bicodon_stats_ORF0.csv")
    assert len(stats) == 3721
    for col in ("ExpectedRate", "SigmaDistance", "Outlier", "NormalizedScore", "DipeptideZ"):
        assert col in stats.columns
