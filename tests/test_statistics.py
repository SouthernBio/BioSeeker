import pandas as pd
import pytest

from bioseeker.analysis.statistics import (
    analyze_pairs,
    calculate_expected_bicodon_rates,
    filter_inhibitory_pairs,
)


def test_calculate_expected_bicodon_rates():
    codon_data = {
        "Codon": ["AAA", "TTT"],
        "ConservationRate": [0.5, 0.8],
        "ReferenceCount": [100, 100],
        "ConservationCount": [50, 80],
    }
    codon_df = pd.DataFrame(codon_data)

    bicodon_data = {
        "CodonPair": ["AAATTT"],
        "Codon1": ["AAA"],
        "Codon2": ["TTT"],
        "ReferenceCount": [100],
        "ConservationCount": [40],
    }
    bicodon_df = pd.DataFrame(bicodon_data)

    result = calculate_expected_bicodon_rates(codon_df, bicodon_df)
    assert result.iloc[0]["ExpectedRate"] == pytest.approx(0.4)


def test_analyze_pairs_regression():
    # Make a dummy dataset where Expected == Observed
    codon_df = pd.DataFrame(
        {
            "Codon": ["AAA", "TTT", "CCC"],
            "ReferenceCount": [100, 100, 100],
            "ConservationCount": [50, 80, 50],
            "ConservationRate": [0.5, 0.8, 0.5],
        }
    )

    bicodon_df = pd.DataFrame(
        {
            "CodonPair": ["AAATTT", "TTTCCC"],
            "Codon1": ["AAA", "TTT"],
            "Codon2": ["TTT", "CCC"],
            "ReferenceCount": [100, 100],
            "ConservationCount": [40, 40],  # exactly Expected
            "ConservationRate": [0.4, 0.4],
        }
    )

    stats_df, fit = analyze_pairs(codon_df, bicodon_df, n_reference_codons=300, sigma_threshold=3.0)

    assert fit.slope == pytest.approx(1.0)
    assert fit.r_squared == pytest.approx(1.0)
    assert stats_df.iloc[0]["Outlier"] == ""


def test_filter_inhibitory_pairs():
    rare = ["AAA"]
    df = pd.DataFrame(
        {
            "CodonPair": ["AAAAAA", "AAATTT", "AAAAAA"],
            "Codon1": ["AAA", "AAA", "AAA"],
            "Codon2": ["AAA", "TTT", "AAA"],
            "Outlier": ["low", "low", "high"],
        }
    )

    res = filter_inhibitory_pairs(df, rare)
    assert len(res) == 1
    assert res.iloc[0]["CodonPair"] == "AAAAAA"
    assert res.iloc[0]["Outlier"] == "high"
