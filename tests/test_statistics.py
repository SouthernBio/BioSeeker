import pandas as pd
import pytest
from bioseeker.analysis.statistics import (
    calculate_expected_bicodon_rates,
    calculate_z_scores,
    filter_inhibitory_pairs
)


def test_calculate_expected_bicodon_rates():
    # Mock data
    # Codon A: Rate 0.5
    # Codon B: Rate 0.8
    codon_data = {
        'codon': ['AAA', 'TTT'],
        'ConservationRate': [0.5, 0.8],
        'ReferenceCount': [100, 100],
        'ConservationCount': [50, 80]
    }
    codon_df = pd.DataFrame(codon_data)
    
    # Bicodon AAATTT (A-B) -> Exp = 0.5 * 0.8 = 0.4
    bicodon_data = {
        'codon_pair': ['AAATTT'],
        'ReferenceCount': [100],
        'ConservationCount': [40] # 0.4 observed
    }
    bicodon_df = pd.DataFrame(bicodon_data)
    
    result = calculate_expected_bicodon_rates(codon_df, bicodon_df)
    assert result.iloc[0]['ExpectedRate'] == pytest.approx(0.4)


def test_calculate_z_scores():
    # Case: Observed >> Expected
    # N=100, Exp=0.1, Obs=0.5
    # Mean = 0.1, Var = 0.1*0.9/100 = 0.0009, Sigma = 0.03
    # Z = (0.5 - 0.1) / 0.03 = 13.33 -> Conserved
    df = pd.DataFrame({
        'codon_pair': ['AAAAAA'],
        'ReferenceCount': [100],
        'ExpectedRate': [0.1],
        'ConservationRate': [0.5]
    })
    
    res = calculate_z_scores(df)
    assert res.iloc[0]['ZScore'] > 3
    assert res.iloc[0]['IsConserved'] == True

def test_filter_inhibitory_pairs():
    # Rare: AAA
    # Pair 1: AAA-AAA (Conserved) -> Inhibitory
    # Pair 2: AAA-TTT (Conserved) -> Not Inhibitory (one rare)
    # Pair 3: AAA-AAA (Not Conserved) -> Not Inhibitory
    
    rare = ['AAA']
    df = pd.DataFrame({
        'codon_pair': ['AAAAAA', 'AAATTT', 'AAAAAA'],
        'IsConserved': [True, True, False]
    })
    
    res = filter_inhibitory_pairs(df, rare)
    assert len(res) == 1
    assert res.iloc[0]['codon_pair'] == 'AAAAAA'
    assert res.iloc[0]['IsConserved'] == True
