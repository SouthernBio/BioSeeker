import pandas as pd
import numpy as np
from typing import Tuple, List


def calculate_expected_bicodon_rates(codon_df: pd.DataFrame, bicodon_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates expected bicodon conservation rates based on constituent codons.
    Rate_exp = Rate(Codon1) * Rate(Codon2)
    """
    # Create a mapping of codon -> rate
    # Handle missing or zero reference counts
    # Rate = ConservationCount / ReferenceCount
    if 'ConservationRate' not in codon_df.columns:
        # Calculate if missing (though assembler should have done it)
        codon_df['ConservationRate'] = codon_df.apply(
            lambda row: row['ConservationCount'] / row['ReferenceCount'] if row['ReferenceCount'] > 0 else 0, axis=1
        )
    
    rate_map = dict(zip(codon_df['codon'], codon_df['ConservationRate']))
    
    # Process bicodon dataframe
    df = bicodon_df.copy()
    if 'ConservationRate' not in df.columns:
        df['ConservationRate'] = df.apply(
             lambda row: row['ConservationCount'] / row['ReferenceCount'] if row['ReferenceCount'] > 0 else 0, axis=1
        )

    def get_expected(bicodon):
        c1 = bicodon[:3]
        c2 = bicodon[3:]
        r1 = rate_map.get(c1, 0)
        r2 = rate_map.get(c2, 0)
        return r1 * r2

    df['ExpectedRate'] = df['codon_pair'].apply(get_expected)
    return df


def perform_regression(bicodon_df: pd.DataFrame) -> Tuple[float, float]:
    """
    Performs linear regression (Observed vs Expected) forcing intercept to 0.
    Returns (slope, r_squared).
    """
    # Filter out cases with 0 expected rate to avoid noise? Or keep all?
    # Usually keep all valid points.
    
    y = bicodon_df['ConservationRate'].values
    x = bicodon_df['ExpectedRate'].values
    
    # Linear regression with intercept=0: y = m*x
    # Least squares: m = sum(x*y) / sum(x^2)
    
    x_new = x[:, np.newaxis]
    m, _, _, _ = np.linalg.lstsq(x_new, y, rcond=None)
    slope = m[0]
    
    # Calculate R^2
    y_pred = slope * x
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    
    return slope, r_squared


def calculate_z_scores(bicodon_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates Z-scores for bicodon conservation.
    Z = (Observed - Expected) / sqrt(Expected * (1-Expected) / N)
    """
    df = bicodon_df.copy()
    
    def get_z(row):
        n = row['ReferenceCount']
        if n == 0:
            return 0
        p_obs = row['ConservationRate']
        p_exp = row['ExpectedRate']
        
        # Avoid division by zero if p_exp is 0 or 1
        if p_exp <= 0 or p_exp >= 1:
            return 0 # logic choice: if expected is 0 and observed is >0, it's infinite, but let's clamp or handle
        
        sigma = np.sqrt(p_exp * (1 - p_exp) / n)
        if sigma == 0:
            return 0
            
        return (p_obs - p_exp) / sigma

    df['ZScore'] = df.apply(get_z, axis=1)
    df['IsConserved'] = df['ZScore'] > 3
    return df


def filter_inhibitory_pairs(bicodon_df: pd.DataFrame, rare_codons: List[str]) -> pd.DataFrame:
    """
    Filters conserved bicodons where both codons are rare.
    """
    rare_set = set(rare_codons)
    
    def is_inhibitory(row):
        if not row['IsConserved']:
            return False
        bicodon = row['codon_pair']
        c1 = bicodon[:3]
        c2 = bicodon[3:]
        return (c1 in rare_set) and (c2 in rare_set)
    
    df = bicodon_df.copy()
    df['IsInhibitory'] = df.apply(is_inhibitory, axis=1)
    return df[df['IsInhibitory']]
