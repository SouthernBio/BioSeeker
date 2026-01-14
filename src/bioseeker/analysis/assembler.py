import pandas as pd
from typing import Tuple, Optional


class ConservationAggregator:
    """aggregates conservation results from multiple files."""
    
    def __init__(self):
        self.codon_total: Optional[pd.DataFrame] = None
        self.bicodon_total: Optional[pd.DataFrame] = None

    def add(self, codon_df: pd.DataFrame, bicodon_df: pd.DataFrame):
        """Adds a single result set to the total."""
        # Ensure indices for alignment
        c_df = codon_df.copy()
        if 'codon' in c_df.columns:
            c_df.set_index('codon', inplace=True)
            
        b_df = bicodon_df.copy()
        if 'codon_pair' in b_df.columns:
            b_df.set_index('codon_pair', inplace=True)

        if self.codon_total is None:
            self.codon_total = c_df
            self.bicodon_total = b_df
        else:
            self.codon_total = self.codon_total.add(c_df, fill_value=0)
            self.bicodon_total = self.bicodon_total.add(b_df, fill_value=0)

    def get_results(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Returns the aggregated dataframes with calculated rates."""
        if self.codon_total is None:
            return pd.DataFrame(), pd.DataFrame()

        # Calculate Rates
        # Avoid division by zero
        c_final = self.codon_total.copy()
        c_final['ConservationRate'] = c_final.apply(
            lambda row: row['ConservationCount'] / row['ReferenceCount'] if row['ReferenceCount'] > 0 else 0, axis=1
        )

        b_final = self.bicodon_total.copy()
        b_final['ConservationRate'] = b_final.apply(
            lambda row: row['ConservationCount'] / row['ReferenceCount'] if row['ReferenceCount'] > 0 else 0, axis=1
        )

        c_final.reset_index(inplace=True)
        b_final.reset_index(inplace=True)

        return c_final, b_final
