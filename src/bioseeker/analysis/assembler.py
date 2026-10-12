"""Aggregation of conservation results over many alignment files."""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
import pandas as pd

from .conservation import ConservationResult


def add_rates(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of ``df`` with ``ConservationRate = ConservationCount / ReferenceCount``.

    The rate is ``NaN`` when ``ReferenceCount`` is 0 (undefined).
    """
    out = df.copy()
    ref = out["ReferenceCount"].to_numpy(dtype=float)
    cons = out["ConservationCount"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        out["ConservationRate"] = np.where(ref > 0, cons / np.where(ref > 0, ref, 1), np.nan)
    return out


class ConservationAggregator:
    """Aggregates :class:`ConservationResult` objects from multiple files.

    All results must have been produced with the same genetic code, so that
    their rows are in the same order.
    """

    def __init__(self) -> None:
        self._codons: Optional[pd.DataFrame] = None
        self._pairs: Optional[pd.DataFrame] = None
        self.n_reference_codons: int = 0
        self.n_alignments: int = 0

    def add(self, result: ConservationResult) -> None:
        """Add a single result to the running totals."""
        if self._codons is None:
            self._codons = result.codons.copy()
            self._pairs = result.pairs.copy()
        else:
            if not (
                self._codons["Codon"].equals(result.codons["Codon"])
                and self._pairs["CodonPair"].equals(result.pairs["CodonPair"])
            ):
                raise ValueError("Cannot aggregate results computed with different genetic codes")
            for total, new in ((self._codons, result.codons), (self._pairs, result.pairs)):
                for col in ("ReferenceCount", "ConservationCount"):
                    total[col] = total[col].to_numpy() + new[col].to_numpy()
        self.n_reference_codons += result.n_reference_codons
        self.n_alignments += 1

    @property
    def empty(self) -> bool:
        return self._codons is None

    def get_results(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Return ``(codon_df, pair_df)`` with aggregated counts and conservation rates.

        Returns two empty DataFrames if nothing was added.
        """
        if self._codons is None:
            return pd.DataFrame(), pd.DataFrame()
        return add_rates(self._codons), add_rates(self._pairs)
