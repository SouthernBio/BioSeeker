"""Codon pair statistics."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable, Tuple

import numpy as np
import pandas as pd

from ..utils.genetic_code import get_genetic_code
from .assembler import add_rates


@dataclass
class RegressionResult:
    """Result of a least-squares fit through the origin."""

    slope: float
    r_squared: float
    sigma: float
    n: int

    def as_dict(self) -> dict:
        return asdict(self)


def _rate_column(df: pd.DataFrame) -> pd.DataFrame:
    return df if "ConservationRate" in df.columns else add_rates(df)


def calculate_expected_bicodon_rates(
    codon_df: pd.DataFrame, bicodon_df: pd.DataFrame
) -> pd.DataFrame:
    """Return a copy of ``bicodon_df`` with ``Codon1``, ``Codon2`` and ``ExpectedRate`` columns.

    ``ExpectedRate`` is the product of the conservation rates of both codons.
    Neither input is modified.
    """
    codons = _rate_column(codon_df)
    rate_map = pd.Series(
        codons["ConservationRate"].to_numpy(dtype=float), index=codons["Codon"].to_numpy()
    )

    df = _rate_column(bicodon_df).copy()
    pairs = df["CodonPair"].astype(str)
    df["Codon1"] = pairs.str[:3]
    df["Codon2"] = pairs.str[3:]
    r1 = df["Codon1"].map(rate_map).to_numpy(dtype=float)
    r2 = df["Codon2"].map(rate_map).to_numpy(dtype=float)
    df["ExpectedRate"] = r1 * r2
    return df


def _fit_mask(df: pd.DataFrame, x: str, y: str) -> np.ndarray:
    mask = np.isfinite(df[x].to_numpy(dtype=float)) & np.isfinite(df[y].to_numpy(dtype=float))
    if "ReferenceCount" in df.columns:
        mask &= df["ReferenceCount"].to_numpy() > 0
    return mask


def perform_regression(
    bicodon_df: pd.DataFrame, x: str = "ExpectedRate", y: str = "ConservationRate"
) -> RegressionResult:
    """Least-squares fit of ``y = m * x`` (intercept fixed at zero).

    Rows with ``ReferenceCount == 0`` or non-finite values are ignored.

    Returns:
        :class:`RegressionResult` with slope, uncentred R², residual standard
        error σ and the number of points used. Values are ``NaN`` when the fit
        is not possible (fewer than 2 points or all ``x`` equal to zero).
    """
    mask = _fit_mask(bicodon_df, x, y)
    xs = bicodon_df[x].to_numpy(dtype=float)[mask]
    ys = bicodon_df[y].to_numpy(dtype=float)[mask]
    n = int(mask.sum())

    sxx = float(np.dot(xs, xs))
    if n < 2 or sxx == 0.0:
        return RegressionResult(slope=float("nan"), r_squared=float("nan"), sigma=float("nan"), n=n)

    slope = float(np.dot(xs, ys) / sxx)
    ss_res = float(np.sum((ys - slope * xs) ** 2))
    syy = float(np.dot(ys, ys))
    r_squared = 1.0 - ss_res / syy if syy > 0 else float("nan")
    sigma = float(np.sqrt(ss_res / (n - 1)))
    return RegressionResult(slope=slope, r_squared=r_squared, sigma=sigma, n=n)


def flag_outliers(
    bicodon_df: pd.DataFrame,
    fit: RegressionResult,
    sigma_threshold: float = 3.0,
    x: str = "ExpectedRate",
    y: str = "ConservationRate",
) -> pd.DataFrame:
    """Return a copy with ``FittedRate``, ``Residual``, ``SigmaDistance`` and ``Outlier`` columns.

    ``Outlier`` is ``"high"`` (more conserved than predicted), ``"low"`` (less
    conserved than predicted) or ``""``. Rows excluded from the fit get NaN.
    """
    df = bicodon_df.copy()
    mask = _fit_mask(df, x, y)
    fitted = fit.slope * df[x].to_numpy(dtype=float)
    residual = df[y].to_numpy(dtype=float) - fitted
    fitted[~mask] = np.nan
    residual[~mask] = np.nan

    if np.isfinite(fit.sigma) and fit.sigma > 0:
        distance = residual / fit.sigma
    else:
        distance = np.full(len(df), np.nan)

    outlier = np.full(len(df), "", dtype=object)
    with np.errstate(invalid="ignore"):
        outlier[distance > sigma_threshold] = "high"
        outlier[distance < -sigma_threshold] = "low"

    df["FittedRate"] = fitted
    df["Residual"] = residual
    df["SigmaDistance"] = distance
    df["Outlier"] = outlier
    return df


def normalized_scores(
    bicodon_df: pd.DataFrame, codon_df: pd.DataFrame, pseudocount: float
) -> pd.DataFrame:
    """Return a copy with the ``NormalizedScore`` column.

    ``X = (CPCR + p) / ((CCR1 + p) * (CCR2 + p))``. Pairs never observed in the
    reference (``ReferenceCount == 0``) get ``NaN``.
    """
    if pseudocount < 0:
        raise ValueError("pseudocount must be >= 0")
    df = (
        bicodon_df
        if "Codon1" in bicodon_df.columns
        else calculate_expected_bicodon_rates(codon_df, bicodon_df)
    )
    df = df.copy()
    codons = _rate_column(codon_df)
    rate_map = pd.Series(
        codons["ConservationRate"].to_numpy(dtype=float), index=codons["Codon"].to_numpy()
    )
    r1 = df["Codon1"].map(rate_map).to_numpy(dtype=float)
    r2 = df["Codon2"].map(rate_map).to_numpy(dtype=float)
    cpcr = df["ConservationRate"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        score = (cpcr + pseudocount) / ((r1 + pseudocount) * (r2 + pseudocount))
    score[df["ReferenceCount"].to_numpy() == 0] = np.nan
    score[~np.isfinite(score)] = np.nan
    df["NormalizedScore"] = score
    return df


def dipeptide_z_scores(bicodon_df: pd.DataFrame, genetic_code: int = 1) -> pd.DataFrame:
    """Return a copy with ``Dipeptide`` and ``DipeptideZ`` columns.

    ``DipeptideZ`` is the z-score of ``log(NormalizedScore)`` among all codon
    pairs encoding the same dipeptide. Families with fewer than two valid
    members, or with zero variance, get ``NaN``.
    """
    code = get_genetic_code(genetic_code)
    df = bicodon_df.copy()
    df["Dipeptide"] = [code.dipeptide(p) for p in df["CodonPair"].astype(str)]
    score = df["NormalizedScore"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_score = np.where(score > 0, np.log(score), np.nan)
    log_s = pd.Series(log_score, index=df.index)
    grouped = log_s.groupby(df["Dipeptide"].to_numpy())
    mean = grouped.transform("mean")
    std = grouped.transform(lambda s: s.std(ddof=1))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = (log_s - mean) / std
    z = z.to_numpy(dtype=float, copy=True)
    z[~np.isfinite(z)] = np.nan
    df["DipeptideZ"] = z
    return df


def analyze_pairs(
    codon_df: pd.DataFrame,
    pair_df: pd.DataFrame,
    n_reference_codons: int,
    *,
    genetic_code: int = 1,
    sigma_threshold: float = 3.0,
) -> Tuple[pd.DataFrame, RegressionResult]:
    """Run the full codon-pair analysis.

    Returns:
        ``(stats_df, fit)`` where ``stats_df`` contains the expected rates,
        regression residuals, outlier flags, normalized scores and dipeptide
        z-scores for every codon pair.
    """
    df = calculate_expected_bicodon_rates(codon_df, pair_df)
    fit = perform_regression(df)
    df = flag_outliers(df, fit, sigma_threshold)
    pseudocount = 1.0 / n_reference_codons if n_reference_codons > 0 else 0.0
    df = normalized_scores(df, codon_df, pseudocount)
    df = dipeptide_z_scores(df, genetic_code)
    return df, fit


def filter_inhibitory_pairs(bicodon_df: pd.DataFrame, rare_codons: Iterable[str]) -> pd.DataFrame:
    """Return the pairs that are high outliers and made of two rare codons."""
    rare = {c.strip().upper() for c in rare_codons}
    pairs = bicodon_df["CodonPair"].astype(str)
    is_rare_pair = pairs.str[:3].isin(rare) & pairs.str[3:].isin(rare)
    is_high = bicodon_df["Outlier"].to_numpy() == "high"
    return bicodon_df[is_rare_pair.to_numpy() & is_high].reset_index(drop=True)


def select_pairs(bicodon_df: pd.DataFrame, pairs: Iterable[str]) -> pd.DataFrame:
    """Return the rows of ``bicodon_df`` for the given codon pairs (e.g. a list of known ICPs).

    Pairs may be written as ``CGACGA`` or ``CGA-CGA``; RNA ``U`` is converted to ``T``.
    """
    wanted = [p.strip().upper().replace("-", "").replace("U", "T") for p in pairs]
    wanted = [p for p in wanted if p]
    return bicodon_df[bicodon_df["CodonPair"].isin(wanted)].reset_index(drop=True)
