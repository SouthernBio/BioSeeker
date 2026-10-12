"""Relative Synonymous Codon Usage (RSCU) analysis."""

from __future__ import annotations

import os
from collections import Counter
from typing import Dict, Iterable, List

import pandas as pd

from ..io.reader import read_fasta
from ..utils.genetic_code import get_genetic_code, is_valid_codon


def count_codons_for_sequence(sequence: str) -> Counter[str]:
    """Count non-overlapping trinucleotides in a sequence.

    Gaps (``-``) are removed before counting. Ambiguous codons (containing
    characters other than ACGT) are ignored.
    """
    seq_clean = sequence.replace("-", "").upper()
    counts: Counter[str] = Counter()
    for i in range(0, len(seq_clean) - len(seq_clean) % 3, 3):
        codon = seq_clean[i : i + 3]
        if is_valid_codon(codon):
            counts[codon] += 1
    return counts


def calculate_rscu_from_counts(
    codon_counts: Counter[str], genetic_code: int = 1
) -> Dict[str, float]:
    """Calculate RSCU values from raw codon counts.

    Returns:
        Mapping from codon to RSCU value. Stop codons are excluded. If an amino
        acid is completely absent, all its synonymous codons get RSCU 0.
    """
    code = get_genetic_code(genetic_code)
    rscu_values: Dict[str, float] = {}

    for aa, syn_codons in code.synonymous_codons().items():
        total_aa_count = sum(codon_counts[c] for c in syn_codons)
        if total_aa_count == 0:
            for codon in syn_codons:
                rscu_values[codon] = 0.0
            continue

        expected = total_aa_count / len(syn_codons)
        for codon in syn_codons:
            rscu_values[codon] = codon_counts[codon] / expected

    return rscu_values


def global_codon_counts(files: Iterable[str | "os.PathLike[str]"]) -> Dict[str, Counter[str]]:
    """Aggregate codon counts per species across multiple alignment files.

    Args:
        files: Paths to FASTA alignment files.

    Returns:
        Mapping from species ID to a ``Counter`` of codon occurrences.
    """
    global_counts: Dict[str, Counter[str]] = {}
    for f in files:
        for record_id, seq in read_fasta(f):
            if record_id not in global_counts:
                global_counts[record_id] = Counter()
            global_counts[record_id].update(count_codons_for_sequence(seq))
    return global_counts


def global_rscu(files: Iterable[str | "os.PathLike[str]"], genetic_code: int = 1) -> pd.DataFrame:
    """Calculate the global RSCU per species across multiple alignment files.

    Returns:
        A DataFrame where the index is the codon (only sense codons) and the
        columns are the species IDs.
    """
    counts = global_codon_counts(files)
    code = get_genetic_code(genetic_code)
    species_list = sorted(counts.keys())

    data = {sp: calculate_rscu_from_counts(counts[sp], genetic_code) for sp in species_list}

    df = pd.DataFrame(data, index=code.sense_codons)
    df.index.name = "Codon"
    return df


def identify_rare_codons(
    rscu_df: pd.DataFrame, threshold: float = 1.0, genetic_code: int = 1
) -> List[str]:
    """Return sense codons where RSCU < threshold for all species.

    Amino acids encoded by a single codon (e.g. Met, Trp in the Standard code)
    are excluded since their RSCU is fixed at 1.0. Stop codons are also excluded.
    """
    code = get_genetic_code(genetic_code)
    # Filter out single-codon amino acids
    syn_fams = code.synonymous_codons()
    multi_codons = {c for fam in syn_fams.values() if len(fam) > 1 for c in fam}

    df = rscu_df[rscu_df.index.isin(multi_codons)]
    rare = df[(df < threshold).all(axis=1)].index.tolist()
    return sorted(rare)
