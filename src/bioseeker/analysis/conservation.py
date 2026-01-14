import numpy as np
import pandas as pd
from typing import List, Tuple
from ..utils.genetic_code import CODON_TUPLE, CODON_PAIRS_TUPLE


def split_into_codons(sequence: str, orf: int) -> List[str]:
    """Splits a sequence into codons based on the ORF."""
    return [sequence[i:i+3] for i in range(orf, len(sequence) - 2, 3)]


def split_into_bicodons(sequence: str, orf: int) -> List[str]:
    """Splits a sequence into bicodons (pairs of codons) based on the ORF."""
    return [sequence[i:i+6] for i in range(orf, len(sequence) - 5, 3)]


def calculate_conservation(sequences: List[str], orf: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculates codon and codon pair conservation rates using the first sequence as reference.
    
    Args:
        sequences (List[str]): List of DNA sequences (aligned).
        orf (int): Reading frame (0, 1, or 2).

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: (codon_df, bicodon_df)
    """
    if orf not in {0, 1, 2}:
        raise ValueError("ORF must be 0, 1, or 2")
    
    if not sequences:
        raise ValueError("No sequences provided")

    # Prepare data structures
    codons_list = [split_into_codons(s, orf) for s in sequences]
    bicodons_list = [split_into_bicodons(s, orf) for s in sequences]

    # Convert to numpy arrays (strings)
    # Assuming all sequences are aligned (same length)
    try:
        codons_arr = np.array(codons_list)
        bicodons_arr = np.array(bicodons_list)
    except ValueError as e:
        # Fallback or error if jagged
        raise ValueError("Sequences must be aligned (same length) for conservation analysis.") from e

    # Reference is the first sequence
    ref_codons = codons_arr[0]
    ref_bicodons = bicodons_arr[0]

    # Calculate conservation: Count matches in column
    # Broadcasting: (N_seqs, LEN) == (LEN,) -> Compare each row to ref
    # We want to compare column-wise to ref[col_idx].
    # codons_arr == ref_codons (Broadcasting effectively matches column-wise)
    
    codon_matches = (codons_arr == ref_codons)
    codon_conserved_counts = codon_matches.sum(axis=0) # Sum across rows (sequences) for each position

    bicodon_matches = (bicodons_arr == ref_bicodons)
    bicodon_conserved_counts = bicodon_matches.sum(axis=0)

    # Threshold for "Conserved" status (original code used > 0.9 ratio)
    n_seqs = len(sequences)
    codon_is_conserved = (codon_conserved_counts / n_seqs) > 0.9
    bicodon_is_conserved = (bicodon_conserved_counts / n_seqs) > 0.9

    # Now we need to aggregate this back to global counts for each Codon Type (TTG, GGC, etc.)
    # We need:
    # - ReferenceCount: How many times each unique codon appears in the REFERENCE sequence.
    # - ConservationCount: How many times that appearance was "conserved" (>90% match in others).

    # Initialize result dicts
    codon_stats = {c: {'ReferenceCount': 0, 'ConservationCount': 0} for c in CODON_TUPLE}
    bicodon_stats = {c: {'ReferenceCount': 0, 'ConservationCount': 0} for c in CODON_PAIRS_TUPLE}

    # Aggregate Codons
    # We zip the reference codons with their conservation status
    for i, codon in enumerate(ref_codons):
        if codon in codon_stats:
            codon_stats[codon]['ReferenceCount'] += 1
            if codon_is_conserved[i]:
                codon_stats[codon]['ConservationCount'] += 1
    
    # Aggregate Bicodons
    for i, bicodon in enumerate(ref_bicodons):
        if bicodon in bicodon_stats:
            bicodon_stats[bicodon]['ReferenceCount'] += 1
            if bicodon_is_conserved[i]:
                bicodon_stats[bicodon]['ConservationCount'] += 1

    # Create DataFrames
    codon_df = pd.DataFrame.from_dict(codon_stats, orient='index')
    codon_df.index.name = 'codon'
    codon_df.reset_index(inplace=True)
    
    bicodon_df = pd.DataFrame.from_dict(bicodon_stats, orient='index')
    bicodon_df.index.name = 'codon_pair'
    bicodon_df.reset_index(inplace=True)

    return codon_df, bicodon_df
