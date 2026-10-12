"""Position-specific conservation of codons and codon pairs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np
import pandas as pd

from ..utils.genetic_code import CODON_INDEX, get_genetic_code

# Lookup table: ASCII byte -> nucleotide index in TCAG order (or -1 if invalid).
_BASE_LUT = np.full(256, -1, dtype=np.int16)
for _i, _b in enumerate("TCAG"):
    _BASE_LUT[ord(_b)] = _i
    _BASE_LUT[ord(_b.lower())] = _i


class ReferenceNotFoundError(ValueError):
    """Raised when the requested reference sequence is not in an alignment."""


@dataclass
class ConservationResult:
    """Conservation counts for one alignment (or an aggregate of several).

    Attributes:
        codons: DataFrame with columns ``codon``, ``ReferenceCount``,
            ``ConservationCount`` (one row per sense codon).
        pairs: DataFrame with columns ``codon_pair``, ``ReferenceCount``,
            ``ConservationCount`` (one row per sense codon pair).
        n_reference_codons: Number of valid sense codons in the reference
            (used to derive pseudocounts).
    """

    codons: pd.DataFrame
    pairs: pd.DataFrame
    n_reference_codons: int


def split_into_codons(sequence: str, orf: int) -> List[str]:
    """Split a sequence into codons starting at offset ``orf``."""
    return [sequence[i : i + 3] for i in range(orf, len(sequence) - 2, 3)]


def split_into_bicodons(sequence: str, orf: int) -> List[str]:
    """Split a sequence into overlapping codon pairs starting at offset ``orf``."""
    return [sequence[i : i + 6] for i in range(orf, len(sequence) - 5, 3)]


def resolve_reference(ids: Sequence[str], reference_id: str | None) -> int:
    """Return the index of ``reference_id`` in ``ids`` (0 if ``reference_id`` is None).

    Raises:
        ReferenceNotFoundError: If ``reference_id`` is not present in ``ids``.
    """
    if reference_id is None:
        return 0
    try:
        return list(ids).index(reference_id)
    except ValueError:
        raise ReferenceNotFoundError(
            f"Reference '{reference_id}' not found (available: {', '.join(ids)})"
        ) from None


def remove_gap_only_columns(sequences: Sequence[str]) -> List[str]:
    """Drop alignment columns in which every sequence has a gap (``-`` or ``.``).

    Such columns carry no information about the sequences being compared. They
    appear when species are removed from a larger multiple alignment (e.g. with
    ``--species``), and if kept they would separate codons that are adjacent in
    every remaining sequence. In a codon-aware alignment gaps come in in-frame
    triplets, so removing them preserves the reading frame.
    """
    seqs = list(sequences)
    if not seqs or not seqs[0]:
        return seqs
    arr = np.frombuffer("".join(seqs).encode("ascii", "replace"), dtype=np.uint8).reshape(
        len(seqs), -1
    )
    keep = ~((arr == ord("-")) | (arr == ord("."))).all(axis=0)
    if keep.all():
        return seqs
    return [row.tobytes().decode("ascii") for row in arr[:, keep]]


def encode_codons(sequences: Sequence[str], orf: int) -> np.ndarray:
    """Encode aligned sequences as a matrix of codon indices.

    Returns:
        ``int16`` array of shape ``(n_sequences, n_codons)`` where each entry is
        the index of the codon in :data:`~bioseeker.utils.genetic_code.ALL_CODONS`
        (0-63), or ``-1`` if the codon contains a gap or a non-ACGT character.
    """
    length = len(sequences[0])
    n_codons = max((length - orf) // 3, 0)
    raw = np.frombuffer(
        "".join(s[orf : orf + 3 * n_codons] for s in sequences).encode("ascii", "replace"),
        dtype=np.uint8,
    ).reshape(len(sequences), n_codons, 3)
    bases = _BASE_LUT[raw]
    idx = bases[..., 0] * 16 + bases[..., 1] * 4 + bases[..., 2]
    idx[(bases < 0).any(axis=2)] = -1
    return idx.astype(np.int16)


def calculate_conservation(
    sequences: Sequence[str],
    orf: int = 0,
    *,
    reference: int = 0,
    threshold: float = 1.0,
    genetic_code: int = 1,
) -> ConservationResult:
    """Calculate codon and codon pair conservation for one alignment.

    Args:
        sequences: Aligned DNA sequences (all of the same length).
        orf: Reading frame offset (0, 1 or 2).
        reference: Index of the reference sequence within ``sequences``.
        threshold: Minimum fraction of sequences (reference included) that must
            share the reference codon for it to count as conserved (0 < t <= 1).
        genetic_code: NCBI translation table ID used to identify stop codons.

    Returns:
        A :class:`ConservationResult`.

    Raises:
        ValueError: On invalid arguments or unaligned sequences.
    """
    if orf not in (0, 1, 2):
        raise ValueError("ORF must be 0, 1, or 2")
    if not sequences:
        raise ValueError("No sequences provided")
    if not 0.0 < threshold <= 1.0:
        raise ValueError("threshold must be in the interval (0, 1]")
    if not 0 <= reference < len(sequences):
        raise ValueError(f"reference index {reference} out of range")
    lengths = {len(s) for s in sequences}
    if len(lengths) != 1:
        raise ValueError(
            "Sequences must be aligned (same length) for conservation analysis; "
            f"found lengths {sorted(lengths)}"
        )

    sequences = remove_gap_only_columns(sequences)

    code = get_genetic_code(genetic_code)
    sense_mask = np.zeros(64, dtype=bool)
    sense_idx = np.array([CODON_INDEX[c] for c in code.sense_codons], dtype=np.int64)
    sense_mask[sense_idx] = True

    idx = encode_codons(sequences, orf)
    n_seqs = idx.shape[0]
    ref = idx[reference]

    # Valid reference positions: ACGT-only and not a stop codon.
    ref_ok = (ref >= 0) & sense_mask[np.clip(ref, 0, 63)]

    # matches[s, i]: sequence s carries the reference codon at position i.
    matches = idx == ref[np.newaxis, :]
    # Small tolerance so that e.g. threshold=0.9 with 10/10 or 9/10 behaves as expected.
    min_matches = threshold * n_seqs - 1e-9

    codon_conserved = matches.sum(axis=0) >= min_matches
    pair_matches = matches[:, :-1] & matches[:, 1:]
    pair_conserved = pair_matches.sum(axis=0) >= min_matches
    pair_ok = ref_ok[:-1] & ref_ok[1:]

    codon_ref = np.bincount(ref[ref_ok], minlength=64)
    codon_cons = np.bincount(ref[ref_ok & codon_conserved], minlength=64)

    pair_index = ref[:-1].astype(np.int32) * 64 + ref[1:].astype(np.int32)
    pair_ref = np.bincount(pair_index[pair_ok], minlength=4096)
    pair_cons = np.bincount(pair_index[pair_ok & pair_conserved], minlength=4096)

    pair_sense_idx = (sense_idx[:, None] * 64 + sense_idx[None, :]).ravel()

    codons = pd.DataFrame(
        {
            "Codon": list(code.sense_codons),
            "ReferenceCount": codon_ref[sense_idx].astype(np.int64),
            "ConservationCount": codon_cons[sense_idx].astype(np.int64),
        }
    )
    pairs = pd.DataFrame(
        {
            "CodonPair": list(code.sense_pairs),
            "ReferenceCount": pair_ref[pair_sense_idx].astype(np.int64),
            "ConservationCount": pair_cons[pair_sense_idx].astype(np.int64),
        }
    )
    return ConservationResult(codons=codons, pairs=pairs, n_reference_codons=int(ref_ok.sum()))


def conservation_for_records(
    records: Sequence[Tuple[str, str]],
    orf: int = 0,
    *,
    reference_id: str | None = None,
    threshold: float = 1.0,
    genetic_code: int = 1,
) -> ConservationResult:
    """Convenience wrapper over :func:`calculate_conservation` for ``(id, seq)`` records."""
    if not records:
        raise ValueError("No sequences provided")
    ids = [r[0] for r in records]
    ref = resolve_reference(ids, reference_id)
    return calculate_conservation(
        [r[1] for r in records],
        orf,
        reference=ref,
        threshold=threshold,
        genetic_code=genetic_code,
    )
