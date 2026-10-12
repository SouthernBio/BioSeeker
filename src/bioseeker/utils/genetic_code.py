"""Genetic code utilities.

Codons are always handled as uppercase DNA strings over the alphabet ``ACGT``.
Translation tables are taken from NCBI (via :mod:`Bio.Data.CodonTable`), so any
standard table ID can be used (1 = Standard, 12 = Alternative Yeast Nuclear, ...).
"""

from __future__ import annotations

from functools import lru_cache
from itertools import product
from typing import Dict, List, Tuple

from Bio.Data import CodonTable

NUCLEOTIDES = "TCAG"
VALID_BASES = frozenset("ACGT")

#: All 64 codons, in canonical TCAG order.
ALL_CODONS: Tuple[str, ...] = tuple("".join(p) for p in product(NUCLEOTIDES, repeat=3))
#: All 4096 codon pairs, in canonical order.
ALL_CODON_PAIRS: Tuple[str, ...] = tuple(a + b for a, b in product(ALL_CODONS, repeat=2))
#: Integer index of every codon in :data:`ALL_CODONS`.
CODON_INDEX: Dict[str, int] = {c: i for i, c in enumerate(ALL_CODONS)}


def is_valid_codon(codon: str) -> bool:
    """Return True if ``codon`` is a 3-letter string made only of A, C, G, T."""
    return len(codon) == 3 and set(codon) <= VALID_BASES


class GeneticCode:
    """A genetic code built from an NCBI translation table.

    Attributes:
        table_id: NCBI translation table ID.
        name: Human readable table name.
        forward: Mapping sense codon -> one-letter amino acid.
        stop_codons: Stop codons of this table.
        sense_codons: Non-terminating codons (61 for the Standard code).
        sense_pairs: All pairs of sense codons (3721 for the Standard code).
    """

    def __init__(self, table_id: int = 1):
        try:
            table = CodonTable.unambiguous_dna_by_id[table_id]
        except KeyError as exc:
            raise ValueError(f"Unknown NCBI genetic code table: {table_id}") from exc

        self.table_id = table_id
        self.name = table.names[0] if table.names else str(table_id)
        self.forward: Dict[str, str] = dict(table.forward_table)
        self.stop_codons: Tuple[str, ...] = tuple(sorted(table.stop_codons))
        self.sense_codons: Tuple[str, ...] = tuple(c for c in ALL_CODONS if c in self.forward)
        self.sense_pairs: Tuple[str, ...] = tuple(
            a + b for a, b in product(self.sense_codons, repeat=2)
        )

    def __repr__(self) -> str:
        return f"GeneticCode(table_id={self.table_id}, name={self.name!r})"

    def translate(self, codon: str) -> str:
        """Translate a codon. Stop codons return ``'*'``; invalid codons raise KeyError."""
        if codon in self.forward:
            return self.forward[codon]
        if codon in self.stop_codons:
            return "*"
        raise KeyError(f"Invalid codon: {codon!r}")

    def is_sense(self, codon: str) -> bool:
        """Return True if ``codon`` is a valid, non-stop codon."""
        return codon in self.forward

    def synonymous_codons(self) -> Dict[str, List[str]]:
        """Return a mapping amino acid -> list of synonymous sense codons."""
        families: Dict[str, List[str]] = {}
        for codon in self.sense_codons:
            families.setdefault(self.forward[codon], []).append(codon)
        return families

    def dipeptide(self, pair: str) -> str:
        """Return the dipeptide (e.g. ``'RR'``) encoded by a 6-nt codon pair."""
        return self.translate(pair[:3]) + self.translate(pair[3:])


@lru_cache(maxsize=None)
def get_genetic_code(table_id: int = 1) -> GeneticCode:
    """Return a cached :class:`GeneticCode` for the given NCBI table ID."""
    return GeneticCode(table_id)
