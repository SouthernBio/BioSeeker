"""Reading of (aligned) FASTA files."""

from __future__ import annotations

import gzip
import os
from pathlib import Path
from typing import IO, List, Tuple, Union

from Bio import SeqIO

PathLike = Union[str, "os.PathLike[str]"]

#: File extensions recognised as alignments (optionally followed by ``.gz``).
ALIGNMENT_EXTENSIONS = (".afa", ".fasta", ".fa", ".fas", ".fna")


def _is_alignment_file(path: Path) -> bool:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return name.endswith(ALIGNMENT_EXTENSIONS)


def _open_text(path: Path) -> IO[str]:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def discover_alignments(path: PathLike) -> List[Path]:
    """Return the alignment files found at ``path``, sorted by name.

    ``path`` can be a single file (returned as-is) or a directory, in which case
    every file with a known alignment extension (see :data:`ALIGNMENT_EXTENSIONS`,
    optionally gzip-compressed) is returned. The order is deterministic.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
    """
    p = Path(path)
    if p.is_file():
        return [p]
    if p.is_dir():
        return sorted(f for f in p.iterdir() if f.is_file() and _is_alignment_file(f))
    raise FileNotFoundError(f"Input path not found: {p}")


def read_fasta(file_path: PathLike) -> List[Tuple[str, str]]:
    """Read sequences from a (possibly gzip-compressed) FASTA file.

    Args:
        file_path: Path to the FASTA file.

    Returns:
        A list of ``(record_id, sequence)`` tuples. Sequences are uppercased.
        ``record_id`` is the first word of the FASTA header (e.g. ``"Dmel"``).
    """
    p = Path(file_path)
    if not p.exists():
        raise FileNotFoundError(f"File not found: {p}")

    with _open_text(p) as handle:
        return [(record.id, str(record.seq).upper()) for record in SeqIO.parse(handle, "fasta")]
