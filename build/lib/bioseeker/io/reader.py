from Bio import SeqIO
from typing import List, Tuple
import os


def read_fasta(file_path: str) -> List[Tuple[str, str]]:
    """
    Reads sequences from a FASTA file.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        List[Tuple[str, str]]: A list of tuples where each tuple contains (description, sequence).
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    sequences = []
    # Biopython handles various FASTA sub-formats automatically
    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # We return a format compatible with the expectations of the analysis logic
            sequences.append((record.description, str(record.seq)))
    
    return sequences
