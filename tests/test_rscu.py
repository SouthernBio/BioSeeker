from pathlib import Path
from tempfile import TemporaryDirectory
from collections import Counter
from bioseeker.analysis.rscu_analysis import (
    calculate_rscu_for_sequence, 
    process_msa_file, 
    get_standard_genetic_code,
    count_codons_for_sequence,
    calculate_rscu_from_counts,
    process_directory
)


def test_genetic_code_excludes_stops():
    aa_map = get_standard_genetic_code()
    all_codons = set(str(c) for codons in aa_map.values() for c in codons)
    stop_codons = {"TAA", "TAG", "TGA"}
    assert not any(stop in all_codons for stop in stop_codons), "Stop codons should be excluded"
    

def test_count_codons():
    seq = "TTTTTCTTCTAA" # 2 TTT, 1 TTC, 1 Stop (excluded in calculate but verify count handles it)
    # The count function cleans and splits. It does count stop codons if they are in triplets.
    # Logic: cleaning removes gaps. 
    counts = count_codons_for_sequence(seq)
    assert counts["TTT"] == 1 # "TTT TT C" -> TTT (1), TTC (1)... wait: TTT TT C T T C T A A
    # "T T T T T C T T C T A A"
    # 0 1 2 3 4 5 6 7 8 9 0 1
    # TTT (1)
    # TTC (1)
    # TTC (1)
    # TAA (1)
    # Correction: "TTT TTC TTC TAA" -> TTT, TTC, TTC, TAA.
    seq_aligned = "TTT-TTC-TTC-TAA"
    counts = count_codons_for_sequence(seq_aligned)
    assert counts["TTT"] == 1
    assert counts["TTC"] == 2
    assert counts["TAA"] == 1


def test_calculate_rscu_values():
    # Phe: TTT, TTC
    # Counts: TTT:2, TTC:0
    counts = Counter({"TTT": 2, "TTC": 0})
    rscu = calculate_rscu_from_counts(counts)
    assert rscu["TTT"] == 2.0
    assert rscu["TTC"] == 0.0


def test_process_msa_file_counts():
    with TemporaryDirectory() as tmpdir:
        input_dir = Path(tmpdir)
        fasta_content = """>SpeciesA
TTTTTT
>SpeciesB
TTCTTC
"""
        fasta_file = input_dir / "test.afa"
        with open(fasta_file, "w") as f:
            f.write(fasta_content)
            
        species_counts = process_msa_file(fasta_file)
        assert "SpeciesA" in species_counts
        assert species_counts["SpeciesA"]["TTT"] == 2
        assert species_counts["SpeciesB"]["TTC"] == 2


def test_process_directory_aggregation():
    with TemporaryDirectory() as tmpdir:
        input_dir = Path(tmpdir) / "input"
        output_file = Path(tmpdir) / "output.csv"
        input_dir.mkdir()
        
        # File 1
        fasta1 = """>SpeciesA
TTT
>SpeciesB
TTC
"""
        with open(input_dir / "1.afa", "w") as f:
            f.write(fasta1)
            
        # File 2
        fasta2 = """>SpeciesA
TTT
>SpeciesB
TTC
"""
        with open(input_dir / "2.afa", "w") as f:
            f.write(fasta2)
            
        process_directory(input_dir, output_file)
        
        assert output_file.exists()
        content = output_file.read_text()
        
        # SpeciesA Total TTT: 2 (1 from each). RSCU(TTT) for SpeciesA should be high/max.
        # SpeciesB Total TTC: 2.
        
        # Verify CSV has columns
        assert "Codon,SpeciesA,SpeciesB" in content or "Codon,SpeciesB,SpeciesA" in content
