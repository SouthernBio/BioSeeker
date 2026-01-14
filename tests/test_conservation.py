from bioseeker.analysis.conservation import calculate_conservation, split_into_codons


def test_split_into_codons():
    seq = "ATGGCC"
    assert split_into_codons(seq, 0) == ["ATG", "GCC"]
    assert split_into_codons(seq, 1) == ["TGG"] # 1:4 TGG, next 4:7 (out)
    assert split_into_codons(seq, 2) == ["GGC"] 

def test_calculate_conservation_identical():
    # 3 identical sequences
    seqs = ["ATGGCC", "ATGGCC", "ATGGCC"]
    c_df, b_df = calculate_conservation(seqs, 0)
    
    # ATG and GCC should be perfectly conserved
    # ATG
    row_atg = c_df[c_df['codon'] == 'ATG'].iloc[0]
    assert row_atg['ReferenceCount'] == 1
    assert row_atg['ConservationCount'] == 1 # 3/3 > 0.9
    
    # GCC
    row_gcc = c_df[c_df['codon'] == 'GCC'].iloc[0]
    assert row_gcc['ReferenceCount'] == 1
    assert row_gcc['ConservationCount'] == 1

def test_calculate_conservation_variation():
    # 1 ref: ATGGCC
    # 2 others: ATGGCC, ATGTCC
    seqs = ["ATGGCC", "ATGGCC", "ATGTCC"]
    # Pos 1: ATG (all match) -> Conserved
    # Pos 2: GCC vs GCC vs TCC. Match 2/3 = 0.66. Not conserved (>0.9)
    
    c_df, b_df = calculate_conservation(seqs, 0)
    
    # ATG
    row_atg = c_df[c_df['codon'] == 'ATG'].iloc[0]
    assert row_atg['ConservationCount'] == 1
    
    # GCC
    row_gcc = c_df[c_df['codon'] == 'GCC'].iloc[0]
    assert row_gcc['ConservationCount'] == 0 # Not conserved
