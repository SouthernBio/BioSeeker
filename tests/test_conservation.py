from bioseeker.analysis.conservation import conservation_for_records


def test_conservation_for_records_identical():
    # 3 identical sequences
    records = [
        ("Seq1", "ATGGCC"),
        ("Seq2", "ATGGCC"),
        ("Seq3", "ATGGCC"),
    ]
    res = conservation_for_records(records, 0, threshold=1.0)
    c_df = res.codons

    # ATG and GCC should be perfectly conserved
    row_atg = c_df[c_df["Codon"] == "ATG"].iloc[0]
    assert row_atg["ReferenceCount"] == 1
    assert row_atg["ConservationCount"] == 1

    row_gcc = c_df[c_df["Codon"] == "GCC"].iloc[0]
    assert row_gcc["ReferenceCount"] == 1
    assert row_gcc["ConservationCount"] == 1


def test_conservation_for_records_variation():
    # 1 ref: ATGGCC
    # 2 others: ATGGCC, ATGTCC
    records = [
        ("Seq1", "ATGGCC"),
        ("Seq2", "ATGGCC"),
        ("Seq3", "ATGTCC"),
    ]

    res = conservation_for_records(records, 0, threshold=0.9)
    c_df = res.codons

    # ATG
    row_atg = c_df[c_df["Codon"] == "ATG"].iloc[0]
    assert row_atg["ConservationCount"] == 1

    # GCC is matched in 2 out of 3 -> 0.66, which is < 0.9
    row_gcc = c_df[c_df["Codon"] == "GCC"].iloc[0]
    assert row_gcc["ConservationCount"] == 0


def test_gap_only_columns_do_not_split_pairs():
    # Species A and B taken from a 3-species MSA in which only C had an insertion.
    # The leftover all-gap codon must not separate CGA-CGA (pairwise alignment
    # of A and B would have them adjacent).
    records = [
        ("A", "ATGCGA---CGAAAA"),
        ("B", "ATGCGA---CGAAAA"),
    ]
    res = conservation_for_records(records, 0)
    row = res.pairs.set_index("CodonPair").loc["CGACGA"]
    assert row["ReferenceCount"] == 1
    assert row["ConservationCount"] == 1
    assert res.n_reference_codons == 4

    # With the inserting species present, the pair is still not conserved.
    res3 = conservation_for_records(records + [("C", "ATGCGAGGGCGAAAA")], 0)
    row3 = res3.pairs.set_index("CodonPair").loc["CGACGA"]
    assert row3["ReferenceCount"] == 0
    assert row3["ConservationCount"] == 0
