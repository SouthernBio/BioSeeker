from utils.GeneticCode import CODON_TUPLE, CODON_PAIRS_TUPLE
import numpy as np 
import pandas as pd 


class InvalidReadingFrame(Exception):
    """Custom exception for when a function receives a reading frame out of range"""
    def __init__(self, message):
        super().__init__(message)



def calculations(sequences_array: list, indicated_index: int, ORF: int):
    """Function to calculate codon and codon pair conservation rates

    Args:
        sequences_array (list): An array of sequences
        indicated_index (int): auxiliary index to identify files
        ORF (int): Reading frame. Can only take values in {0, 1, 2}

    Raises:
        InvalidReadingFrame: invalid reading frame.

    Returns:
        pandas.DataFrame: codon dataframe
        pandas.DataFrame: codon pair dataframe
    """
    if ORF not in {0, 1 ,2}:
        raise InvalidReadingFrame("Error: 'ORF' can only take values in {0, 1, 2}")

    codons = []
    bicodons = []
    CODON_SIZE = 3

    #   Using the obtained array in the previous function, we proceed to divide the sequences in fragments of 3 and 6 bases
    #   Groups of codons and bicodons are assembled with an "CODON_SIZE" step of 3 nucleotides
    print(f"Creating codon and codon pair arrays at ORF+{ORF}...")
    for _,k in sequences_array:
        for j in k:
            codon = [j[i:i+CODON_SIZE] for i in range(ORF, len(j), CODON_SIZE)]
            bicodon = [j[i:i+2*CODON_SIZE] for i in range(ORF, len(j), CODON_SIZE)]
            codons.append(codon)
            bicodons.append(bicodon)

    #   Converting to Numpy array
    codons = np.asarray(codons)
    bicodons = np.asarray(bicodons) # Numpy: VisibleDeprecationWarning, dtype = object. Applying the suggested solution just breaks the code.

    #   Reference lists
    print(f"Establishing reference sequences at ORF+{ORF}...")
    codon_reference = codons[0]
    bicodon_reference = bicodons[0]

    #   Creating accessory lists to store the calculations
    number_of_reference_codons = len(codon_reference)
    number_of_reference_bicodons = len(bicodon_reference)

    codon_history = np.zeros(number_of_reference_codons, dtype="i")
    bicodon_history = np.zeros(number_of_reference_bicodons, dtype="i")

    reference_codon_history = np.zeros(61, dtype="i")
    reference_bicodon_history = np.zeros(3721, dtype="i")

    conserved_codon_history = np.zeros(number_of_reference_codons, dtype="i")
    conserved_bicodon_history = np.zeros(number_of_reference_bicodons, dtype="i")

    codon_conservation = np.zeros(61, dtype="i")
    bicodon_conservation = np.zeros(3721, dtype="i")

    #   Calculations
    print(f"Computing conservation for codon_history and bicodon_history at ORF+{ORF}...")
    for count, i in enumerate(codon_reference, start = 0):
        slicer = codons[:,count]
        for j in slicer:
            if j == i:
                codon_history[count] += 1

    for count, i in enumerate(bicodon_reference, start = 0):
        slicer = bicodons[:,count]
        for j in slicer:
            if j == i:
                bicodon_history[count] += 1

    print("Computation successful.")
    print(f"Computing conservation rates for conserved_codon_history and conserved_bicodon_history at ORF+{ORF}...")

    aux1 = len(codons[:,0])
    aux2 = len(bicodons[:,0])

    for count, x in enumerate(codon_history, start = 0):
        if x/aux1 > 0.9:
            conserved_codon_history[count] += 1

    for count, x in enumerate(bicodon_history, start = 0):
        if x/aux2 > 0.9:
            conserved_bicodon_history[count] += 1

    print("Computation successful.")
    print(f"Assembling conservation matrix(cod_ref, hisc) for ORF+{ORF}...") # This part of the code is computationally expensive. We have to find a way to optimize it


    conservation_matrix_reference_codons = np.column_stack((np.asarray(codon_reference), np.asarray(conserved_codon_history)))
    conservation_matrix_reference_bicodons = np.column_stack((np.asarray(bicodon_reference), np.asarray(conserved_bicodon_history)))

    for count, codon in enumerate(CODON_TUPLE, start = 0):
        for i in conservation_matrix_reference_codons:
            if (i[1] != '0') and (i[0] == codon):
                codon_conservation[count] += 1

    for count, bicodon in enumerate(CODON_PAIRS_TUPLE, start = 0):
        for i in conservation_matrix_reference_bicodons:
            if (i[1] != '0') and (i[0] == bicodon):
                bicodon_conservation[count] += 1

    for count, cod in enumerate(CODON_TUPLE, start = 0):
        for i in codon_reference:
            if cod == i:
                reference_codon_history[count] += 1

    for count, bicod in enumerate(CODON_PAIRS_TUPLE, start = 0):
        for i in bicodon_reference:
            if bicod == i:
                reference_bicodon_history[count] += 1

    codon_array = np.column_stack((list(CODON_TUPLE), reference_codon_history, codon_conservation))
    bicodon_array = np.column_stack((list(CODON_PAIRS_TUPLE), reference_bicodon_history, bicodon_conservation))

    print(f"The assembly has been successful. Creating dataframes for ORF+{ORF}...")

    #   Creating dataframes and changing column names
    codon_dataframe = pd.DataFrame(codon_array)
    bicodon_dataframe = pd.DataFrame(bicodon_array)

    codon_dataframe = codon_dataframe.rename(columns = {0:'codon', 1:'ReferenceCount', 2:'ConservationCount'})
    bicodon_dataframe = bicodon_dataframe.rename(columns = {0:'codon_pair', 1:'ReferenceCount', 2:'ConservationCount'})

    #   Exporting dataframes to CSV files
    dataframe_codons = codon_dataframe.to_csv("history_codons" + "_" + str(indicated_index) + "_" + str(ORF) + ".csv", index = False)
    dataframe_bicodons = bicodon_dataframe.to_csv("history_bicodons" + "_" + str(indicated_index) + "_" + str(ORF) + ".csv", index = False)

    return dataframe_codons, dataframe_bicodons
