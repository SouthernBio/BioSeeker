#   3   -   Function to calculate codon and codon pair conservation rates
#           from an array of sequences
#           where ORF can take values in {0, 1, 2} depending on the reading frame
def calculations(sequences_array, indicated_index, ORF):
    if ORF not in {0, 1 ,2}:
        raise Exception("Error: 'ORF' can only take values in {0, 1, 2}")

    import numpy as np
    import pandas as pd
    
    codons = [] 
    bicodons = []
    n = 3

    #   Using the obtained array in the previous function, we proceed to divide the sequences in fragments of 3 and 6 bases
    #   Groups of codons and bicodons are assembled with an "n" step of 3 nucleotides
    print(f"Creando array de codons y bicodons en el ORF+{ORF}...")
    for _,k in sequences_array:
        for j in k:   
            codon = [j[i:i+n] for i in range(ORF, len(j), n)]
            bicodon = [j[i:i+2*n] for i in range(ORF, len(j), n)]
            codons.append(codon)
            bicodons.append(bicodon)  

    #   Converting to Numpy array
    codons = np.asarray(codons) 
    bicodons = np.asarray(bicodons) # Numpy: VisibleDeprecationWarning, dtype = object. Applying the suggested solution just breaks the code.

    #   Pattern lists of codons and codon pairs
    #   The lists, obviously, do not contemplate stop codons
    cod1 = ["AAA", "AAT", "AAC", "AAG", "ATA", "ATT", "ATC", "ATG", "ACA", "ACT", "ACC", "ACG", "AGA", "AGT", "AGC", "AGG", "TAT", "TAC", "TTA", "TTT", "TTC", "TTG", "TCA", "TCT", "TCC", "TCG", "TGT", "TGC", "TGG", "CAA", "CAT", "CAC", "CAG", "CTA", "CTT", "CTC", "CTG", "CCA", "CCT", "CCC", "CCG", "CGA", "CGT", "CGC", "CGG", "GAA", "GAT", "GAC", "GAG", "GTA", "GTT", "GTC", "GTG", "GCA", "GCT", "GCC", "GCG", "GGA", "GGT", "GGC", "GGG"]
    cod2 = cod1
    codcod = []

    for i in cod1:  
        for j in cod2:
            sumaComponentes = i + j
            codcod.append(sumaComponentes)

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
    count = 0
    for i in codon_reference:
        slicer = codons[:,count]
        for j in slicer:
            if j == i:
                codon_history[count] += 1 
        count += 1

    count = 0
    for i in bicodon_reference:
        slicer = bicodons[:,count]
        for j in slicer:
            if j == i:
                bicodon_history[count] += 1
        count += 1

    print("Computation successful.")
    print(f"Computing conservation rates for conserved_codon_history and conserved_bicodon_history at ORF+{ORF}...")

    count1 = 0
    count2 = 0
    aux1 = len(codons[:,0])
    aux2 = len(bicodons[:,0])

    for x in codon_history:
        if x/aux1 > 0.9:
            conserved_codon_history[count1] += 1
        count1 += 1

    for x in bicodon_history:
        if x/aux2 > 0.9:
            conserved_bicodon_history[count2] += 1
        count2 += 1
    
    print("Computation successful.")
    print(f"Assembling conservartion matrix(cod_ref, hisc) for ORF+{ORF}...")

    conservation_matrix_reference_codons = np.column_stack((np.asarray(codon_reference), np.asarray(conserved_codon_history)))
    conservation_matrix_reference_bicodons = np.column_stack((np.asarray(bicodon_reference), np.asarray(conserved_bicodon_history)))

    count = 0
    for codon in cod1:
        for i in conservation_matrix_reference_codons:
            if (i[1] != '0') and (i[0] == codon):
                codon_conservation[count] += 1
        count += 1

    count = 0
    for bicodon in codcod:
        for i in conservation_matrix_reference_bicodons:
            if (i[1] != '0') and (i[0] == bicodon):
                bicodon_conservation[count] += 1
        count += 1    

    count = 0
    for cod in cod1:
        for i in codon_reference:
            if cod == i:
                reference_codon_history[count] += 1
        count += 1

    count = 0
    for bicod in codcod:
        for i in bicodon_reference:
            if bicod == i:
                reference_bicodon_history[count] += 1
        count += 1

    codon_array = np.column_stack((cod1, reference_codon_history, codon_conservation))
    bicodon_array = np.column_stack((codcod, reference_bicodon_history, bicodon_conservation))
    
    print(f"Assembly has been successful. Creating dataframes for ORF+{ORF}...")

    #   Creating dataframes and changing column names
    codon_dataframe = pd.DataFrame(codon_array)
    bicodon_dataframe = pd.DataFrame(bicodon_array)

    codon_dataframe = codon_dataframe.rename(columns = {0:'codon', 1:'reference', 2:'conservation'})
    bicodon_dataframe = bicodon_dataframe.rename(columns = {0:'bicodon', 1:'reference', 2:'conservation'})

    #   Exporting dataframes to CSV files
    dataframe_codons = codon_dataframe.to_csv("history_codons" + "_" + str(indicated_index) + "_" + str(ORF) + ".csv", index = False), 
    dataframe_bicodons = bicodon_dataframe.to_csv("history_bicodons" + "_" + str(indicated_index) + "_" + str(ORF) + ".csv", index = False)

    return dataframe_codons, dataframe_bicodons