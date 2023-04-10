# Assembler

def assembler(directorio, ORF):
    if ORF not in {0,1,2}:
        raise Exception("Error. ORF must take values in {0,1,2}, but " + str(ORF) + " was given.")
        
    import pandas as pd
    import os
    import numpy as np

    cod1 = ["AAA", "AAT", "AAC", "AAG", "ATA", "ATT", "ATC", "ATG", "ACA", "ACT", "ACC", "ACG", "AGA", "AGT", "AGC", "AGG", "TAT", "TAC", "TTA", "TTT", "TTC", "TTG", "TCA", "TCT", "TCC", "TCG", "TGT", "TGC", "TGG", "CAA", "CAT", "CAC", "CAG", "CTA", "CTT", "CTC", "CTG", "CCA", "CCT", "CCC", "CCG", "CGA", "CGT", "CGC", "CGG", "GAA", "GAT", "GAC", "GAG", "GTA", "GTT", "GTC", "GTG", "GCA", "GCT", "GCC", "GCG", "GGA", "GGT", "GGC", "GGG"]
    cod2 = cod1
    codcod = []
    for i in cod1:  
        for j in cod2:
            suma_componentes = i + j
            codcod.append(suma_componentes)

    ceros1 = np.zeros(61)
    ceros2 = np.zeros(3721)

    df_codons = pd.DataFrame({'referencia' : ceros1, 'conservacion' : ceros1}, index = cod1)
    df_bicodons = pd.DataFrame({'referencia' :ceros2 , 'conservacion' : ceros2}, index = codcod)
    

    contents = os.listdir(directorio)

    for file in contents:
        if os.path.isfile(os.path.join(directorio, file)) and file.startswith("historial_codones") and file.endswith(f'{ORF}.csv'):
            dataframe = pd.read_csv(file, header = 0, index_col = 0)
            df_codons = df_codons + dataframe
        elif os.path.isfile(os.path.join(directorio, file)) and file.startswith("historial_bicodones") and file.endswith(f'{ORF}.csv'):
            dataframe = pd.read_csv(file, header = 0, index_col = 0)
            df_bicodons = df_bicodons + dataframe

    datos_codones = df_codons.to_csv(f'codon_data_ORF+{ORF}.csv')
    datos_bicodones = df_bicodons.to_csv(f'bicodon_data_ORF+{ORF}.csv')

    return datos_codones, datos_bicodones