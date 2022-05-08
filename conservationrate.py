#   3   -   Función de cálculo de tasas de conservación de codones y bicodones
#           a partir de un array de secuencias
#           donde ORF puede tomar valores {0, 1, 2} según el marco de lectura a considerar
def calculos_de_conservacion(array_de_secuencias, indice, ORF):
    if ORF not in {0, 1 ,2}:
        raise Exception("Error: 'ORF' solo puede tomar valores en {0, 1, 2}")

    import numpy as np, pandas as pd
    codones, bicodones = [], []
    n = 3

    #   Usando el array obtenido en la función anterior, se procede a dividir las secuencias en fragmentos de 3 y 6 bases
    #   Los pares de codones (bicodones) se arman con un salto "n" de 3 nucleótidos
    print(f"Creando array de codones y bicodones en el ORF+{ORF}...")
    for _,k in array_de_secuencias:
        for j in k:   
            codon = [j[i:i+n] for i in range(ORF, len(j), n)]
            bicodon = [j[i:i+2*n] for i in range(ORF, len(j), n)]
            codones.append(codon)
            bicodones.append(bicodon)  

    #   Conversión a Numpy Array
    codones, bicodones = np.asarray(codones), np.asarray(bicodones)

    #   Listas patrón de codones y bicodones
    #   Esta lista no contempla, obviamente, los codones de stop
    cod1 = ["AAA", "AAT", "AAC", "AAG", "ATA", "ATT", "ATC", "ATG", "ACA", "ACT", "ACC", "ACG", "AGA", "AGT", "AGC", "AGG", "TAT", "TAC", "TTA", "TTT", "TTC", "TTG", "TCA", "TCT", "TCC", "TCG", "TGT", "TGC", "TGG", "CAA", "CAT", "CAC", "CAG", "CTA", "CTT", "CTC", "CTG", "CCA", "CCT", "CCC", "CCG", "CGA", "CGT", "CGC", "CGG", "GAA", "GAT", "GAC", "GAG", "GTA", "GTT", "GTC", "GTG", "GCA", "GCT", "GCC", "GCG", "GGA", "GGT", "GGC", "GGG"]
    cod2 = cod1
    codcod = []

    for i in cod1:  
        for j in cod2:
            sumaComponentes = i + j
            codcod.append(sumaComponentes)

    #   Listas de referencia
    print(f"Estableciendo secuencia de referencia en el ORF+{ORF}...")
    cod3_ref = codones[0]
    cod6_ref = bicodones[0]

    #   Creación de las listas accesorias para almacenar los cálculos
    largo3, largo6 = len(cod3_ref), len(cod6_ref)
    historial3, historial6 = np.zeros(largo3, dtype="i"), np.zeros(largo6, dtype="i")
    his3, his6 = np.zeros(61, dtype="i"), np.zeros(3721, dtype="i")
    his3c, his6c = np.zeros(largo3, dtype="i"), np.zeros(largo6, dtype="i")
    conserva3, conserva6 = np.zeros(61, dtype="i"), np.zeros(3721, dtype="i")

    #   Cálculos
    print(f"Realizando cálculos de conservación para historial3 e historial6 en el ORF+{ORF}...")
    contador = 0
    for i in cod3_ref:
        slicer = codones[:,contador]
        for j in slicer:
            if j == i:
                historial3[contador] += 1 
        contador += 1

    contador = 0
    for i in cod6_ref:
        slicer = bicodones[:,contador]
        for j in slicer:
            if j == i:
                historial6[contador] += 1
        contador += 1

    print("Cálculo exitoso.")
    print(f"Realizando cálculos de conservación para his3c e his6c en el ORF+{ORF}...")

    contador1, contador2 = 0, 0
    aux1, aux2 = len(codones[:,0]), len(bicodones[:,0])

    for x in historial3:
        if x/aux1 > 0.9:
            his3c[contador1] += 1
        contador1 += 1

    for x in historial6:
        if x/aux2 > 0.9:
            his6c[contador2] += 1
        contador2 += 1
    
    print("Cálculo exitoso.")
    print(f"Ensamblando matriz de conservación (cod_ref, hisc) para el ORF+{ORF}...")

    matriz_conservacion_cod3ref = np.column_stack((np.asarray(cod3_ref), np.asarray(his3c)))
    matriz_conservacion_cod6ref = np.column_stack((np.asarray(cod6_ref), np.asarray(his6c)))

    contador = 0
    for codon in cod1:
        for i in matriz_conservacion_cod3ref:
            if (i[1] != '0') and (i[0] == codon):
                conserva3[contador] += 1
        contador += 1

    contador = 0
    for bicodon in codcod:
        for i in matriz_conservacion_cod6ref:
            if (i[1] != '0') and (i[0] == bicodon):
                conserva6[contador] += 1
        contador += 1    

    contador = 0
    for cod in cod1:
        for i in cod3_ref:
            if cod == i:
                his3[contador] += 1
        contador += 1

    contador = 0
    for bicod in codcod:
        for i in cod6_ref:
            if bicod == i:
                his6[contador] += 1
        contador += 1

    array_codones = np.column_stack((cod1, his3, conserva3))
    array_bicodones = np.column_stack((codcod, his6, conserva6))
    
    print(f"El ensamblaje ha sido exitoso. Creando dataframes para el ORF+{ORF}...")

    #   Creación de dataframes y cambio de nombre de las columnas
    df_codones, df_bicodones = pd.DataFrame(array_codones), pd.DataFrame(array_bicodones)
    df_codones = df_codones.rename(columns = {0:'codon', 1:'referencia', 2:'conservacion'})
    df_bicodones = df_bicodones.rename(columns = {0:'bicodon', 1:'referencia', 2:'conservacion'})

    #   Exportación de dataframes a .CSV
    dataframe_cod3 = df_codones.to_csv("historial_codones" + "_" + str(indice) + "_" + str(ORF) + ".csv", index = False), 
    dataframe_cod6 = df_bicodones.to_csv("historial_bicodones" + "_" + str(indice) + "_" + str(ORF) + ".csv", index = False)

    return dataframe_cod3, dataframe_cod6