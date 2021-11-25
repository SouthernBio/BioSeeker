#   1   -   Función de lectura y almacenamiento de archivos .afa
def lista_de_archivos(directorio_de_trabajo):
    import os
    directorio = directorio_de_trabajo
    contenido = os.listdir(directorio)

    secuencias = []
    for fichero in contenido:
        if os.path.isfile(os.path.join(directorio, fichero)) and fichero.endswith('.afa'):
            secuencias.append(fichero)

    return secuencias

#   2   -   Función de extracción de secuencias a partir de archivo .afa
def extraer_secuencias(archivo):
    import os
    textlist, arrX, arrY = [], [], []
    auxvar = True

    archivo_abierto = open(os.path.expanduser(archivo), "r")
    for line in archivo_abierto:
        stripped_line = line.strip()
        line_list = stripped_line.split()
        textlist.append(line_list)
            
    for k in textlist:
        if auxvar == True:
            arrX.append(k)
            auxvar = False    
        else:
            arrY.append(k)
            auxvar = True

    arrXY = []
    for i in range(len(arrX)):
        arrXY.append([arrX[i], arrY[i]])
    
    return arrXY

#   3   -   Función de cálculo de tasas de conservación de codones y bicodones
#           a partir de un array de secuencias
def calculos_de_conservacion(array_de_secuencias, indice):
    import numpy as np, pandas as pd
    codones, bicodones = [], []
    n = 3

    for _,k in array_de_secuencias:
        for j in k:   
            codon = [j[i:i+n] for i in range(0, len(j), n)]
            bicodon = [j[i:i+2*n] for i in range(0, len(j), n)]
            codones.append(codon)
            bicodones.append(bicodon)  

    codones, bicodones = np.asarray(codones), np.asarray(bicodones)

    cod1 = ["AAA", "AAT", "AAC", "AAG", "ATA", "ATT", "ATC", "ATG", "ACA", "ACT", "ACC", "ACG", "AGA", "AGT", "AGC", "AGG", "TAT", "TAC", "TTA", "TTT", "TTC", "TTG", "TCA", "TCT", "TCC", "TCG", "TGT", "TGC", "TGG", "CAA", "CAT", "CAC", "CAG", "CTA", "CTT", "CTC", "CTG", "CCA", "CCT", "CCC", "CCG", "CGA", "CGT", "CGC", "CGG", "GAA", "GAT", "GAC", "GAG", "GTA", "GTT", "GTC", "GTG", "GCA", "GCT", "GCC", "GCG", "GGA", "GGT", "GGC", "GGG"]
    cod2 = cod1
    codcod = []

    for i in cod1:  
        for j in cod2:
            sumaComponentes = i + j
            codcod.append(sumaComponentes)

    cod3_ref = codones[0]
    cod6_ref = bicodones[0]

    largo3, largo6 = len(cod3_ref), len(cod6_ref)
    historial3, historial6 = np.zeros(largo3, dtype="i"), np.zeros(largo6, dtype="i")
    his3, his6 = np.zeros(61, dtype="i"), np.zeros(3721, dtype="i")
    his3c, his6c = np.zeros(largo3, dtype="i"), np.zeros(largo6, dtype="i")
    conserva3, conserva6 = np.zeros(61, dtype="i"), np.zeros(3721, dtype="i")

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
    
    df_codones, df_bicodones = pd.DataFrame(array_codones), pd.DataFrame(array_bicodones)

    return df_codones.to_csv("historial_codones" + "_" + indice + ".csv"), df_bicodones.to_csv("historial_bicodones" + "_" + indice + ".csv")

#   4   -   Generación de la lista de archivos
lista = lista_de_archivos("C:\\Users\\Usuario\\Desktop\\Tesis\\droso12spp\\archivos")

#   5   -   Creación de los dataframes
indice = 0
for archivo in lista:
    array_de_secuencias = extraer_secuencias(archivo)
    calculos_de_conservacion(array_de_secuencias)
    indice += 1

