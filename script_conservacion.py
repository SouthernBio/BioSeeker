#   1   -   Función de lectura y almacenamiento de archivos .afa
#           En el argumento de la función se indica el directorio donde se encuentran los archivos de alineamiento
def lista_de_archivos(directorio_de_trabajo):
    import os
    directorio = directorio_de_trabajo
    contenido = os.listdir(directorio)

    #   Se crea una lista de archivos de extensión .afa, los cuales luego serán procesados mediante otra función
    secuencias = []
    print("Programa iniciado. Creando lista de archivos AFA en el directorio " + directorio_de_trabajo)
    for fichero in contenido:
        if os.path.isfile(os.path.join(directorio, fichero)) and fichero.endswith('.afa'):
            secuencias.append(fichero)
    print("Se ha creado la lista de archivos de forma exitosa.")

    return secuencias

#   2   -   Función de extracción de secuencias a partir de archivo .afa
#           Esta función toma como argumento un archivo de la lista creada con lista_de_archivos()
#           para abrirlo y hacer un sorting de las líneas presentes en el mismo
def extraer_secuencias(archivo):
    import os
    textlist, arrX, arrY = [], [], []
    auxvar = True

    #   Almacenamiento de cada línea del archivo en "textlist"
    print("Extrayendo alineamientos del archivo " + archivo)
    archivo_abierto = open(os.path.expanduser(archivo), "r")
    for line in archivo_abierto:
        stripped_line = line.strip()
        line_list = stripped_line.split()
        textlist.append(line_list)
    print("Extracción exitosa a partir del archivo " + archivo)

    #   Sorting de las líneas mediante el booleano auxvar
    #   Este booleano hace las veces de switch, y envía las líneas a una lista u otra
    #   La especie irá hacia arrX y la secuencia irá hacia arrY
    print("Creando array de alineamiento...")
    for k in textlist:
        if auxvar == True:
            arrX.append(k)
            auxvar = False    
        else:
            arrY.append(k)
            auxvar = True

    #   Creación de array a partir de las listas anteriores
    arrXY = []
    for i in range(len(arrX)):
        arrXY.append([arrX[i], arrY[i]])
    print("La creación del array de alineamiento ha sido exitosa para el archivo " + archivo)
    
    return arrXY

#   3   -   Función de cálculo de tasas de conservación de codones y bicodones
#           a partir de un array de secuencias
def calculos_de_conservacion(array_de_secuencias, indice):
    import numpy as np, pandas as pd
    codones, bicodones = [], []
    n = 3

    #   Usando el array obtenido en la función anterior, se procede a dividir las secuencias en fragmentos de 3 y 6 bases
    #   Los pares de codones (bicodones) se arman con un salto "n" de 3 nucleótidos
    print("Creando array de codones y bicodones...")
    for _,k in array_de_secuencias:
        for j in k:   
            codon = [j[i:i+n] for i in range(0, len(j), n)]
            bicodon = [j[i:i+2*n] for i in range(0, len(j), n)]
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
    print("Estableciendo secuencia de referencia...")
    cod3_ref = codones[0]
    cod6_ref = bicodones[0]

    #   Creación de las listas accesorias para almacenar los cálculos
    largo3, largo6 = len(cod3_ref), len(cod6_ref)
    historial3, historial6 = np.zeros(largo3, dtype="i"), np.zeros(largo6, dtype="i")
    his3, his6 = np.zeros(61, dtype="i"), np.zeros(3721, dtype="i")
    his3c, his6c = np.zeros(largo3, dtype="i"), np.zeros(largo6, dtype="i")
    conserva3, conserva6 = np.zeros(61, dtype="i"), np.zeros(3721, dtype="i")

    #   Cálculos
    print("Realizando cálculos de conservación para historial3 e historial6...")
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
    print("Realizando cálculos de conservación para his3c e his6c...")

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
    print("Ensamblando matriz de conservación (cod_ref, hisc)...")

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
    
    print("El ensamblaje ha sido exitoso. Creando dataframes...")

    #   Creación de dataframes y cambio de nombre de las columnas
    df_codones, df_bicodones = pd.DataFrame(array_codones), pd.DataFrame(array_bicodones)
    df_codones = df_codones.rename(columns = {0:'codon', 1:'referencia', 2:'conservacion'})
    df_bicodones = df_bicodones.rename(columns = {0:'bicodon', 1:'referencia', 2:'conservacion'})

    #   Exportación de dataframes a .CSV
    dataframe_cod3 = df_codones.to_csv("historial_codones" + "_" + str(indice) + ".csv", index = False), 
    dataframe_cod6 = df_bicodones.to_csv("historial_bicodones" + "_" + str(indice) + ".csv", index = False)

    return dataframe_cod3, dataframe_cod6

#   4   -   Generación de la lista de archivos
lista = lista_de_archivos("/home/usuario/Documentos/GitHub/Conservacion-de-codones-raros")

#   5   -   Creación de los dataframes
lista_de_dataframes_cod3 , lista_de_dataframes_cod6 = [], []
indice_indicado = 0
for archivo in lista:
    array_de_secuencias = extraer_secuencias(archivo)
    calculos_de_conservacion(array_de_secuencias, indice_indicado)
    indice_indicado += 1

#   6   -   Mover dataframes a otra carpeta para procesarlos con R
import os, shutil
directorio_dataframes = "/home/usuario/Documentos/GitHub/Conservacion-de-codones-raros/"
contenido = os.listdir(directorio_dataframes)

for fichero in contenido:
    if os.path.isfile(os.path.join(directorio_dataframes, fichero)) and fichero.startswith("historial") and fichero.endswith('.csv'):
        shutil.move(directorio_dataframes + fichero, "/home/usuario/Documentos/GitHub/Conservacion-de-codones-raros/dataframes/" + fichero)
