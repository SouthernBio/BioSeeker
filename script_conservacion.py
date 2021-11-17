import numpy as np

#   Lectura de archivos .afa y almacenamiento en una lista
import os
directorio = '/home/usuario/Documentos/GitHub/Conservacion-de-codones-raros'
contenido = os.listdir(directorio)

secuencias = []
for fichero in contenido:
    if os.path.isfile(os.path.join(directorio, fichero)) and fichero.endswith('.afa'):
        secuencias.append(fichero)

#   Array(arrX, arrY) donde arrX=especies, arrY=secuencias
#   El booleano "auxvar" hace las veces de switch para hacer un sorting de las líneas

textlist, arrX, arrY = [], [], []
auxvar = True

for i in range(len(secuencias)):
    archivo = open(os.path.expanduser(secuencias[i]), "r")
    for line in archivo:
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


#   Separación en codones y dicodones
#   la lista de dicodones tiene en cuenta la superposición de codones
#   n = 3 corresponde al tamaño del salto a la hora de definir la nueva secuencia a añadir a las listas
codones, bicodones = [], []
n = 3

for _,k in arrXY:
    for j in k:   
        codon = [j[i:i+n] for i in range(0, len(j), n)]
        bicodon = [j[i:i+2*n] for i in range(0, len(j), n)]
        codones.append(codon)
        bicodones.append(bicodon)   #   Problema a resolver: queda un codon al final de la lista...

codones, bicodones = np.asarray(codones), np.asarray(bicodones)

#   Listas de codones y pares de codones
#   Generación de lista de pares de codones mediante un bucle

cod1 = ["AAA", "AAT", "AAC", "AAG", "ATA", "ATT", "ATC", "ATG", "ACA", "ACT", "ACC", "ACG", "AGA", "AGT", "AGC", "AGG", "TAT", "TAC", "TTA", "TTT", "TTC", "TTG", "TCA", "TCT", "TCC", "TCG", "TGT", "TGC", "TGG", "CAA", "CAT", "CAC", "CAG", "CTA", "CTT", "CTC", "CTG", "CCA", "CCT", "CCC", "CCG", "CGA", "CGT", "CGC", "CGG", "GAA", "GAT", "GAC", "GAG", "GTA", "GTT", "GTC", "GTG", "GCA", "GCT", "GCC", "GCG", "GGA", "GGT", "GGC", "GGG"]
cod2 = cod1
codcod = []

for i in cod1:  
    for j in cod2:
        sumaComponentes = i + j
        codcod.append(sumaComponentes)

#   Cálculo de las frecuencias de conservación:

cod3_ref = codones[0]
cod6_ref = bicodones[0]

largo3, largo6 = len(cod3_ref), len(cod6_ref)
historial3, historial6 = np.zeros(largo3, dtype="i"), np.zeros(largo6, dtype="i")
his3, his6 = np.zeros(61, dtype="i"), np.zeros(3721, dtype="i")
his3c, his6c = np.zeros(largo3, dtype="i"), np.zeros(largo6, dtype="i")

#   Codones:
#   Acá se cuenta cuántas veces aparece el codón de referencia en dicha posición
#   historial3 almacena la cantidad de veces que aparece ese codón
contador = 0
for i in cod3_ref:
    slicer = codones[:,contador]
    for j in slicer:
        if j == i:
            historial3[contador] += 1 
    contador += 1

#   Bicodones
#   Comportamiento similar al bloque anterior, pero para bicodones
contador = 0
for i in cod6_ref:
    slicer = bicodones[:,contador]
    for j in slicer:
        if j == i:
            historial6[contador] += 1
    contador += 1

#   Frec. de conservación:
#   Este bloque cuenta cuántas veces el codón (o bicodón) fue conservado en dicha posición
#   Si la frecuencia de conservación es mayor al 90%, entonces suma 1 en la posición de his3c
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

print("Acá va historial3:")
print(historial3)

print("Acá va historial6: ")
print(historial6)

print("Acá va his3c:")
print(his3c)

#   Cálculo de his3/his6
#   Este bloque cuenta cuántas veces un cierto codón (o bicodón) apareció en la secuencia de referencia
contador = 0
for cod in cod1:
    for i in cod3_ref:
        if cod == i:
            his3[contador] += 1
    contador += 1

print("Acá va his3:")
print(his3)

"""
contador = 0
for bicod in codcod:
    for i in cod6_ref:
        if bicod == i:
            his6[contador] += 1 #IndexError: index 399 is out of bounds for axis 0 with size 396
    contador += 1
"""