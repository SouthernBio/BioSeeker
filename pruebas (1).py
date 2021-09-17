#   SCRIPT PARA CÁLCULO DE FRECUENCIAS DE CONSERVACIÓN DE CODONES Y PARES DE CODONES
#   Autor: Facundo Martínez
#   Año 2021

#   Lectura de archivos .afa y almacenamiento en una lista
import os
directorio = '/home/usuario/Escritorio/Tesis/droso12spp/'
contenido = os.listdir(directorio)

secuencias = []
for fichero in contenido:
    if os.path.isfile(os.path.join(directorio, fichero)) and fichero.endswith('.afa'):
        secuencias.append(fichero)

#   Array(arrX, arrY) donde arrX=especies, arrY=secuencias
#   El booleano "auxvar" hace las veces de switch para hacer un sorting de las líneas

textlist = []
arrX = []
arrY = []
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

array = []
for i in range(len(arrX)):
    array.append([arrX[i], arrY[i]])

"""
#   Separación en codones y dicodones
#   la lista de dicodones tiene en cuenta la superposición de codones
#   n = 3 corresponde al tamaño del salto a la hora de definir la nueva secuencia a añadir a las listas
codones = []
bicodones = []
n = 3

for _,k in array:
    for j in k:   
        codon = [j[i:i+n] for i in range(0, len(j), n)]
        dicodon = [j[i:i+2*n] for i in range(0, len(j), n)]
        codones.append(codon)
        bicodones.append(dicodon)   #   Problema a resolver: queda un codon al final de la lista...
"""
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

"""
his3 = []
his6 = []

for i in cod1:
    contador = 0
    for j in codones:
        for k in j:
            if k == i:
                contador += 1
    his3.append(contador)

for i in codcod:
    contador = 0
    for j in bicodones:
        for k in j:
            if k == i:
                contador += 1
    his6.append(contador)

print(his3)
"""
codones = []
contador = 0

for i in arrX:
    if arrX[contador] == '>Dmel':
        print(arrX[contador])
        contador += 1
        
#   Para iterar respecto al primer item de codones/bicodones puedo usar un bucle "while"
#   for i in codones: for j in i: while i "menor que cierto valor"
#   me servirá un bucle while?


"""
########### CÓDIGO ADICIONAL ###########

#   Visualizar longitud de codones y bicodones
element = 0
contador = 0
for i in codones/bicodones:
    largo = len(codones/bicodones[element])
    element += 1
    contador += largo
print(contador)

#   Ver tamaño de objeto en memoria:
from sys import getsizeof
print(getsizeof(objeto))"""