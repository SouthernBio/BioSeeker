import listfiles as lf, seqxtract as sq, conservationrate as cr

#   Generación de la lista de archivos
lista = lf.lista_de_archivos("/home/usuario/Documentos/GitHub/Conservacion-de-codones-raros")

#   Creación de los dataframes
lista_de_dataframes_cod3 , lista_de_dataframes_cod6 = [], []
indice_indicado = 0
for archivo in lista:
    array_de_secuencias = sq.extraer_secuencias(archivo)
    cr.calculos_de_conservacion(array_de_secuencias, indice_indicado, 1)
    cr.calculos_de_conservacion(array_de_secuencias, indice_indicado, 2)
    indice_indicado += 1

#   Mover dataframes a otra carpeta para procesarlos con R
import os, shutil
directorio_dataframes = "/home/usuario/Documentos/GitHub/Conservacion-de-codones-raros/"
contenido = os.listdir(directorio_dataframes)

for fichero in contenido:
    if os.path.isfile(os.path.join(directorio_dataframes, fichero)) and fichero.startswith("historial") and fichero.endswith('.csv'):
        shutil.move(directorio_dataframes + fichero, "/home/usuario/Documentos/GitHub/Conservacion-de-codones-raros/dataframes/" + fichero)
