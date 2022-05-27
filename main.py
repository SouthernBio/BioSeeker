import listfiles as lf
import seqxtract as sq
import conservationrate as cr
from assembler import assembler
import subprocess

#   Generating list of files
lista = lf.lista_de_archivos(r"C:\Users\marti\Documents\GitHub\Conservacion-de-codones-raros")

#   Creating dataframes
lista_de_dataframes_cod3 , lista_de_dataframes_cod6, archivos_no_leidos = [], [], []
indice_indicado = 0
for archivo in lista:
    try:
        array_de_secuencias = sq.extraer_secuencias(archivo)
        cr.calculos_de_conservacion(array_de_secuencias, indice_indicado, 1)
        cr.calculos_de_conservacion(array_de_secuencias, indice_indicado, 2)
        indice_indicado += 1
    except:
        print(f"Error al analizar el archivo {archivo}")
        archivos_no_leidos.append(archivo)

print("Análisis terminado. Se han omitido los siguientes archivos en el análisis:\n", archivos_no_leidos)

#   Unreadable MSA files
textfile = open("unreadable.txt", "w")
for archivo in archivos_no_leidos:
    textfile.write(archivo + "\n")
textfile.close()

#   Dataframe assembly
directorio = "C:\\Users\\marti\\Documents\\GitHub\\Conservacion-de-codones-raros\\"
#assembler(directorio, 0)   # ORF+0
assembler(directorio, 1)    # ORF+1
assembler(directorio, 2)    # ORF+2

# Delete unnecessary CSV files
subprocess.call([r'C:\Users\marti\Documents\GitHub\Conservacion-de-codones-raros\delete_files.bat'])