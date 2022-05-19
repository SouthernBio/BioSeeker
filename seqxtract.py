#   2   -   Función de extracción de secuencias a partir de archivo .afa
#           Esta función toma como argumento un archivo de la lista creada con lista_de_archivos()
#           para abrirlo y hacer un sorting de las líneas presentes en el mismo
def extraer_secuencias(archivo):
    # Chequeo de extensión
    if archivo.endswith('.afa') == False:
        raise Exception('Error: extensión de archivo inválida, debe tener extensión ".afa". ¿Has modificado el script "listfiles.py"?')
    
    # Importaciones y elementos necesarios
    import os
    textlist, arrX, arrY = [], [], []
    auxvar = True

    #   Almacenamiento de cada línea del archivo en "textlist"
    print(f"Extrayendo alineamientos del archivo {archivo}")
    archivo_abierto = open(os.path.expanduser(archivo), "r")
    for line in archivo_abierto:
        stripped_line = line.strip()
        line_list = stripped_line.split()
        textlist.append(line_list)
    print(f"Extracción exitosa a partir del archivo {archivo}")

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
    print(f"La creación del array de alineamiento ha sido exitosa para el archivo {archivo}")
    
    return arrXY