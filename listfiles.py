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