# Herramientas de análisis de conservación de codones y pares de codones

Este proyecto facilita el cálculo de la tasa de conservación de codones y pares de codones para un género. Forma parte de mi actual trabajo de tesina de ingeniería.
### 1. ¿Cómo funciona?

El archivo `script_conservacion.py` toma como input un fichero (o conjunto de ficheros) que contiene genes homólogos previamente alineados, en formato FASTA. Luego de un procesamiento mínimo del fichero para extraer los datos, se crea una matriz mediante la librería `numpy`, y así poder iterar sobre porciones de la misma. La información obtenida (conteo de codones en la secuencia de referencia, y cantidad de veces que se conservó el codon) se almacena en un archivo CSV, el cual se crea utilizando la librería `pandas`. Por cada archivo de alineamiento múltiple se generan dos dataframes: uno para codones, y otro para pares de codones. En la carpeta `/dataframes` podrán encontrar dos scripts de R: `dataframes.r` y `analisis.r`. El primer script (`dataframes.r`) realiza un ensamblado de los dataframes localizados en esa carpeta, y devuelve tres nuevos dataframes - uno para codones, otro para pares de codones, y un tercer dataframe opcional que comprende únicamente los pares de codones que han aparecido en las secuencias de referencia y se han conservado al menos una vez. El segundo script (`analisis.r`) procesa la información de los dataframes generados para así determinar qué codones y pares de codones poseen una tasa de conservación superior a la estimada.

### 2. Dependencias

Para ejecutar este script debe instalar el lenguaje Python y el gestor de paquetes `pip`, en Ubuntu puede realizarlo a través de la terminal:
```bash
$ sudo apt-get update
$ sudo apt-get install python3 python3-pip
```
Una vez instalados Python y su gestor de paquetes, puede proceder a installar `numpy` y `pandas` mediante el siguiente comando: 
```bash
$ pip3 install numpy pandas
```
Para el ensamblado de los dataframes generados y su posterior análisis se usa R. Para instalar el lenguaje se usan los siguientes comandos:
```bash
$ sudo apt update
$ sudo apt -y install r-base gdebi-core
```
Luego de descargar el software RStudio (que funciona en base a R) [del sitio oficial](https://rstudio.com/products/rstudio/download/#download), se procede a instalarlo. 

```bash
$ ls
rstudio-1.2.5019-amd64.deb # Ejemplo de archivo .deb descargado, el nombre puede variar según la versión
$ sudo gdebi rstudio-1.2.5019-amd64.deb
```

### 3. ¿Cómo se usa?
El archivo `script_conservacion.py` debe colocarse en el directorio donde se encuentra el conjunto de archivos de alineamiento múltiple que se desean analizar. El script hace uso de la librería `os`, por lo que debe modificar manualmente el mismo para añadir el directorio de trabajo donde se encuentran las secuencias a analizar. Este directorio puede obtenerlo automáticamente abriendo una terminal desde la carpeta de las secuencias y tipeando `pwd`. En el mismo directorio donde se encuentran los archivos de alineamiento múltiple, debe crear una carpeta llamada "dataframes", en donde se localizan dos scripts de R que realizan el ensamblado de los dataframes generados y el análisis correspondiente.
