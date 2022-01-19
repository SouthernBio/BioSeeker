# Calculadora de tasa de conservación de codones y pares de codones

Este proyecto facilita el cálculo de la tasa de conservación de codones y pares de codones para un género. Forma parte de mi actual trabajo de tesina de ingeniería.
### 1. ¿Cómo funciona?

El script toma como input un archivo (o conjunto de archivos) que contiene genes homólogos previamente alineados, en formato FASTA. Luego de un procesamiento mínimo del archivo para extraer los datos, se crea una matriz mediante la librería `numpy`, y así poder iterar sobre porciones de la misma. La información obtenida (conteo de codones en la secuencia de referencia, y cantidad de veces que se conservó el codon) se almacena en un archivo CSV, el cual se crea utilizando la librería `pandas`

### 2. Dependencias

Para ejecutar este script debe instalar el lenguaje Python y el gestor de paquetes `pip`, en Ubuntu puede realizarlo a través de la terminal:
```bash
$ sudo apt-get update
$ sudo apt-get install python3 python3-pip
```
Una vez instalados Python y su gestor de paquetes, puede proceder a installar `numpy` y `pandas` mediante los siguientes comandos: 
```bash
$ pip3 install numpy
$ pip3 install pandas
```
Para el ensamblado de los dataframes (y su posterior análisis, el cual no está contemplado en este repositorio por el momento) generados se usa R. Para instalar el lenguaje se usan los siguientes comandos:
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
El archivo `script_conservacion.py` debe colocarse en el directorio donde se encuentra el conjunto de archivos de alineamiento múltiple que se desean analizar. El script hace uso de la librería `os`, por lo que debe modificar manualmente el mismo para añadir el directorio de trabajo donde se encuentran las secuencias a analizar. Este directorio puede obtenerlo automáticamente abriendo una terminal desde la carpeta de las secuencias y tipeando `pwd`. En el mismo directorio donde se encuentran los archivos de alineamiento múltiple, debe crear una carpeta llamada "dataframes", en donde se localiza un script de R que realiza el ensamblado de los dataframes generados.
