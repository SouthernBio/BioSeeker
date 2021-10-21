# Calculadora de tasa de conservación de codones

***+++ TRABAJO EN PROGRESO +++***

Este script facilita el cálculo de la tasa de conservación de codones. 

### 1. ¿Cómo funciona?

El script toma como input un archivo que contiene genes homólogos previamente alineados, en formato FASTA. Luego de un procesamiento mínimo del archivo para extraer los datos, se crea una matriz mediante la librería `numpy`, y así poder iterar sobre porciones de la misma.

### 2. Dependencias

Para ejecutar este script debe instalar el lenguaje Python y el gestor de paquetes `pip`, en Ubuntu puede realizarlo a través de la terminal:
```bash
$ sudo apt-get update
$ sudo apt-get install python3 python3-pip
```
Una vez instalados Python y su gestor de paquetes, puede proceder a installar `numpy` mediante el siguiente comando:
```bash
$ pip3 install numpy
```

### 3. ¿Cómo se usa?
El archivo `script_conservacion.py` debe colocarse en el directorio donde se encuentra el archivo FASTA que desea analizar. El script hace uso de la librería `os`, por lo que debe modificar manualmente el mismo para añadir el directorio de trabajo donde se encuentran las secuencias a analizar. Este directorio puede obtenerlo automáticamente abriendo una terminal desde la carpeta de las secuencias y tipeando `pwd`.  

