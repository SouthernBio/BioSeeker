# Algorithm for the analysis of codon/bicodon conservation rates across linked species

This project facilitates calculing codon and bicodon conservation rates for a given genus. 

### 1. How does it work?

In this repo you will find a Python script called `main.py`. It takes a file (or a group of files) as input, which contains homologous genes previously aligned, in FASTA format. After parsing the file(s) for data extraction, it creates a matrix using `Numpy`, in order to iterate across matrix slices. The obtained information (codon count from reference sequence, and number of times that said codon was conserved across species) is stored in a CSV file, which is created using 
`Pandas`. For each MSA file two types of dataframes will be generated - one for codons, and another for codon pairs. The algorithm calculates codon/bicodon conservation rates across all 3 ORFs. So, there will be a total of 6 CSV files that will be generated. 

### 2. Dependencies

To execute this script you must install Python and its package manager, `pip`. You can do it on Ubuntu through the terminal:
```bash
$ sudo apt-get update
$ sudo apt-get install python3 python3-pip
```
Once you have installed Python and its package manager, you can proceed to install `numpy` and `pandas`:
```bash
$ pip3 install numpy pandas
```

### 3. Additional details
`main.py` must be edited to include the directory in which the MSA files are located. After parsing the files and calculating conservation rates, it will also generate a file called `unreadable.txt` which stores the names of MSA files that could not be parsed. Then, it will assemble all individual dataframes into 6 different dataframes that contain all the information across linked species. After the analysis is over, unnecessary dataframes will be deleted using `delete_files.bat`.

### 4. Contributing
Contributions are always welcome!
