# BioSeeker: Python library for the analysis of codon/bicodon conservation rates across linked species

![BioSeeker (c) 2023 All rights reserved](https://github.com/fx-biocoder/BioSeeker/blob/main/utils/bioseeker.png)
This project facilitates calculating codon and bicodon conservation rates for a given genus.

### Useful Links:
* [Project backlog (in progress)](https://github.com/SouthernBio/BioSeeker/projects)
* [Wiki (in progress)](https://github.com/SouthernBio/BioSeeker/wiki)
* [Paper: Conservation of location of several specific inhibitory codon pairs in the Saccharomyces sensu stricto yeasts reveals translational selection (Ghoneim et al., 2019)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6379720/)

### 1. About the project

Open-source bioinformatics project licensed under GPL v3.0 The inspiration for this project came from [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6379720/), which I've tried to (partially) replicate using Drosophila's alignments from [FlyDIVaS](https://flydivas.info). Feel free to make as many additions as you'd like.

### 2. How does it work?

In this repo you will find a Python script called `bioseeker.py`. It takes a file (or a group of files) as input, which contains homologous genes previously aligned, in FASTA format. After parsing the file(s) for data extraction, it creates a matrix using `NumPy`, in order to iterate across matrix slices. The obtained information (codon count from reference sequence, and number of times that said codon was conserved across species) is stored in a CSV file, which is created using `Pandas`. For each MSA file two types of dataframes will be generated - one for codons, and another for codon pairs. The algorithm calculates codon/bicodon conservation rates across all 3 reading frames. So, there will be a total of 6 CSV files that will be generated (2 for each reading frame). 

### 3. Installing and running the program

Start by cloning the repository:

```bash
$ git clone https://github.com/fx-biocoder/BioSeeker
```

Copy and paste the MSA FASTA files on the directory where `bioseeker.py` is located. Then, you can just type `./bioseeker` in the terminal.

If you want to test how the program works before using it on your data, you will find alignment files on `FASTA_files/`. 

### 4. Dependencies

To execute this script you must install Python, Git (if you want to clone the repo with Git) and its package manager, `pip`. You can do it on Ubuntu through the terminal:
```bash
$ sudo apt-get update
$ sudo apt-get install python3 python3-pip git-all
```
Once you have installed Python and its package manager, you can proceed to install `pipenv`:
```bash
$ pip3 install pipenv
```
Then, you can activate the virtual environment and run BioSeeker:
```bash
$ pipenv shell
$ bioseeker
```

### 5. Additional details
After parsing the files and calculating conservation rates, it will also generate a file called `unreadable.txt` which stores the names of MSA files that could not be parsed. Then, it will assemble all individual dataframes into 6 different dataframes that contain all the information across linked species.

## 💙 Support this project

Your contribution would help SouthernBio in improving the quality of this project and adding additional features. If you find this project useful and/or interesting, please consider offering your support on Github Sponsors, Ko-Fi or PayPal

[![Github-sponsors](https://img.shields.io/badge/sponsor-30363D?style=for-the-badge&logo=GitHub-Sponsors&logoColor=#EA4AAA)](https://github.com/sponsors/fx-biocoder) [![Ko-Fi](https://img.shields.io/badge/Ko--fi-F16061?style=for-the-badge&logo=ko-fi&logoColor=white)](https://ko-fi.com/biocoder) [![PayPal](https://img.shields.io/badge/PayPal-00457C?style=for-the-badge&logo=paypal&logoColor=white)](https://paypal.me/facumartinez680)
