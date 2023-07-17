# BioSeeker: Algorithm for the analysis of codon/bicodon conservation rates across linked species

This project facilitates calculating codon and bicodon conservation rates for a given genus.

### Useful Links:
* [Project backlog (in progress)](https://github.com/fx-biocoder/codon-conservation-rate/projects)
* [Wiki (in progress)](https://github.com/fx-biocoder/codon-conservation-rate/wiki)
* [Paper: Conservation of location of several specific inhibitory codon pairs in the Saccharomyces sensu stricto yeasts reveals translational selection (Ghoneim et al., 2019)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6379720/)

### 1. About the project

This project started in 2021 as part of my thesis to obtain an engineering degree, but due to personal reasons I've decided to change my career path and the project was left abandoned. I think it would be beneficial for others to make use of my work, so I'm now releasing it as an open-source bioinformatics project licensed under GPL v3.0, everyone is welcome to participate in it. The inspiration for this project came from [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6379720/), which I've tried to (partially) replicate using Drosophila's alignments from [FlyDIVaS](https://flydivas.info). Feel free to make as many additions as you'd like.

### 2. How does it work?

In this repo you will find a Python script called `main.py`. It takes a file (or a group of files) as input, which contains homologous genes previously aligned, in FASTA format. After parsing the file(s) for data extraction, it creates a matrix using `Numpy`, in order to iterate across matrix slices. The obtained information (codon count from reference sequence, and number of times that said codon was conserved across species) is stored in a CSV file, which is created using `Pandas`. For each MSA file two types of dataframes will be generated - one for codons, and another for codon pairs. The algorithm calculates codon/bicodon conservation rates across all 3 reading frames. So, there will be a total of 6 CSV files that will be generated (2 for each reading frame). 

### 3. Running the algorithm

Copy and paste the MSA FASTA files on the directory where `main.py` is located. Then, you can just type `./main.py` in the terminal.

### 4. Dependencies

To execute this script you must install Python and its package manager, `pip`. You can do it on Ubuntu through the terminal:
```bash
$ sudo apt-get update
$ sudo apt-get install python3 python3-pip
```
Once you have installed Python and its package manager, you can proceed to install `pipenv` and run the virtual environment:
```bash
$ pip3 install pipenv
$ pipenv install
$ pipenv shell
```

### 5. Additional details
After parsing the files and calculating conservation rates, it will also generate a file called `unreadable.txt` which stores the names of MSA files that could not be parsed. Then, it will assemble all individual dataframes into 6 different dataframes that contain all the information across linked species. After the analysis is over, unnecessary dataframes will be deleted using `delete_files.bat` or `delete_files.sh` depending on your current OS (NOTE: you must edit the shell scripts to include your current working directory).

## 💙 Contributing

Your contribution would help me a lot in finishing this project. If you find this project useful and/or interesting, please consider supporting me on Github Sponsors, Ko-Fi or PayPal

[![Github-sponsors](https://img.shields.io/badge/sponsor-30363D?style=for-the-badge&logo=GitHub-Sponsors&logoColor=#EA4AAA)](https://github.com/sponsors/fx-biocoder) [![Ko-Fi](https://img.shields.io/badge/Ko--fi-F16061?style=for-the-badge&logo=ko-fi&logoColor=white)](https://ko-fi.com/biocoder) [![PayPal](https://img.shields.io/badge/PayPal-00457C?style=for-the-badge&logo=paypal&logoColor=white)](https://paypal.me/facumartinez680)
