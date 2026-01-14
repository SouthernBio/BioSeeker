# BioSeeker

![BioSeeker](https://github.com/fx-biocoder/BioSeeker/blob/main/utils/bioseeker.png)

BioSeeker is a Python library and CLI tool for the analysis of codon and bicodon conservation rates across linked species.

This project facilitates calculating codon and bicodon conservation rates for a given genus. It is inspired by the paper [Conservation of location of several specific inhibitory codon pairs in the Saccharomyces sensu stricto yeasts reveals translational selection (Ghoneim et al., 2019)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6379720/).

## Features

- **Robust Parsing**: Uses Biopython to handle standard FASTA/AFA formats.
- **Efficient Analysis**: Vectorized calculations using NumPy for performance.
- **Conservation Metrics**: Calculates conservation rates for both single codons and codon pairs (bicodons) across three reading frames (ORF 0, 1, 2).
- **RSCU Analysis**: Calculates Relative Synonymous Codon Usage (RSCU) for each species across the entire dataset.
- **Aggregation**: Automatically aggregates results from multiple alignment files.
- **Dockerized**: specific container for easy deployment and usage.

## Installation

### From Source

Requires Python 3.9+.

```bash
git clone https://github.com/SouthernBio/BioSeeker.git
cd BioSeeker
pip install .
```

### Development

```bash
python3 -m venv venv
source venv/bin/activate
pip install -e .
pip install pytest
```

## Usage

### Command Line

BioSeeker provides a simple CLI. Place your aligned FASTA files (extension `.afa`, `.fasta`, or `.fa`) in a directory (e.g., `FASTA_files`).

```bash
bioseeker --input FASTA_files --output results
```

**Arguments:**
- `--input`, `-i`: Directory containing input alignment files (default: `FASTA_files`).
- `--output`, `-o`: Directory to save CSV results (default: `results`).

### RSCU Analysis

To calculate the global RSCU values:

```bash
python src/bioseeker/analysis/rscu_analysis.py FASTA_files results/global_rscu.csv
```

### Docker

You can run BioSeeker without installing Python dependencies using Docker.

**Build the image:**
```bash
docker build -t bioseeker .
```

**Run the container:**
Assuming your data is in the current directory under `FASTA_files`:

```bash
docker run -v $(pwd)/FASTA_files:/data/input -v $(pwd)/results:/data/output bioseeker --input /data/input --output /data/output
```

## Output

The tool generates CSV files in the output directory:
- `codon_conservation_ORF0.csv`, `ORF1.csv`, `ORF2.csv`
- `bicodon_conservation_ORF0.csv`, `ORF1.csv`, `ORF2.csv`

Each file contains:
- `codon` / `codon_pair`: The sequence.
- `ReferenceCount`: Occurrences in the reference sequence.
- `ConservationCount`: Number of times the codon was conserved (>90% match) across alignments.
- `ConservationRate`: `ConservationCount / ReferenceCount`.

## License

GNU General Public License v3.0. See [LICENSE](LICENSE.md).

## 💙 Support this project

Your contribution would help SouthernBio in improving the quality of this project and adding additional features. If you find this project useful and/or interesting, please consider offering your support on Github Sponsors, Ko-Fi or PayPal

[![Github-sponsors](https://img.shields.io/badge/sponsor-30363D?style=for-the-badge&logo=GitHub-Sponsors&logoColor=#EA4AAA)](https://github.com/sponsors/fx-biocoder) [![Ko-Fi](https://img.shields.io/badge/Ko--fi-F16061?style=for-the-badge&logo=ko-fi&logoColor=white)](https://ko-fi.com/biocoder) [![PayPal](https://img.shields.io/badge/PayPal-00457C?style=for-the-badge&logo=paypal&logoColor=white)](https://paypal.me/facumartinez680)
