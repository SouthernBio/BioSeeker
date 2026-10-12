# BioSeeker

[![PyPI Version](https://img.shields.io/pypi/v/bioseeker.svg)](https://pypi.org/project/bioseeker/)
[![Ubuntu PPA](https://img.shields.io/badge/Ubuntu%20PPA-ppa%3Afx--biocoder%2Fbioseeker-E95420.svg?logo=ubuntu&logoColor=white)](https://launchpad.net/~fx-biocoder/+archive/ubuntu/bioseeker)
[![CI](https://github.com/SouthernBio/BioSeeker/actions/workflows/ci.yml/badge.svg)](https://github.com/SouthernBio/BioSeeker/actions)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.md)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](pyproject.toml)

![BioSeeker](https://github.com/SouthernBio/BioSeeker/raw/main/logo.png)

BioSeeker is a Python library and CLI tool for the analysis of codon and bicodon conservation rates across linked species.

With this project, you can calculate codon and bicodon conservation rates for a given genus. It is based on the math from [Conservation of location of several specific inhibitory codon pairs in the Saccharomyces sensu stricto yeasts reveals translational selection (Ghoneim et al., 2019)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6379720/).

## What it does

BioSeeker calculates conservation rates for single codons and codon pairs (bicodons) across three reading frames (ORF 0, 1, 2) by reading standard FASTA formats with Biopython. It automatically aggregates results from multiple alignment files and uses NumPy to speed up the math.

Beyond simple counting, it runs linear regressions, finds outliers, and computes dipeptide-based z-scores to help you isolate translational selection signals. It can also calculate the Relative Synonymous Codon Usage (RSCU) for your dataset. If you prefer to avoid installing local dependencies, there is a Docker image available.

You can read the exact mathematical definitions of these metrics in the [Methodology Documentation](https://github.com/SouthernBio/BioSeeker/blob/main/docs/METHODOLOGY.md).

## Installation

### From Ubuntu PPA (Recommended for Ubuntu/Debian)

```bash
sudo add-apt-repository ppa:fx-biocoder/bioseeker
sudo apt update
sudo apt install bioseeker
```

### From PyPI

```bash
pip install bioseeker
```

### From Source

Requires Python 3.11+.

```bash
git clone https://github.com/SouthernBio/BioSeeker.git
cd BioSeeker
pip install .
```

### Development

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pytest
```

## Usage

### Command Line

BioSeeker provides a command line interface with a few subcommands. Place your aligned FASTA files in a directory (e.g., `FASTA_files`).

**Full Analysis Pipeline:**
The `analyze` subcommand runs the full pipeline (counts, rates, regressions, and plots).

```bash
bioseeker analyze --input FASTA_files --output results --species Scer Spar Smik Skud Sbay
```

**Key Arguments for `analyze`:**
- `--input`, `-i`: Directory containing input alignment files (default: `FASTA_files`).
- `--output`, `-o`: Directory to save CSV results and plots (default: `results`).
- `--species`: Only process alignments containing all specified species (e.g., `--species Scer Spar Smik Skud Sbay`). Incomplete alignments are discarded.
- `--reference`: Reference sequence ID (default: first sequence in the file).
- `--threshold`: Conservation threshold (default: `1.0`, meaning 100% identity required).
- `--genetic-code`: NCBI translation table ID (default: `1` for Standard). This is necessary because different organisms may use different translation tables (e.g., *Candida* uses table 12 where CTG translates to Serine instead of Leucine, and different tables have different stop codons).
- `--orfs`: Reading frames to analyze (default: `0 1 2`).
- `--sigma-threshold`: Threshold (in σ) for outlier detection (default: `3.0`).

### RSCU Analysis

To calculate the global RSCU values across your alignments:

```bash
bioseeker rscu --input FASTA_files --output results
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

The tool generates CSV files and plots in the output directory for each analyzed ORF (e.g., `ORF0`):
- `codon_conservation_ORF0.csv`: Raw codon counts.
- `bicodon_conservation_ORF0.csv`: Raw bicodon counts.
- `bicodon_stats_ORF0.csv`: Full metrics including ExpectedRate, FittedRate, Outlier status, NormalizedScore, and DipeptideZ.
- `conservation_regression_ORF0.png`: Plot of observed vs expected conservation rates.

For a full explanation of the outputs, refer to the [Methodology Documentation](https://github.com/SouthernBio/BioSeeker/blob/main/docs/METHODOLOGY.md).

## License

GNU General Public License v3.0. See [LICENSE](LICENSE.md).

## 💙 Support this project

Your contribution would help SouthernBio in improving the quality of this project and adding additional features. If you find this project useful and/or interesting, please consider offering your support on Github Sponsors, Ko-Fi or PayPal

[![Github-sponsors](https://img.shields.io/badge/sponsor-30363D?style=for-the-badge&logo=GitHub-Sponsors&logoColor=#EA4AAA)](https://github.com/sponsors/fx-biocoder) [![Ko-Fi](https://img.shields.io/badge/Ko--fi-F16061?style=for-the-badge&logo=ko-fi&logoColor=white)](https://ko-fi.com/biocoder) [![PayPal](https://img.shields.io/badge/PayPal-00457C?style=for-the-badge&logo=paypal&logoColor=white)](https://paypal.me/facumartinez680)
