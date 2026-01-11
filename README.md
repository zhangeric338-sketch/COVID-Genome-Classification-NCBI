# COVID-Genome-Classification-NCBI

A robust pipeline for downloading, processing, and visualizing SARS-CoV-2 genome sequences from NCBI for machine learning classification tasks.

## Overview

This project provides a complete data preparation pipeline for COVID-19 genomic research. It downloads coronavirus genome sequences from the National Center for Biotechnology Information (NCBI) database, organizes them by variant strain, splits them into training/test sets with balanced representation, and generates visualizations for analysis.

**Key Use Cases:**
- Machine learning research on viral variant classification
- Genomic data analysis and visualization
- Educational projects on COVID-19 variants
- NCBI API/CLI workflow demonstrations

## Features

- **Multi-Variant Dataset** - 140 pre-selected accessions across 7 SARS-CoV-2 strains
- **Flexible Downloads** - Partial (size-based) or full dataset options
- **Parallel Processing** - Configurable worker count for concurrent downloads
- **Balanced Splitting** - Proportional train/test splits preserving variant distribution
- **Rich Visualizations** - Automated charts for dataset composition and analysis
- **Dual Download Methods** - NCBI CLI (fast) with API fallback
- **Reproducible** - Deterministic sampling with configurable random seed
- **Colab Ready** - Automatic Google Drive integration when running in Colab

## Quick Start

```bash
# Clone the repository
git clone https://github.com/zhangeric338-sketch/COVID-Genome-Classification-NCBI.git
cd COVID-Genome-Classification-NCBI

# Install dependencies
pip install -r requirements.txt

# Download a 50MB sample dataset and generate visualizations
python main.py
```

## Installation

### Requirements

- Python 3.7+
- pip

### Python Dependencies

```bash
pip install -r requirements.txt
```

This installs:
- `ncbi-datasets-pylib` - NCBI Datasets Python library
- `requests` - HTTP requests for API fallback
- `matplotlib` - Data visualization
- `pandas` - Data analysis and JSON handling

### Optional: NCBI Datasets CLI

For faster downloads, install the NCBI Datasets command-line tool:

```bash
# macOS (Homebrew)
brew install ncbi-datasets-cli

# Linux
curl -o datasets https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
chmod +x datasets
sudo mv datasets /usr/local/bin/

# Verify installation
datasets --version
```

See [NCBI Datasets CLI documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) for more options.

> **Note:** The pipeline works without the CLI by falling back to the NCBI API, but downloads will be slower.

## Usage

### Basic Commands

```bash
# Download default 50MB dataset
python main.py

# Query NCBI for available data size (no download)
python main.py --query-size

# Download full dataset (all 140 accessions, ~840MB)
python main.py --full-dataset

# Download specific size
python main.py --size-gb 0.5

# Custom output directory
python main.py --output-dir ./my_data

# Adjust parallel workers
python main.py --workers 8

# Set random seed for reproducibility
python main.py --seed 123
```

### Command-Line Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--full-dataset` | Download all 140 accessions | `False` |
| `--size-gb` | Target dataset size in GB | `0.05` (50MB) |
| `--output-dir` | Output directory path | `data` |
| `--seed` | Random seed for sampling | `42` |
| `--workers` | Parallel download workers | `4` |
| `--query-size` | Query total available size, then exit | `False` |

### Standalone Scripts

```bash
# Download only (no visualization)
python tools/download_data.py --output ./covid_data

# Visualize existing dataset
python src/visualization.py
```

## Project Structure

```
COVID-Genome-Classification-NCBI/
├── main.py                 # Main entry point - orchestrates the pipeline
├── requirements.txt        # Python dependencies
├── LICENSE                 # MIT License
├── README.md              # This file
│
├── src/
│   ├── download.py        # Basic download module (legacy)
│   └── visualization.py   # Dataset visualization and train/test splitting
│
├── tools/
│   └── download_data.py   # Advanced parallel download with 140 accessions
│
└── data/                  # Generated output directory
    ├── *.zip              # Downloaded genome files
    ├── train/             # Training set (80%)
    ├── test/              # Test set (20%)
    └── *.png              # Visualization outputs
```

## Data Pipeline

```
┌─────────────────────────────────────────────────────────────────┐
│                          main.py                                │
│         Parse arguments, orchestrate pipeline                   │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                   tools/download_data.py                        │
│  • Sample accessions based on --size-gb or use all (--full)    │
│  • Download via NCBI CLI or API fallback                       │
│  • Parallel processing with configurable workers               │
│  • Save .zip files to output directory                         │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────────┐
│                   src/visualization.py                          │
│  • Load and analyze NCBI metadata                              │
│  • Generate dataset composition charts                         │
│  • Split into train/test with balanced variants                │
│  • Create train/test distribution visualizations               │
└────────────────────────────┬────────────────────────────────────┘
                             │
                             ▼
                    Output: data/ directory
              (genomes, train/test splits, PNG charts)
```

## Dataset Details

### Variant Strains

The dataset includes 140 pre-selected SARS-CoV-2 genome accessions, 20 from each variant:

| Strain | Description | Example Accessions |
|--------|-------------|-------------------|
| **Reference** | Early pandemic genomes | MT123290, MT123291, ... |
| **Alpha** | B.1.1.7 lineage | MT188341, MT326110, ... |
| **Beta** | B.1.351 lineage | MT291826, MT291827, ... |
| **Gamma** | P.1 lineage | MW633477-MW633496 |
| **Delta** | B.1.617.2 lineage | MW633497-MW633516 |
| **Omicron** | B.1.1.529 lineage | OM095411, OM095412, ... |
| **Recent** | Latest submissions | ON563414, OP912844, ... |

### Size Estimates

| Download Option | Accessions | Estimated Size |
|-----------------|------------|----------------|
| Default (50MB) | ~8 | 50 MB |
| `--size-gb 0.5` | ~83 | 500 MB |
| `--full-dataset` | 140 | ~840 MB |

> **Note:** Individual genome files average ~6MB each.

### Data Source

- **Database:** NCBI GenBank
- **Organism:** SARS-CoV-2
- **Host:** Human
- **Quality:** Complete genomes only

## Train/Test Splitting

The pipeline automatically splits downloaded data into training and test sets:

- **Default Ratio:** 80% train / 20% test
- **Balanced Splitting:** Each variant strain is proportionally represented in both sets
- **Deterministic:** Uses random seed for reproducible splits across runs
- **Output:** Files copied to `train/` and `test/` subdirectories

### Example Distribution

For a full dataset download:
```
Train (112 genomes):
  - Reference: 16, Alpha: 16, Beta: 16, Gamma: 16
  - Delta: 16, Omicron: 16, Recent: 16

Test (28 genomes):
  - Reference: 4, Alpha: 4, Beta: 4, Gamma: 4
  - Delta: 4, Omicron: 4, Recent: 4
```

## Output Files

After running the pipeline, the output directory contains:

| File/Directory | Description |
|----------------|-------------|
| `*.zip` | Downloaded genome files (NCBI format) |
| `train/` | Training set genomes |
| `test/` | Test set genomes |
| `download_summary.png` | Overview of downloaded data |
| `train_test_split.png` | Train/test distribution chart |
| `dataset_composition.png` | Detailed 4-panel variant analysis |

## Examples

### Example 1: Quick Research Dataset

Download a small dataset for initial experimentation:

```bash
python main.py --size-gb 0.1 --output-dir ./quick_test
```

### Example 2: Full Reproducible Dataset

Download everything with a specific seed for reproducibility:

```bash
python main.py --full-dataset --seed 42 --output-dir ./full_data
```

### Example 3: Fast Download with More Workers

Maximize download speed on a good connection:

```bash
python main.py --full-dataset --workers 16
```

### Example 4: Check Available Data First

See what's available before committing to a download:

```bash
python main.py --query-size
```

## Google Colab

The pipeline automatically detects Google Colab environments and saves to Google Drive:

```python
# In a Colab notebook
!git clone https://github.com/zhangeric338-sketch/COVID-Genome-Classification-NCBI.git
%cd COVID-Genome-Classification-NCBI
!pip install -r requirements.txt
!python main.py --size-gb 0.1
```

Data will be saved to `/content/drive/MyDrive/` if Drive is mounted.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'feat: add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/) for providing genomic data access
- The global research community contributing SARS-CoV-2 sequences to GenBank
