No hace falta que me pases el `README.md` anterior. Tengo el registro completo de la arquitectura de **Aether**, el flujo de streaming $O(1)$, la lógica de las tres fases, las banderas de la CLI, la suite de pruebas con `pytest`, la integración de DevOps (Docker y Snakemake) y la estrategia de *Targeted Clinical Screening*.

Aquí tienes una versión profesional, completa y estructurada del `README.md` lista para colocar en la raíz de tu repositorio.

---

### `README.md` para Aether

```markdown
# 🧬 Aether: Production-Grade Bioinformatic Engine

[![Python Version](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![Testing](https://img.shields.io/badge/QA-Pytest-brightgreen.svg)](https://docs.pytest.org/)
[![Docker](https://img.shields.io/badge/Container-Docker-2496ED.svg)](https://www.docker.com/)
[![Workflow](https://img.shields.io/badge/Orchestration-Snakemake-green.svg)](https://snakemake.github.io/)

**Aether** is a high-performance, production-grade genomic data processor designed for real-time Quality Control (QC), sliding-window read trimming, and targeted clinical variant screening. 

Built with an $O(1)$ memory streaming pipeline, Aether processes raw or `gzip`-compressed FASTQ files without uncompressing them to disk or exhausting system RAM, making it suitable for high-throughput sequencing (HTS) datasets of any scale.

---

## 🚀 Key Features

* **$O(1)$ Memory Efficiency:** Native `gzip` stream processing ensures memory footprint remains flat (~40–50 MB) even when processing multi-gigabyte FASTQ files.
* **Dynamic Quality Trimming:** 3' sliding-window algorithm that removes low-quality nucleotide tails based on Phred scores.
* **Aether Purity Score:** Custom quantitative metric evaluating sequence accuracy and ambiguous base penalty ($N$-content) per read.
* **Targeted Clinical Screening:** Deterministic sequence matching against a target database of clinical phenotypes and mutations.
* **Rich Terminal Interface:** Real-time progress indicators, formatted logging, and tabular execution summaries powered by `rich`.
* **DevOps Ready:** Built-in suite of unit tests (`pytest`), multi-stage lightweight `Dockerfile`, and `Snakemake` workflow orchestration.

---

## 🛠️ Targeted Clinical Screening Strategy

Aether uses a targeted screening model for variant detection rather than computationally heavy genome-wide alignments ($O(N)$ space/time overhead):

1. **Target Catalog:** Known mutations or sequence signatures are stored in `data/metadata/clinical_variants.csv`.
2. **In-Memory Matching:** Reads are screened in real time against the cached catalog (`lru_cache`).
3. **Automated Reporting:** Matched variants and their associated *Aether Purity Scores* are exported directly to an output CSV report.

> **Note:** A default clinical template targeting oncogenic mutations (e.g., *BRAF V600E*, *EGFR L858R*) is provided in `data/metadata/clinical_variants.csv`. Users can extend this file for custom pathogen or resistance gene panels.

---

## 📂 Project Structure

```text
Aether/
├── src/
│   ├── io_handler.py       # Streaming FASTQ reader (gzip / uncompressed)
│   ├── filter_engine.py    # Sliding window trimming & Aether Purity Score
│   ├── variant_caller.py   # Clinical DB loader & sequence context matching
│   └── reporter.py         # CSV report generator & summary table rendering
├── tests/
│   ├── conftest.py         # Pytest configuration & environment hooks
│   └── test_core.py        # Unit test suite for core logic
├── data/
│   ├── raw/                # Input FASTQ / FASTQ.GZ samples (git-ignored)
│   ├── processed/          # Generated CSV reports (git-ignored)
│   └── metadata/           # Target clinical variant database
├── main.py                 # CLI interface & orchestration entrypoint
├── Dockerfile              # Multi-stage lightweight container setup
├── Snakefile               # Snakemake workflow specification
├── requirements.txt        # Production & testing dependencies
├── conftest.py             # Root test discovery configuration
└── README.md               # Project documentation

```

---

## ⚡ Quickstart & Usage

### 1. Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/caster-biotech/Aether.git
cd Aether
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install -r requirements.txt

```

### 2. Running the Pipeline via CLI

```bash
python main.py \
    --input data/raw/sample_1.fastq.gz \
    --output data/processed/clinical_report.csv \
    --db data/metadata/clinical_variants.csv \
    --min-phred 30.0 \
    --min-length 35

```

#### CLI Options:

* `-i, --input` **(Required)**: Path to raw or compressed FASTQ file (`.fastq`, `.fastq.gz`).
* `-o, --output`: Path for the output CSV report (Default: `data/processed/clinical_report.csv`).
* `--db`: Path to the targeted clinical variants database (Default: `data/metadata/clinical_variants.csv`).
* `--min-phred`: Minimum average Phred quality threshold for sliding-window trimming (Default: `20.0`).
* `--min-length`: Minimum read length retained after trimming (Default: `35`).

---

## 🧪 Quality Assurance (QA)

Run automated unit tests covering IO handling, trimming edge cases, purity score calculations, and database cache degradation:

```bash
pytest

```

---

## 🐳 Docker Deployment

Build and run Aether inside a isolated container:

```bash
# Build Docker image
docker build -t aether:latest .

# Run processing inside container
docker run --rm \
  -v $(pwd)/data:/app/data \
  aether:latest -i data/raw/sample_1.fastq.gz -o data/processed/report.csv

```

---

## 🔀 Workflow Orchestration (Snakemake)

Execute Aether as part of a Snakemake pipeline:

```bash
snakemake --cores 1

```

