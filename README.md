# High-Throughput Genomic Variant & Clinical Quality Control Processor (BioProcessor Core)

A production-grade, memory-efficient Python pipeline designed to parse raw Next-Generation Sequencing (NGS) FASTQ data, enforce strict multi-tier quality control (QC) criteria, and map high-quality reads against a clinical biomarker registry using Pandas.

## 🧬 Architectural Overview & Design Patterns

The system implements a strict **Layered Architecture (Separation of Concerns)** using pure functions and generator patterns. This architecture guarantees an **O(1) memory footprint**, making it exceptionally optimized for restricted execution environments, low-RAM edge computing, or mobile-linux terminal infrastructures (e.g., Termux/Android environments).

## Pipeline 

[Raw FASTQ Stream] → (Layer #1: io_handler via Generators) → (Layer #2: filter_engine) ── [Drop if Q < 30 / Ns > 5%]  →

(High-Quality Reads Only) → (Layer #3: variant_caller via Lookup Hashing) → (Layer #4: reporter via Pandas) -> [Processed CSV Report]

### Core Components

1. **Layer #1: I/O Handler (`src/io_handler.py`)**: Implements streaming data processing via Python Generators (`yield`). It parses multi-gigabyte genomic sequences lazily, preventing heap overflow and RAM spikes.
2. **Layer #2: Quality Control Engine (`src/filter_engine.py`)**: Enforces three concurrent analytical filters:
   - Average Phred Score Threshold ($\ge$ Q30).
   - Strict local per-base quality verification ($\ge$ Q20).
   - Ambiguous base composition analysis ($N \le 5\%$).
3. **Layer #3: Clinical Variant Caller (`src/variant_caller.py`)**: Dynamically decouples data from logic by loading an external CSV metadata registry with **Pandas**, casting lookups into an optimized hash map for real-time sequence alignment.
4. **Layer #4: Analytical Reporter (`src/reporter.py`)**: Leverages Pandas structures to sanitize, normalize, and enforce regulatory compliance mapping formats before exporting finalized epidemiological CSV data.

## 🚀 Local Installation & Execution

### Prerequisites
- Python 3.10+
- Production Dependencies (specified in `requirements.txt`):
  ```text
  biopython==1.83
  pandas==2.2.2
