"""
Minimal Snakemake workflow for Aether BioProcessor Core integration.
Demonstrates Aether orchestrated as a processing node in a workflow manager.
"""

# Default target rule
rule all:
    input:
        "data/processed/clinical_report.csv"

# Rule to execute the Aether QC & Variant Caller pipeline
rule aether_qc_variant_call:
    input:
        fastq="data/raw/sample.fastq",
        db="data/metadata/clinical_variants.csv"
    output:
        report="data/processed/clinical_report.csv"
    log:
        "logs/aether_execution.log"
    shell:
        """
        python main.py \
            --input {input.fastq} \
            --db {input.db} \
            --output {output.report} \
            --min-phred 30.0 \
            --min-length 35 > {log} 2>&1
        """
