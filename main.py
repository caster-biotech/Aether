"""
Main Orchestrator for the Aether - Production Ready.

"""
import os
import argparse
from src.io_handler import stream_genomic_records
from src.filter_engine import is_high_quality
from src.variant_caller import scan_for_clinical_variants
from src.reporter import generate_clinical_report

def main():
    parser = argparse.ArgumentParser(description="BioProcessor Core - High-Throughput QC & Variant Caller")
    parser.add_argument("--input", required=True, help="Path to input FASTQ or FASTQ.GZ file")
    parser.add_argument("--output", default=os.path.join("data", "processed", "clinical_report.csv"), help="Path for generated CSV report")
    
    args = parser.parse_args()

    # Asegurar que la carpeta de destino exista
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    print("[+] Initializing Production-Ready BioProcessor Core...")
    
    raw_reads = stream_genomic_records(args.input)
    detected_variants = []
    
    for read in raw_reads:
        if is_high_quality(read):
            variant_report = scan_for_clinical_variants(read)
            if variant_report:
                print(f"    [ALERT] Clinical Variant Found in {read.id}!")
                detected_variants.append(variant_report)
                
    success = generate_clinical_report(detected_variants, args.output)
    if success:
        print(f"\n[SUCCESS] Pipeline completed. Clinical report saved at: {args.output}")
    else:
        print("\n[INFO] Pipeline completed. No critical anomalies detected.")

if __name__ == "__main__":
    main()