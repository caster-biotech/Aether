"""
Main Orchestrator for the BioProcessor Core - High Quality Filter Step.
"""
"""
Main Orchestrator for the BioProcessor Core - Production Ready.
"""

import os
from src.io_handler import stream_genomic_records as stream_fastq_records
from src.filter_engine import is_high_quality
from src.variant_caller import scan_for_clinical_variants
from src.reporter import generate_clinical_report

def main():
    # Configuración de rutas de archivos
    input_fastq = os.path.join("data", "raw", "sample.fastq")
    output_report = os.path.join("data", "processed", "clinical_report.csv")
    
    # Asegurar que las carpetas de salida existan
    os.makedirs(os.path.dirname(output_report), exist_ok=True)
    os.makedirs(os.path.dirname(input_fastq), exist_ok=True)

    # =========================================================================
    # GENERADOR DE DATOS SINTÉTICOS (Para pruebas locales/Termux offline)
    # =========================================================================
    if not os.path.exists(input_fastq):
        print("[*] Input FASTQ not found. Generating synthetic clinical test reads...")
        with open(input_fastq, "w") as f:
            # Caso 1: Secuencia de alta calidad que CONTIENE la mutación "ACCGTTTGA" (S450L de Rifampicina)
            f.write("@PATIENT_001_RIF_RESISTANT\nACCGTTTGATCG\n+\nIIIIIIIIIIII\n") 
            # Caso 2: Falla por base individual mala (Tiene un '!' que equivale a Phred Q0)
            f.write("@READ_002_BAD_BASE\nATCGATCGATCG\n+\nIIII!IIIIIII\n") 
            # Caso 3: Falla por exceso de ambigüedad (Tiene demasiadas letras 'N')
            f.write("@READ_003_TOO_MANY_NS\nATNNNTNNNNCG\n+\nIIIIIIIIIIII\n")
    # =========================================================================

    print("[+] Initializing Production-Ready BioProcessor Core...")
    
    # 1. Capa de Lectura (Streaming)
    raw_reads = stream_fastq_records(input_fastq)
    detected_variants = []
    
    for read in raw_reads:
        # 2. Capa de Control de Calidad (Filtros estrictos)
        if is_high_quality(read):
            # 3. Capa de Identificación (Búsqueda en base de datos externa)
            variant_report = scan_for_clinical_variants(read)
            if variant_report:
                print(f"    [ALERT] Clinical Variant Found in {read.id}!")
                detected_variants.append(variant_report)
                
    # 4. Capa de Reporte (Exportación con Pandas)
    success = generate_clinical_report(detected_variants, output_report)
    if success:
        print(f"\n[SUCCESS] Pipeline completed. Clinical report saved at: {output_report}")
    else:
        print("\n[INFO] Pipeline completed. No critical anomalies detected.")

if __name__ == "__main__":
    main()