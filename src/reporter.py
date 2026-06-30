"""
Module: reporter
Responsibility: Standardize, normalize and export clinical findings using Pandas.
"""

import pandas as pd
from typing import List, Dict

def generate_clinical_report(detected_variants: List[Dict], output_path: str) -> bool:
    """
    Takes a list of variant dictionaries, normalizes them into a DataFrame 
    and exports a structured CSV report.
    """
    if not detected_variants:
        print("[INFO] No clinical variants to report. Output file skipped.")
        return False
        
    # Pandas entra en acción para estructurar la matriz de datos
    df = pd.DataFrame(detected_variants)
    
    # Forzamos un orden de columnas limpio e internacional para auditorías
    column_order = ["Sample_ID", "Variant_Found", "Clinical_Phenotype", "Sequence_Context"]
    df = df[column_order]
    
    # Exportamos sin el índice numérico de Pandas para mantener el archivo limpio
    df.to_csv(output_path, index=False)
    return True