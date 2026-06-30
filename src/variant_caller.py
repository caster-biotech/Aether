"""
Module: variant_caller
Responsibility: Identify specific clinical mutations using an external database.
"""

import os
import pandas as pd
from Bio.SeqRecord import SeqRecord
from typing import Dict, Optional

# Definimos la ruta al CSV de metadatos
DB_PATH = os.path.join("data", "metadata", "clinical_variants.csv")

def load_clinical_db() -> Dict:
    """Loads the clinical variants CSV and converts it into an optimized lookup dictionary."""
    if not os.path.exists(DB_PATH):
        # Si no existe el archivo, devolvemos un dict vacío para evitar que el script colapse
        print(f"[WARNING] Clinical database not found at {DB_PATH}. Running without biomarkers.")
        return {}
    
    # Tu línea senior: Cargamos con Pandas, indexamos por la secuencia y convertimos a dict
    df = pd.read_csv(DB_PATH)
    return df.set_index("Sequence").to_dict(orient="index")

def scan_for_clinical_variants(record: SeqRecord) -> Optional[Dict]:
    """Scans a high-quality read for known clinical mutation biomarkers."""
    sequence_str = str(record.seq).upper()
    
    # Cargamos la base de datos dinámicamente
    clinical_db = load_clinical_db()
    
    # Búsqueda eficiente en el diccionario cargado desde el CSV
    for marker, info in clinical_db.items():
        if marker in sequence_str:
            return {
                "Sample_ID": record.id,
                "Variant_Found": info["Mutation"],
                "Clinical_Phenotype": info["Phenotype"],
                "Sequence_Context": marker
            }
            
    return None