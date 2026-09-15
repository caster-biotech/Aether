"""
Module: variant_caller
Responsibility: Identify specific clinical mutations using an external database.
"""

import os
import pandas as pd
from functools import lru_cache
from Bio.SeqRecord import SeqRecord
from typing import Dict, Optional

# Define metadata CSV route 
DB_PATH = os.path.join("data", "metadata", "clinical_variants.csv")

@lru_cache(maxsize=1)
def load_clinical_db() -> Dict:
    """Loads the clinical variants CSV and converts it into an optimized lookup dictionary.
    
    Decorated with @lru_cache(maxsize=1) to guarantee an O(1) disk I/O footprint. 
    Without caching, calling this function inside a high-throughput streaming loop 
    (e.g., millions of FASTQ reads) causes severe I/O bottlenecks by re-reading 
    the CSV file from disk on every single record. 
    
    The decorator intercepts execution after the first run and returns the in-memory 
    dictionary directly for all subsequent calls.
        
    """
    if not os.path.exists(DB_PATH):
        # If file doesn't exists, return empty dict to avoid script collapse 
        print(f"[WARNING] Clinical database not found at {DB_PATH}. Running without biomarkers.")
        return {}
    
    # Load with Pandas, index secuence and convert to dict
    df = pd.read_csv(DB_PATH)
    return df.set_index("Sequence").to_dict(orient="index")

def scan_for_clinical_variants(record: SeqRecord) -> Optional[Dict]:
    """Scans a high-quality read for known clinical mutation biomarkers."""
    sequence_str = str(record.seq).upper()
    
    # Loads the database dynamically
    clinical_db = load_clinical_db()
    
    # Efficient search in the dictionary loaded from the CSV.
    for marker, info in clinical_db.items():
        if marker in sequence_str:
            return {
                "Sample_ID": record.id,
                "Variant_Found": info["Mutation"],
                "Clinical_Phenotype": info["Phenotype"],
                "Sequence_Context": marker
            }
            
    return None