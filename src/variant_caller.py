"""
Module: variant_caller
Responsibility: Identify specific clinical mutations using an external database.
"""

import os
from functools import lru_cache
from typing import Any, Dict, Optional
import pandas as pd
from Bio.SeqRecord import SeqRecord

# Define default path to the clinical metadata CSV registry
DB_PATH = os.path.join("data", "metadata", "clinical_variants.csv")


@lru_cache(maxsize=1)
def load_clinical_db(db_path: str = DB_PATH) -> Dict[str, Dict[str, Any]]:
    """
    Loads the clinical variants CSV file and converts it into an in-memory hash map.

    Decorated with @lru_cache(maxsize=1) to ensure an O(1) disk I/O profile.
    In a production-grade genomic streaming pipeline processing millions of FASTQ
    reads, querying an uncached database on each read would force repeated disk reads
    and redundant CSV parsing, introducing severe disk I/O bottlenecks that stall the
    streaming engine.

    By caching the resulting dictionary after the first read, all subsequent calls
    bypass disk I/O entirely and retrieve the in-memory lookup table in O(1) time complexity.

    Args:
        db_path: Filepath to the clinical variants CSV registry.

    Returns:
        Dict[str, Dict[str, Any]]: Mapping of mutation biomarker sequence patterns to their
        corresponding clinical metadata (mutation name, phenotype, etc.).
    """
    if not os.path.exists(db_path):
        print(f"[WARNING] Clinical database not found at {db_path}. Running without biomarkers.")
        return {}

    # Load with Pandas, index by biomarker Sequence, and transform into a lookup dictionary
    df = pd.read_csv(db_path)
    return df.set_index("Sequence").to_dict(orient="index")


def scan_for_clinical_variants(record: SeqRecord, db_path: str = DB_PATH) -> Optional[Dict[str, Any]]:
    """
    Scans a high-quality read for known clinical mutation biomarkers.

    Args:
        record: BioPython SeqRecord representing a high-quality genomic read.
        db_path: Filepath to the clinical variants CSV registry.

    Returns:
        Optional[Dict[str, Any]]: Variant finding metadata if a biomarker is detected,
        None otherwise.
    """
    sequence_str = str(record.seq).upper()

    # Retrieve cached clinical variant lookup dictionary
    clinical_db = load_clinical_db(db_path)

    # Perform linear sequence matching against biomarker keys
    for marker, info in clinical_db.items():
        if marker in sequence_str:
            return {
                "Sample_ID": record.id,
                "Variant_Found": info["Mutation"],
                "Clinical_Phenotype": info["Phenotype"],
                "Sequence_Context": marker,
            }

    return None