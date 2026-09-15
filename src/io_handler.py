"""
Module: io_handler
Responsibility: Handle low-memory footprint I/O operations for genomic data.
"""
import gzip
from Bio import SeqIO
from typing import Generator, Any

def stream_genomic_records(file_path: str, file_format: str = "fastq") -> Generator[Any, None, None]: #
    """
    Streams genomic records from plain text or gzipped files (.gz).
    Maintains O(1) RAM footprint
    
    """
    # Dynamically selects the appropriate opening function
    open_fn = gzip.open if file_path.endswith(".gz") else open
    try:
        with open_fn(file_path, "rt") as handle:
            for record in SeqIO.parse(handle, file_format):
                yield record
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Target file not found at: {file_path}") from e
    except Exception as e:
        raise RuntimeError(f"Error streaming biological data: {str(e)}") from e
    
""" si quisieramos usar una lista se veria asi: 

# VERSIÓN CON LISTA (Carga todo de golpe en RAM)
def read_all_genomic_records(file_path: str, file_format: str = "fastq"):
    # SeqIO.parse sigue leyendo, pero list() lo obliga a meter TODO en la RAM ya mismo
    records_list = list(SeqIO.parse(file_path, file_format))
    return records_list  # Devuelve la lista completa y termina """
