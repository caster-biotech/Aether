"""
Module: io_handler
Responsibility: Handle low-memory footprint I/O operations for genomic data.
"""

from Bio import SeqIO
from typing import Generator, Any

def stream_genomic_records(file_path: str, file_format: str = "fastq") -> Generator[Any, None, None]: #
    """
    Yields single records from a genomic file. 
    Defaults to 'fastq' format.

    Streams genomic records one by one using a generator to optimize memory allocation.
    Safe for low-RAM environments (e.g., Termux/Mobile environments).
    
    Args:
        file_path (str): Path to the raw genomic file.
        file_format (str): Biopython standard format identifier (e.g., 'fastq', 'fasta').
        
    Yields:
        Bio.SeqRecord: A single genomic record object.
    """
    try:
        # SeqIO.parse returns an iterator, keeping memory usage at O(1)
        for record in SeqIO.parse(file_path, file_format):
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
