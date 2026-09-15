"""
Module: io_handler
Responsibility: Handle low-memory footprint streaming I/O operations for genomic data.
"""

import gzip
import os
from typing import Generator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

COMPRESSED_EXTENSIONS = (".fastq.gz", ".fq.gz")
UNCOMPRESSED_EXTENSIONS = (".fastq", ".fq")


def stream_genomic_records(
    file_path: str, 
    file_format: str = "fastq"
) -> Generator[SeqRecord, None, None]:
    """
    Streams genomic records lazily from compressed (.fastq.gz, .fq.gz) or uncompressed
    (.fastq, .fq) FASTQ files using Python's native gzip module and Biopython SeqIO.

    Maintains a strict O(1) memory (RAM) footprint by yielding one record at a time
    via generator evaluation, preventing memory exhaustion when processing large-scale
    sequencing datasets.

    Args:
        file_path: Path to the target genomic FASTQ file.
        file_format: File format expected by Biopython SeqIO (default: 'fastq').

    Yields:
        SeqRecord: Individual genomic records parsed sequentially.

    Raises:
        FileNotFoundError: If the specified file does not exist on disk.
        ValueError: If the file extension is not recognized as a valid FASTQ format.
        RuntimeError: If an error occurs during file streaming or record parsing.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Target file not found at: {file_path}")

    lower_path = file_path.lower()
    if lower_path.endswith(COMPRESSED_EXTENSIONS):
        open_fn = gzip.open
    elif lower_path.endswith(UNCOMPRESSED_EXTENSIONS):
        open_fn = open
    else:
        raise ValueError(
            f"Unsupported file format for '{file_path}'. "
            f"Supported extensions: {COMPRESSED_EXTENSIONS + UNCOMPRESSED_EXTENSIONS}"
        )

    try:
        with open_fn(file_path, mode="rt", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, file_format):
                yield record
    except Exception as e:
        raise RuntimeError(f"Error streaming biological data from '{file_path}': {str(e)}") from e
