"""
Module: filter_engine
Responsibility: Apply production-grade Quality Control (QC) filters on raw genomic reads.
"""

from Bio.SeqRecord import SeqRecord

def has_acceptable_ambiguity(record: SeqRecord, max_n_ratio: float = 0.05) -> bool:
    """Checks if the percentage of ambiguous bases ('N') is below the threshold."""
    sequence_str = str(record.seq).upper()
    if not sequence_str:
        return False
    n_count = sequence_str.count('N')
    return (n_count / len(sequence_str)) <= max_n_ratio

def passes_per_base_quality(record: SeqRecord, min_base_quality: int = 20) -> bool:
    """Ensures no single base falls below a critical technical quality score (e.g., Q20)."""
    quality_scores = record.letter_annotations.get("phred_quality", [])
    if not quality_scores:
        return False
    # Si alguna base es menor a min_base_quality, la lectura se rechaza
    return all(score >= min_base_quality for score in quality_scores)

def is_high_quality(
    record: SeqRecord, 
    min_phred_avg: float = 30.0, 
    min_base_quality: int = 20, 
    max_n_ratio: float = 0.05
) -> bool:
    """
    Orchestrates multiple QC filters on a single SeqRecord.
    Composes pure functions to guarantee clean code and easy debugging.
    """
    # Filtro 1: Calidad promedio global
    quality_scores = record.letter_annotations.get("phred_quality", [])
    if not quality_scores:
        return False
    
    avg_score = sum(quality_scores) / len(quality_scores)
    if avg_score < min_phred_avg:
        return False
        
    # Filtro 2: Calidad mínima por base
    if not passes_per_base_quality(record, min_base_quality):
        return False
        
    # Filtro 3: Control de bases ambiguas (Ns)
    if not has_acceptable_ambiguity(record, max_n_ratio):
        return False
        
    return True