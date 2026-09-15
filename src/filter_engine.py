"""
Module: filter_engine
Responsibility: Apply production-grade Quality Control (QC) filters and trimming on raw genomic reads.
"""

from Bio.SeqRecord import SeqRecord


def has_acceptable_ambiguity(record: SeqRecord, max_n_ratio: float = 0.05) -> bool:
    """
    Checks if the proportion of ambiguous nucleotide bases ('N') is within the allowed threshold.

    Args:
        record: BioPython SeqRecord to evaluate.
        max_n_ratio: Maximum acceptable ratio of 'N' bases (default: 0.05).

    Returns:
        bool: True if the ratio of 'N' bases does not exceed max_n_ratio, False otherwise.
    """
    sequence_str = str(record.seq).upper()
    if not sequence_str:
        return False
    n_count = sequence_str.count("N")
    return (n_count / len(sequence_str)) <= max_n_ratio


def passes_per_base_quality(record: SeqRecord, min_base_quality: int = 20) -> bool:
    """
    Ensures no individual nucleotide base falls below a critical technical quality threshold (e.g., Q20).

    Args:
        record: BioPython SeqRecord containing Phred quality scores.
        min_base_quality: Minimum allowable quality score per base (default: 20).

    Returns:
        bool: True if every base score meets or exceeds min_base_quality, False otherwise.
    """
    quality_scores = record.letter_annotations.get("phred_quality", [])
    if not quality_scores:
        return False
    # If any base quality is below min_base_quality, the read is rejected
    return all(score >= min_base_quality for score in quality_scores)


def trim_low_quality_tails(
    record: SeqRecord,
    window_size: int = 4,
    min_quality: float = 20.0,
    min_length: int = 35,
) -> SeqRecord:
    """
    Trims low-quality bases from the 3' end of a SeqRecord using a sliding window approach.

    Scans the read starting from the 3' terminal toward the 5' end using a window
    of size `window_size`. If the average Phred quality of the window falls below
    `min_quality`, the window bases are trimmed. This process continues iteratively
    until the current terminal window satisfies QC or the remaining sequence length
    drops below `min_length`.

    Args:
        record: BioPython SeqRecord containing Phred quality scores.
        window_size: Number of bases evaluated in each sliding window (default: 4).
        min_quality: Minimum mean Phred quality threshold required for the window (default: 20.0).
        min_length: Minimum allowed sequence length threshold (default: 35).

    Returns:
        SeqRecord: The trimmed SeqRecord preserving associated annotations, or the original
        record if no trimming was required.
    """
    quality_scores = record.letter_annotations.get("phred_quality", [])
    if not quality_scores:
        return record

    current_end = len(quality_scores)
    while current_end >= window_size and current_end >= min_length:
        window = quality_scores[current_end - window_size : current_end]
        avg_quality = sum(window) / float(window_size)
        if avg_quality < min_quality:
            current_end -= window_size
            if current_end < min_length:
                break
        else:
            break

    if current_end < len(record):
        return record[:max(0, current_end)]
    return record


def calculate_aether_purity_score(record: SeqRecord) -> float:
    """
    Calculates the normalized Aether Purity Score for a genomic record.

    Computes a normalized float score in the range [0.0, 100.0] derived from
    the mean Phred quality score and the proportion of ambiguous nucleotide bases ('N').
    Formula: (Mean_Phred / 40.0) * (1 - N_ratio) * 100

    Args:
        record: BioPython SeqRecord to evaluate.

    Returns:
        float: Normalized purity score between 0.0 and 100.0.
    """
    sequence_str = str(record.seq).upper()
    if not sequence_str:
        return 0.0

    quality_scores = record.letter_annotations.get("phred_quality", [])
    if not quality_scores:
        return 0.0

    mean_phred = sum(quality_scores) / len(quality_scores)
    n_ratio = sequence_str.count("N") / len(sequence_str)

    score = (mean_phred / 40.0) * (1.0 - n_ratio) * 100.0
    return max(0.0, min(100.0, float(score)))


def is_high_quality(
    record: SeqRecord,
    min_phred_avg: float = 30.0,
    min_base_quality: int = 20,
    max_n_ratio: float = 0.05,
    min_length: int = 35,
    window_size: int = 4,
    trim_min_quality: float = 20.0,
) -> bool:
    """
    Orchestrates production-grade quality control filters on an individual SeqRecord.

    Applies 3' sliding window trimming first, followed by strict validation of:
      1. Minimum post-trimming sequence length.
      2. Global average Phred quality score.
      3. Minimum per-base technical quality threshold.
      4. Ambiguous nucleotide ('N') composition ratio.

    Args:
        record: Raw input SeqRecord to filter.
        min_phred_avg: Minimum mean Phred quality score threshold (default: 30.0).
        min_base_quality: Minimum allowable quality score for each base (default: 20).
        max_n_ratio: Maximum allowable proportion of ambiguous 'N' bases (default: 0.05).
        min_length: Minimum sequence length required after trimming (default: 35).
        window_size: Sliding window size for 3' tail trimming (default: 4).
        trim_min_quality: Minimum average Phred quality required during trimming (default: 20.0).

    Returns:
        bool: True if the read passes all QC criteria, False otherwise.
    """
    # Filter 1: Apply 3' sliding window trimming
    trimmed_record = trim_low_quality_tails(
        record,
        window_size=window_size,
        min_quality=trim_min_quality,
        min_length=min_length,
    )

    # Filter 2: Verify post-trimming sequence length
    if len(trimmed_record) < min_length:
        return False

    # Filter 3: Verify global mean Phred quality
    quality_scores = trimmed_record.letter_annotations.get("phred_quality", [])
    if not quality_scores:
        return False

    avg_score = sum(quality_scores) / len(quality_scores)
    if avg_score < min_phred_avg:
        return False

    # Filter 4: Verify per-base quality threshold
    if not passes_per_base_quality(trimmed_record, min_base_quality):
        return False

    # Filter 5: Verify ambiguous nucleotide ('N') ratio
    if not has_acceptable_ambiguity(trimmed_record, max_n_ratio):
        return False

    return True