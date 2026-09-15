import os
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Dict, Any

from src.io_handler import stream_genomic_records
from src.filter_engine import trim_low_quality_tails, calculate_aether_purity_score
from src.variant_caller import load_clinical_db

def test_stream_genomic_records_file_not_found():
    """Verifies that attempting to stream a missing file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        list(stream_genomic_records("non_existent_file.fastq"))

def test_sliding_window_trimming():
    """Tests 3' sliding window trimming on a synthetic SeqRecord with a low-quality tail."""
    # 40 high-quality bases (Q30) followed by 10 low-quality bases (Q10)
    seq = Seq("A" * 40 + "C" * 8)
    qual = [30] * 40 + [10] * 8
    record = SeqRecord(seq, id="test_read", letter_annotations={"phred_quality": qual})
    
    # Trim with window_size=4, min_quality=20.0, min_length=35
    trimmed = trim_low_quality_tails(record, window_size=4, min_quality=20.0, min_length=35)
    
    assert len(trimmed) == 40
    assert str(trimmed.seq) == "A" * 40
    assert len(trimmed.letter_annotations["phred_quality"]) == 40

def test_aether_purity_score_calculation():
    """Tests Aether Purity Score boundary cases and expected output."""
    # Perfect score: Q40 and 0% N
    seq_perfect = Seq("A" * 50)
    qual_perfect = [40] * 50
    record_perfect = SeqRecord(seq_perfect, id="perfect", letter_annotations={"phred_quality": qual_perfect})
    assert calculate_aether_purity_score(record_perfect) == 100.0

    # Mixed score: Q40 and 10% N -> (40/40) * (1 - 0.1) * 100 = 90.0
    seq_mixed = Seq("A" * 45 + "N" * 5)
    record_mixed = SeqRecord(seq_mixed, id="mixed", letter_annotations={"phred_quality": qual_perfect})
    assert abs(calculate_aether_purity_score(record_mixed) - 90.0) < 1e-6

    # Empty sequence -> 0.0
    record_empty = SeqRecord(Seq(""), id="empty")
    assert calculate_aether_purity_score(record_empty) == 0.0

def test_lru_cache_clinical_db(tmp_path):
    """Verifies lru_cache behavior and graceful degradation when CSV is missing."""
    # Verify cache returns the identical object reference on subsequent calls
    db1 = load_clinical_db()
    db2 = load_clinical_db()
    assert db1 is db2, "lru_cache failed to return the cached object instance"

    # Test graceful degradation with a missing database path
    missing_db_path = str(tmp_path / "missing.csv")
    empty_db = load_clinical_db(missing_db_path)
    assert empty_db == {}, "Expected empty dictionary when database file is missing"
