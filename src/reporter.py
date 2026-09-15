"""
Module: reporter
Responsibility: Standardize, normalize and export clinical findings and execution metrics using Pandas.
"""

import pandas as pd
from typing import List, Dict, Any


def generate_clinical_report(detected_variants: List[Dict[str, Any]], output_path: str) -> bool:
    """
    Takes a list of variant dictionaries, normalizes them into a DataFrame, 
    and exports a structured CSV report.
    
    Args:
        detected_variants: List of dictionaries containing variant details.
        output_path: Path where the CSV report will be saved.
        
    Returns:
        bool: True if the report was generated successfully, False if there were no variants to report.
    """
    if not detected_variants:
        return False
        
    # Use Pandas to structure the data matrix
    df = pd.DataFrame(detected_variants)
    
    # Enforce a clean, international column order for audits
    column_order = [
        "Sample_ID", 
        "Variant_Found", 
        "Clinical_Phenotype", 
        "Sequence_Context", 
        "Aether_Purity_Score"
    ]
    
    # Ensure all required columns exist even if data is missing, filling with NA
    for col in column_order:
        if col not in df.columns:
            df[col] = pd.NA
            
    df = df[column_order]
    
    # Export without the Pandas numeric index to keep the file clean
    df.to_csv(output_path, index=False)
    return True


def compute_execution_summary(
    total_reads: int, 
    passed_reads: int, 
    dropped_reads: int, 
    detected_variants_list: List[Dict[str, Any]], 
    purity_scores_list: List[float]
) -> Dict[str, Any]:
    """
    Computes normalized pipeline statistics based on raw execution metrics.
    
    Args:
        total_reads: Total number of raw reads processed.
        passed_reads: Number of reads that passed all QC filters.
        dropped_reads: Number of reads that failed QC.
        detected_variants_list: List of identified clinical variants.
        purity_scores_list: List of Aether Purity Scores for all passed reads.
        
    Returns:
        Dict[str, Any]: Normalized pipeline statistics including pass rate, drop rate,
                        and mean purity score.
    """
    pass_rate = (passed_reads / total_reads * 100.0) if total_reads > 0 else 0.0
    drop_rate = (dropped_reads / total_reads * 100.0) if total_reads > 0 else 0.0
    
    mean_purity = (
        sum(purity_scores_list) / len(purity_scores_list) 
        if purity_scores_list else 0.0
    )
    
    return {
        "Total_Reads": total_reads,
        "Passed_Reads": passed_reads,
        "Dropped_Reads": dropped_reads,
        "Pass_Rate_Percent": pass_rate,
        "Drop_Rate_Percent": drop_rate,
        "Variants_Identified": len(detected_variants_list),
        "Mean_Aether_Purity_Score": mean_purity
    }