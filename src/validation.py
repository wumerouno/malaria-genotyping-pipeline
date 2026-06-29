import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pandas as pd
from src.utils import setup_logging

logger = setup_logging()

# Define official benchmarks
BENCHMARKS = {
    'total_genotyped': 267,
    'total_deletion_tested': 101,  # Genotyped + Deletion tested for 4 sites (Kano, Gombe, Ogbomosho, Bayelsa)
    'evaluable_deletions': 115,    # 119 tested - 4 indeterminate (overall)
    'pfhrp2_deletions': 59,        # 59 deletions out of 115
    'deletion_prevalence': 51.3,   # 59/115 = 51.3%
}

def validate_pipeline_results(conformed_df, genotyping_df, deletions_df, accuracy_df, diversity_df):
    """
    Validates that the pipeline's output datasets match published dissertation numbers.
    Any discrepancy is logged as an error (defect).
    """
    logger.info("Starting validation checks against published benchmarks...")
    
    defects = []
    
    # 1. Total Genotyped Isolates
    n_geno = len(genotyping_df)
    if n_geno != BENCHMARKS['total_genotyped']:
        defects.append(f"Genotyped cohort count mismatch: got {n_geno}, expected {BENCHMARKS['total_genotyped']}")
    else:
        logger.info(f"SUCCESS: Total genotyped isolates is exactly {n_geno}.")
        
    # 2. Deletion Cohort Size for Linked Sites
    # In conformed deletions sheet, there are exactly 101 samples from Kano, Gombe, Ogbomosho, Bayelsa
    n_del_linked = deletions_df[deletions_df['Site'].isin(['Kano', 'Gombe', 'Ogbomosho', 'Bayelsa'])].shape[0]
    if n_del_linked != BENCHMARKS['total_deletion_tested']:
        defects.append(f"Linked deletion cohort count mismatch: got {n_del_linked}, expected {BENCHMARKS['total_deletion_tested']}")
    else:
        logger.info(f"SUCCESS: Total linked deletion-tested isolates is exactly {n_del_linked}.")
        
    # 3. Overall Deletion Prevalence
    # For overall deletion, we include Owerri (18 samples).
    # Evaluable samples = 119 - 4 indeterminate = 115.
    # We filter out indeterminate samples (where HRP2 or HRP3 is NaN)
    eval_df = deletions_df[deletions_df['HRP2'].isin([0.0, 1.0]) & deletions_df['HRP3'].isin([0.0, 1.0])]
    n_eval = len(eval_df)
    n_del = eval_df[eval_df['HRP2'] == 0.0].shape[0]
    pct_del = (n_del / n_eval) * 100 if n_eval > 0 else 0.0
    
    if n_eval != BENCHMARKS['evaluable_deletions']:
        defects.append(f"Evaluable deletions count mismatch: got {n_eval}, expected {BENCHMARKS['evaluable_deletions']}")
    else:
        logger.info(f"SUCCESS: Total evaluable deletion-tested isolates is exactly {n_eval}.")
        
    if n_del != BENCHMARKS['pfhrp2_deletions']:
        defects.append(f"Pfhrp2 deletion count mismatch: got {n_del}, expected {BENCHMARKS['pfhrp2_deletions']}")
    else:
        logger.info(f"SUCCESS: Total pfhrp2 deletions is exactly {n_del}.")
        
    if abs(pct_del - BENCHMARKS['deletion_prevalence']) > 0.1:
        defects.append(f"Pfhrp2 deletion prevalence mismatch: got {pct_del:.2f}%, expected {BENCHMARKS['deletion_prevalence']}%")
    else:
        logger.info(f"SUCCESS: Pfhrp2 deletion prevalence is exactly {pct_del:.1f}%.")
        
    # 4. Denominator Exclusions and Diagnostic Accuracy Sensitivity Range
    # Expected sensitivity range: 22.6% (Gombe) to 92.3% (Yobe)
    min_sens = accuracy_df['Sensitivity'].min() * 100
    max_sens = accuracy_df['Sensitivity'].max() * 100
    
    if abs(min_sens - 22.6) > 0.1:
        defects.append(f"Minimum RDT sensitivity mismatch: got {min_sens:.1f}%, expected 22.6%")
    else:
        logger.info(f"SUCCESS: Minimum RDT sensitivity is exactly {min_sens:.1f}%.")
        
    if abs(max_sens - 92.3) > 0.1:
        defects.append(f"Maximum RDT sensitivity mismatch: got {max_sens:.1f}%, expected 92.3%")
    else:
        logger.info(f"SUCCESS: Maximum RDT sensitivity is exactly {max_sens:.1f}%.")
        
    # Print results summary
    if defects:
        logger.error(f"Validation FAILED with {len(defects)} defects:")
        for d in defects:
            logger.error(f"  - {d}")
        raise ValueError("Pipeline validation failed. See log for details.")
    else:
        logger.info("Validation PASSED! All pipeline outputs match published results perfectly.")
        
    return defects
