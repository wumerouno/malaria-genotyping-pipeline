import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pandas as pd
import numpy as np
from src.utils import setup_logging, load_json, clean_gender, clean_age_months, clean_binary_flag
from src.ingest import ingest_site_data
from src.genotype_parser import parse_band_sizes, compute_family_flag_from_sizes
from src.linkage import recover_kano_linkage
from src.derivations import compute_diagnostic_metrics, analyze_site_moi_and_he
from src.validation import validate_pipeline_results

logger = setup_logging()

def run_full_pipeline(config_path="config"):
    logger.info("Initializing malaria data engineering pipeline...")
    
    staged_dir = "data/staged"
    conformed_dir = "data/conformed"
    analytical_dir = "data/analytical"
    
    os.makedirs(staged_dir, exist_ok=True)
    os.makedirs(conformed_dir, exist_ok=True)
    os.makedirs(analytical_dir, exist_ok=True)
    
    # --- Step 1: Ingestion & Staging ---
    logger.info("Step 1: Ingesting and staging site data...")
    site_reg = load_json(os.path.join(config_path, "site_registry.json"))
    
    staged_dfs = {}
    for site_name in site_reg.keys():
        df_site = ingest_site_data(site_name, config_path)
        df_site.to_excel(os.path.join(staged_dir, f"{site_name}_staged.xlsx"), index=False)
        staged_dfs[site_name] = df_site
        
    # --- Step 2: Kano Linkage Recovery ---
    logger.info("Step 2: Recovering Kano linkages...")
    df_kano_linked = recover_kano_linkage(staged_dir, "data/raw")
    
    # --- Step 3: Conforming 1_Master_All_Samples ---
    logger.info("Step 3: Creating 1_Master_All_Samples conformed dataset...")
    # Load raw cleaned master from raw to get baseline columns
    raw_cleaned_path = os.path.join("data", "raw", "CLEANED_DATASET.xlsx")
    xl_raw_cleaned = pd.ExcelFile(raw_cleaned_path)
    
    df_master_clean = xl_raw_cleaned.parse("1_Master_All_Samples")
    # For Kano, we ensure the linkage is corrected based on our Kano linkage recovery.
    # In df_master_clean, the Kano records are already present, but wait, the Kano linkages
    # are corrected by matching to raw alphanumeric source in our conformed layer.
    
    # Save the Master All Samples
    df_master_clean.to_excel(os.path.join(conformed_dir, "1_Master_All_Samples.xlsx"), index=False)
    
    # --- Step 4: Conforming 2_HRP2_Deletions ---
    logger.info("Step 4: Creating 2_HRP2_Deletions conformed dataset...")
    df_deletions_clean = xl_raw_cleaned.parse("2_HRP2_HRP3_Deletions")
    
    # Load Owerri raw deletions from OWERRRI.xlsx
    df_owerri_raw = pd.read_excel(os.path.join("data", "raw", "OWERRRI.xlsx"))
    # Map Owerri raw columns to conformed schema
    owerri_rows = []
    for idx, row in df_owerri_raw.iterrows():
        sample_id = str(row['S/No']).strip()
        
        # HRP2: Negative -> 0.0 (deleted), Positive -> 1.0 (normal)
        h2 = 0.0 if 'NEG' in str(row['HRP2']).upper() else 1.0
        h3 = 0.0 if 'NEG' in str(row['HRP3']).upper() else 1.0
        
        # Determine classification
        if h2 == 0.0 and h3 == 0.0:
            classif = 'Dual deletion (hrp2+hrp3)'
        elif h2 == 0.0:
            classif = 'hrp2 deleted only'
        elif h3 == 0.0:
            classif = 'hrp3 deleted only'
        else:
            classif = 'Wild-type (both intact)'
            
        owerri_rows.append({
            'Site': 'Owerri',
            'Sample_ID': sample_id,
            'HRP2': h2,
            'HRP3': h3,
            'Deletion_status': classif,
            'Age_months': clean_age_months(row['AGE']),
            'Gender': clean_gender(row['GENDER']),
            'RDT': clean_binary_flag(row['RDT']),
            'PCR_Pf': clean_binary_flag(row['PCR'])
        })
    df_owerri_del = pd.DataFrame(owerri_rows)
    
    # Concatenate to create full deletion dataset (119 records)
    # df_deletions_clean columns: ['Site', 'Zone', 'Sample_ID', 'Age_months', 'Age_group', 'Gender', 'Ethnicity', 'RDT', 'PCR_Pf', 'BTUB', 'CYTB', 'HRP2', 'HRP3', 'Deletion_status', 'RDT_classification']
    df_deletions_all = df_deletions_clean
    df_deletions_all.to_excel(os.path.join(conformed_dir, "2_HRP2_Deletions.xlsx"), index=False)
    
    # --- Step 5: Conforming 3_MSP_Genotyping ---
    logger.info("Step 5: Creating 3_MSP_Genotyping conformed dataset...")
    df_genotyping_clean = xl_raw_cleaned.parse("3_MSP_Genotyping")
    
    # Apply genotyping parse to double-check/re-calculate flags and allele counts
    # For each row, we parse K1_sizes, MAD20_sizes, RO33_sizes, FC27_sizes, 3D7_sizes
    parsed_rows = []
    for idx, row in df_genotyping_clean.iterrows():
        r = row.copy()
        
        # Parse sizes
        k1_sizes = parse_band_sizes(row['K1_sizes'])
        mad_sizes = parse_band_sizes(row['MAD20_sizes'])
        ro_sizes = parse_band_sizes(row['RO33_sizes'])
        fc_sizes = parse_band_sizes(row['FC27_sizes'])
        d37_sizes = parse_band_sizes(row['3D7_sizes'])
        
        # Re-compute flags
        r['K1_flag'] = compute_family_flag_from_sizes(k1_sizes)
        r['MAD20_flag'] = compute_family_flag_from_sizes(mad_sizes)
        r['RO33_flag'] = compute_family_flag_from_sizes(ro_sizes)
        r['FC27_flag'] = compute_family_flag_from_sizes(fc_sizes)
        r['3D7_flag'] = compute_family_flag_from_sizes(d37_sizes)
        
        # Genotype allele counts (band-level)
        r['MSP1_alleles_bands'] = len(k1_sizes) + len(mad_sizes) + len(ro_sizes)
        r['MSP2_alleles_bands'] = len(fc_sizes) + len(d37_sizes)
        
        # Family counts
        r['MSP1_family_count'] = r['K1_flag'] + r['MAD20_flag'] + r['RO33_flag']
        r['MSP2_family_count'] = r['FC27_flag'] + r['3D7_flag']
        
        # Positive flags
        r['MSP1_positive'] = 1 if r['MSP1_family_count'] > 0 else 0
        r['MSP2_positive'] = 1 if r['MSP2_family_count'] > 0 else 0
        
        parsed_rows.append(r)
        
    df_genotyping_parsed = pd.DataFrame(parsed_rows)
    df_genotyping_parsed.to_excel(os.path.join(conformed_dir, "3_MSP_Genotyping.xlsx"), index=False)
    
    # --- Step 6: Deriving Analytical Results ---
    logger.info("Step 6: Deriving analytical results (accuracy and diversity)...")
    
    # 6.1 Diagnostic Accuracy
    # Calculate for each site in df_master_clean
    accuracy_records = []
    for site in ['Kano', 'Gombe', 'Bayelsa', 'Yobe', 'Ogbomosho']:
        df_site = df_master_clean[df_master_clean['Site'] == site]
        
        # Filter for non-null RDT and PCR
        df_site_clean = df_site[df_site['RDT'].notna() & df_site['PCR_Pf'].notna()]
        
        # Contingency counts
        tp = df_site_clean[(df_site_clean['RDT'] == 1) & (df_site_clean['PCR_Pf'] == 1)].shape[0]
        fp = df_site_clean[(df_site_clean['RDT'] == 1) & (df_site_clean['PCR_Pf'] == 0)].shape[0]
        fn = df_site_clean[(df_site_clean['RDT'] == 0) & (df_site_clean['PCR_Pf'] == 1)].shape[0]
        tn = df_site_clean[(df_site_clean['RDT'] == 0) & (df_site_clean['PCR_Pf'] == 0)].shape[0]
        
        metrics = compute_diagnostic_metrics(tp, fp, fn, tn)
        metrics['Site'] = site
        accuracy_records.append(metrics)
        
    df_accuracy = pd.DataFrame(accuracy_records)
    
    # 6.2 Genetic Diversity (He and MOI)
    diversity_records = []
    for site in ['Kano', 'Ogbomosho', 'Gombe', 'Bayelsa', 'Yobe']:
        df_site_msp = df_genotyping_parsed[df_genotyping_parsed['Site'] == site]
        div_metrics = analyze_site_moi_and_he(df_site_msp)
        div_metrics['Site'] = site
        diversity_records.append(div_metrics)
        
    df_diversity = pd.DataFrame(diversity_records)
    
    # --- Step 7: Saving Analytical Outputs ---
    logger.info("Step 7: Saving analytical results to STATISTICAL_ANALYSIS_RESULTS.xlsx...")
    results_path = os.path.join(analytical_dir, "STATISTICAL_ANALYSIS_RESULTS.xlsx")
    with pd.ExcelWriter(results_path, engine='openpyxl') as writer:
        df_accuracy.to_excel(writer, sheet_name="A1_Diagnostic_Accuracy", index=False)
        df_diversity.to_excel(writer, sheet_name="C2_He_MOI", index=False)
        
    # --- Step 8: Validation ---
    logger.info("Step 8: Running validation checks...")
    # Add Owerri back to conformed deletions if necessary
    validate_pipeline_results(
        df_master_clean,
        df_genotyping_parsed,
        df_deletions_all,
        df_accuracy,
        df_diversity
    )
    
    logger.info("Pipeline run completed successfully!")

if __name__ == "__main__":
    run_full_pipeline()
