import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pandas as pd
from src.utils import setup_logging, load_json, clean_gender, clean_age_months, clean_binary_flag

logger = setup_logging()

def ingest_site_data(site_name, config_path="config"):
    # Load config files
    site_reg = load_json(os.path.join(config_path, "site_registry.json"))
    col_map = load_json(os.path.join(config_path, "column_mapping.json"))
    
    site_info = site_reg[site_name]
    mapping = col_map[site_name]
    
    raw_dir = os.path.join("data", "raw")
    filepath = os.path.join(raw_dir, site_info["raw_file"])
    
    logger.info(f"Ingesting {site_name} from {filepath}")
    
    # Load raw excel
    df = pd.read_excel(filepath, sheet_name=site_info["raw_sheet"])
    logger.info(f"Raw shape: {df.shape}")
    
    # Rename columns using mapping
    inv_mapping = {v: k for k, v in mapping.items()}
    df_renamed = df.rename(columns=inv_mapping)
    
    # Keep only conformed columns that exist in the renamed dataframe
    conformed_cols = list(mapping.keys())
    # Identify which conformed columns actually mapped
    available_cols = [c for c in conformed_cols if c in df_renamed.columns]
    df_conformed = df_renamed[available_cols].copy()
    
    # Standard cleanings
    if 'Gender' in df_conformed.columns:
        df_conformed['Gender_clean'] = df_conformed['Gender'].apply(clean_gender)
    if 'Age' in df_conformed.columns:
        df_conformed['Age_months_clean'] = df_conformed['Age'].apply(clean_age_months)
    if 'RDT' in df_conformed.columns:
        df_conformed['RDT_clean'] = df_conformed['RDT'].apply(clean_binary_flag)
    if 'PCR' in df_conformed.columns:
        df_conformed['PCR_clean'] = df_conformed['PCR'].apply(clean_binary_flag)
        
    # Site label
    df_conformed['Site'] = site_name
    
    # Exclusions
    # Note: We do NOT apply exclusions directly to the conformed dataset. We keep all rows
    # but add a flag or apply exclusions in derivations and validation to produce denominators.
    
    logger.info(f"Conformed shape: {df_conformed.shape}")
    return df_conformed

def ingest_all_sites(config_path="config"):
    site_reg = load_json(os.path.join(config_path, "site_registry.json"))
    staged_dir = os.path.join("data", "staged")
    os.makedirs(staged_dir, exist_ok=True)
    
    for site_name in site_reg.keys():
        df_site = ingest_site_data(site_name, config_path)
        output_path = os.path.join(staged_dir, f"{site_name}_staged.xlsx")
        df_site.to_excel(output_path, index=False)
        logger.info(f"Saved conformed data to {output_path}")

if __name__ == "__main__":
    ingest_all_sites()
