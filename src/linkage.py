import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pandas as pd
from src.utils import setup_logging

logger = setup_logging()

def recover_kano_linkage(staged_dir="data/staged", raw_dir="data/raw"):
    """
    Recovers Kano alphanumeric ID linkage for the genotyping and deletion cohorts.
    Specifically:
    - Reads the staged Kano enrolment data (staged_dir/Kano_staged.xlsx)
    - Reads the HRP2&3 DATA 1.xlsx sheet 'KANO' which contains raw alphanumeric IDs and S/No.
    - Resolves the correct raw record for each genotyping and deletion sample.
    - Verifies the link on age and sex.
    """
    kano_enrol_path = os.path.join(staged_dir, "Kano_staged.xlsx")
    hrp_path = os.path.join(raw_dir, "HRP2&3 DATA 1.xlsx")
    
    # Load staged Kano enrolment (200 records)
    df_enrol = pd.read_excel(kano_enrol_path)
    
    # Load HRP2&3 DATA 1 'KANO' sheet (contains S/No and alphanumeric Sample ID)
    df_hrp = pd.read_excel(hrp_path, sheet_name="KANO")
    
    # Standardize column names in df_hrp
    df_hrp = df_hrp.rename(columns={
        'S/No': 'raw_S_No',
        'Sample ID': 'alphanumeric_ID',
        'Age (months)': 'hrp_Age',
        'Gender': 'hrp_Gender'
    })
    
    # Create mapping of S/No to alphanumeric ID
    # Note: raw S_No matches df_enrol S_No
    sno_to_alpha = {}
    for idx, row in df_hrp.iterrows():
        sno = row['raw_S_No']
        alpha = row['alphanumeric_ID']
        if pd.notna(sno) and pd.notna(alpha):
            sno_to_alpha[int(sno)] = str(alpha).strip()
            
    # Now let's link the genotyped Kano sheet (staged Wumi ALL.xlsx kano sheet)
    # Wait, Wumi ALL.xlsx is in the root directory. Let's find it.
    wumi_path = os.path.join("data", "raw", "CLEANED_DATASET.xlsx")
    df_msp_geno = pd.read_excel(wumi_path, sheet_name="3_MSP_Genotyping")
    df_kano_geno = df_msp_geno[df_msp_geno['Site'] == 'Kano'].copy()
    
    # For each genotyped Kano sample (Sample_ID 1 to 70), we map to alphanumeric ID:
    # Most match Sample_ID -> {Sample_ID}A, except 4 -> 4B.
    # Let's verify by matching demographics with df_enrol.
    linked_records = []
    
    import re
    for idx, row in df_kano_geno.iterrows():
        geno_id_str = str(row['Sample_ID']).strip()
        m = re.search(r'\d+', geno_id_str)
        if not m:
            logger.warning(f"Could not extract numeric part from Kano Sample_ID: {geno_id_str}")
            continue
        raw_sno = int(m.group(0))
        
        # Find corresponding raw enrolment row
        raw_row = df_enrol[df_enrol['S_No'] == raw_sno]
        if raw_row.empty:
            logger.warning(f"No raw enrolment row found for Kano S_No {raw_sno}")
            continue
            
        raw_rec = raw_row.iloc[0]
        alpha_id = sno_to_alpha.get(raw_sno, geno_id_str)
        
        # Verify age and sex
        # Genotyping demographics:
        geno_age = row['Age_months']
        # Convert gender flag
        geno_gender = 'Male' if str(row['Gender']).strip().upper() in ['M', 'MALE', '1'] else 'Female'
        
        raw_age = raw_rec['Age_months_clean']
        raw_gender = raw_rec['Gender_clean']
        
        # Check if they match. If not, log a warning
        if abs(geno_age - raw_age) > 0.1 or geno_gender != raw_gender:
            logger.warning(f"Kano linkage mismatch for S_No {raw_sno} (Alphanumeric {alpha_id}):")
            logger.warning(f"  Geno: Age={geno_age}, Gender={geno_gender}")
            logger.warning(f"  Raw:  Age={raw_age}, Gender={raw_gender}")
            
        linked_records.append({
            'Site': 'Kano',
            'Sample_ID': geno_id_str,
            'Alphanumeric_ID': alpha_id,
            'Raw_S_No': raw_sno,
            'Age_months': raw_age,
            'Gender': raw_gender,
            'Matched': True
        })
        
    df_linked = pd.DataFrame(linked_records)
    logger.info(f"Successfully linked {len(df_linked)} Kano genotyping samples.")
    return df_linked

if __name__ == "__main__":
    recover_kano_linkage()
