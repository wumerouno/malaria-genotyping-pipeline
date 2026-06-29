import numpy as np
import pandas as pd
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.utils import setup_logging
from src.genotype_parser import parse_band_sizes, compute_family_flag_from_sizes

logger = setup_logging()

def wilson_score_interval(p, n, conf=0.95):
    """
    Computes Wilson score confidence interval for a proportion p with sample size n.
    """
    if n == 0:
        return 0.0, 0.0
    z = 1.96  # For 95% confidence
    denom = 1 + (z**2) / n
    center = p + (z**2) / (2 * n)
    spread = z * np.sqrt((p * (1 - p)) / n + (z**2) / (4 * n**2))
    lower = (center - spread) / denom
    upper = (center + spread) / denom
    return max(0.0, lower), min(1.0, upper)

def compute_diagnostic_metrics(tp, fp, fn, tn):
    """
    Computes sensitivity, specificity, PPV, NPV, and Cohen's Kappa,
    along with 95% Wilson score confidence intervals.
    """
    n = tp + fp + fn + tn
    
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0.0
    
    sens_ci = wilson_score_interval(sens, tp + fn)
    spec_ci = wilson_score_interval(spec, tn + fp)
    ppv_ci = wilson_score_interval(ppv, tp + fp)
    npv_ci = wilson_score_interval(npv, tn + fn)
    
    # Cohen's Kappa
    po = (tp + tn) / n if n > 0 else 0.0
    pe = ((tp + fn) * (tp + fp) + (fn + tn) * (fp + tn)) / (n**2) if n > 0 else 0.0
    kappa = (po - pe) / (1 - pe) if (1 - pe) > 0 else 0.0
    
    return {
        'n': n,
        'TP': tp, 'FP': fp, 'FN': fn, 'TN': tn,
        'Sensitivity': sens, 'Sens_Lower': sens_ci[0], 'Sens_Upper': sens_ci[1],
        'Specificity': spec, 'Spec_Lower': spec_ci[0], 'Spec_Upper': spec_ci[1],
        'PPV': ppv, 'PPV_Lower': ppv_ci[0], 'PPV_Upper': ppv_ci[1],
        'NPV': npv, 'NPV_Lower': npv_ci[0], 'NPV_Upper': npv_ci[1],
        'Kappa': kappa
    }

def compute_expected_heterozygosity(counts, n):
    """
    Computes Expected Heterozygosity (He) using Nei's formula on normalized frequencies:
    He = [n / (n - 1)] * [1 - sum(pi_norm^2)]
    """
    if n <= 1:
        return 0.0
    c_sum = sum(counts)
    if c_sum == 0:
        return 0.0
    pi_norm = [c / c_sum for c in counts]
    sum_pi_sq = sum(p**2 for p in pi_norm)
    he = (n / (n - 1)) * (1.0 - sum_pi_sq)
    return max(0.0, he)

def analyze_site_moi_and_he(df_site_msp):
    """
    Computes site-level MOI and Expected Heterozygosity (He) for MSP1 and MSP2,
    matching the published results workbook C2_He_MOI sheet.
    """
    n_total = len(df_site_msp)
    if n_total == 0:
        return {}
    
    # Sum of positive flags
    c_K1 = df_site_msp['K1_flag'].sum()
    c_MAD = df_site_msp['MAD20_flag'].sum()
    c_RO = df_site_msp['RO33_flag'].sum()
    c_FC = df_site_msp['FC27_flag'].sum()
    c_3D7 = df_site_msp['3D7_flag'].sum()
    
    # Heterozygosity (He)
    he_msp1 = compute_expected_heterozygosity([c_K1, c_MAD, c_RO], n_total)
    he_msp2 = compute_expected_heterozygosity([c_FC, c_3D7], n_total)
    
    # MOI denominators (number of positive samples for each marker)
    n_msp1_pos = df_site_msp[df_site_msp['MSP1_positive'] == 1].shape[0]
    n_msp2_pos = df_site_msp[df_site_msp['MSP2_positive'] == 1].shape[0]
    n_any_pos = df_site_msp[(df_site_msp['MSP1_positive'] == 1) | (df_site_msp['MSP2_positive'] == 1)].shape[0]
    
    # MOI values
    # MOI = Total alleles detected / Number of positive samples
    moi_msp1 = (c_K1 + c_MAD + c_RO) / n_msp1_pos if n_msp1_pos > 0 else 1.0
    moi_msp2 = (c_FC + c_3D7) / n_msp2_pos if n_msp2_pos > 0 else 1.0
    
    # Combined MOI
    # Total alleles detected across both loci / number of samples positive for either
    moi_combined = (c_K1 + c_MAD + c_RO + c_FC + c_3D7) / n_any_pos if n_any_pos > 0 else 1.0
    
    # Polyclonal %
    # Polyclonal if family flag sum > 1
    poly_msp1 = df_site_msp[df_site_msp['MSP1_family_count'] > 1].shape[0] / n_total if n_total > 0 else 0.0
    poly_msp2 = df_site_msp[df_site_msp['MSP2_family_count'] > 1].shape[0] / n_total if n_total > 0 else 0.0
    
    return {
        'n': n_total,
        'He_MSP1': he_msp1,
        'He_MSP2': he_msp2,
        'MOI_MSP1': moi_msp1,
        'MOI_MSP2': moi_msp2,
        'MOI_Combined': moi_combined,
        'MSP1_Polyclonal_Pct': poly_msp1,
        'MSP2_Polyclonal_Pct': poly_msp2
    }
