import re
import pandas as pd

def parse_band_sizes(val):
    """
    Parses band sizes from string or numeric input.
    E.g. '180/260' -> [180, 260]
         '180' -> [180]
         180 -> [180]
         '0' -> []
         NaN -> []
    """
    if pd.isna(val):
        return []
    
    # Clean and split string
    s = str(val).strip()
    if s == '0' or s == '0.0' or s.upper() == 'NEG' or s.upper() == 'NEGATIVE':
        return []
    
    # Replace non-digit/slash delimiters with commas
    s = s.replace('/', ',').replace(';', ',').replace('&', ',')
    
    # Extract all digits
    bands = []
    parts = s.split(',')
    for p in parts:
        p_clean = ''.join(c for c in p if c.isdigit())
        if p_clean:
            b_val = int(p_clean)
            if b_val > 0:
                bands.append(b_val)
    return bands

def compute_family_flag_from_sizes(sizes_list):
    """
    If size list is non-empty, the family is present (flag = 1).
    Otherwise flag = 0.
    """
    return 1 if len(sizes_list) > 0 else 0
