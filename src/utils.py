import logging
import os
import json
import re
import pandas as pd

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger("malaria_pipeline")

def load_json(filepath):
    with open(filepath, 'r') as f:
        return json.load(f)

def clean_gender(val):
    if pd.isna(val):
        return None
    s = str(val).strip().upper()
    if s.startswith('M') or s == '1' or s == '1.0':
        return 'Male'
    if s.startswith('F') or s == '0' or s == '0.0':
        return 'Female'
    return None

def clean_age_months(val):
    """
    Standardize age to months.
    Handles numeric, string like '18 M', '1 Y', '36 D', etc.
    """
    if pd.isna(val):
        return None
    s = str(val).strip().upper()
    # Check if pure numeric
    try:
        return float(s)
    except ValueError:
        pass
    
    # Try parsing patterns
    # Years: '10 Y', '10Y', '9 y'
    m_year = re.search(r'(\d+(?:\.\d+)?)\s*Y', s)
    if m_year:
        return float(m_year.group(1)) * 12.0
    
    # Months: '18 M', '18M', '11 M '
    m_month = re.search(r'(\d+(?:\.\d+)?)\s*M', s)
    if m_month:
        return float(m_month.group(1))
    
    # Days: '36 D', '39 D'
    m_day = re.search(r'(\d+(?:\.\d+)?)\s*D', s)
    if m_day:
        return float(m_day.group(1)) / 30.0
    
    # Check for simple lowercase 'y' or other format
    m_y_simple = re.search(r'(\d+(?:\.\d+)?)\s*y', str(val).strip())
    if m_y_simple:
        return float(m_y_simple.group(1)) * 12.0
        
    return None

def clean_binary_flag(val):
    if pd.isna(val):
        return None
    s = str(val).strip().upper()
    if s in ['1', '1.0', 'POSITIVE', 'POS', 'YES', 'TRUE', 'T']:
        return 1
    if s in ['0', '0.0', 'NEGATIVE', 'NEG', 'NO', 'FALSE', 'F']:
        return 0
    return None
