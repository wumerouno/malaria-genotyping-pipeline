import os
import sys
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.utils import clean_gender, clean_age_months, clean_binary_flag
from src.genotype_parser import parse_band_sizes, compute_family_flag_from_sizes
from src.derivations import compute_diagnostic_metrics, compute_expected_heterozygosity

def approx(val1, val2, tol=1e-3):
    if val1 is None or val2 is None:
        return val1 == val2
    return abs(val1 - val2) <= tol

def test_demographic_cleaning():
    # Test clean_gender
    assert clean_gender("M") == "Male"
    assert clean_gender("Female") == "Female"
    assert clean_gender("f") == "Female"
    assert clean_gender("  m  ") == "Male"
    assert clean_gender(None) is None
    
    # Test clean_age_months
    assert approx(clean_age_months("36 D"), 1.2, tol=0.1)
    assert clean_age_months("18 M") == 18.0
    assert clean_age_months("12 Y") == 144.0
    assert clean_age_months(15) == 15.0
    assert clean_age_months(None) is None
    assert clean_age_months("1.5yrs") == 18.0
    assert clean_age_months("3yrs") == 36.0

def test_binary_flag_cleaning():
    assert clean_binary_flag("Positive") == 1
    assert clean_binary_flag("Negative") == 0
    assert clean_binary_flag("pos") == 1
    assert clean_binary_flag("neg") == 0
    assert clean_binary_flag(1) == 1
    assert clean_binary_flag(0) == 0
    assert clean_binary_flag(None) is None
    assert clean_binary_flag("indeterminate") is None

def test_genotype_parser():
    # Test parse_band_sizes
    assert parse_band_sizes("180/260") == [180, 260]
    assert parse_band_sizes("300, 320") == [300, 320]
    assert parse_band_sizes("220") == [220]
    assert parse_band_sizes("Negative") == []
    assert parse_band_sizes("") == []
    assert parse_band_sizes(None) == []
    
    # Test compute_family_flag_from_sizes
    assert compute_family_flag_from_sizes([180]) == 1
    assert compute_family_flag_from_sizes([]) == 0

def test_nei_heterozygosity():
    # Test compute_expected_heterozygosity
    counts = [10, 10]  # Even proportions, n = 20
    he = compute_expected_heterozygosity(counts, 20)
    # Expected: 20/19 * (1 - (0.5^2 + 0.5^2)) = 20/19 * (1 - 0.5) = 10/19 = 0.5263
    assert approx(he, 0.5263, tol=1e-4)
    
    # Test he with all zero counts
    assert compute_expected_heterozygosity([0, 0], 0) == 0.0

def test_diagnostic_metrics():
    # Test diagnostic_metrics
    # TP=28, FP=24, FN=24, TN=23
    m = compute_diagnostic_metrics(28, 24, 24, 23)
    assert approx(m['Sensitivity'], 0.5385, tol=1e-4)
    assert approx(m['Specificity'], 0.4894, tol=1e-4)
    assert approx(m['PPV'], 0.5385, tol=1e-4)
    assert approx(m['NPV'], 0.4894, tol=1e-4)
    assert approx(m['Kappa'], 0.0278, tol=1e-4)

if __name__ == "__main__":
    print("Running automated test suite manually...")
    try:
        test_demographic_cleaning()
        print(" - test_demographic_cleaning: PASSED")
        test_binary_flag_cleaning()
        print(" - test_binary_flag_cleaning: PASSED")
        test_genotype_parser()
        print(" - test_genotype_parser: PASSED")
        test_nei_heterozygosity()
        print(" - test_nei_heterozygosity: PASSED")
        test_diagnostic_metrics()
        print(" - test_diagnostic_metrics: PASSED")
        print("\nAll tests PASSED successfully!")
    except AssertionError as e:
        print(f"\nTest FAILED: AssertionError occurred.")
        sys.exit(1)
    except Exception as e:
        print(f"\nTest FAILED with error: {e}")
        sys.exit(1)
