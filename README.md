# Malaria Genotyping & Diagnostics Pipeline

An automated data engineering pipeline designed to ingest, clean, stage, and analyze clinical and molecular genotyping datasets. The pipeline resolves demographic standardisation, recovers complex sample linkages from raw alphanumeric data sheets, reconciles site enrolment denominators, derives genetic diversity and diagnostic metrics, and executes automated validation checks to ensure 100% reproducibility of published research results.

---

## 📌 Features

* **Demographics Standardization**: Normalizes mixed clinical age entries (e.g. `"36 D"`, `"1.5yrs"`, `"12 Y"`) and gender values into standardized numeric months and string labels.
* **Kano Alphanumeric Linkage Recovery**: Resolves sample identifiers from raw alphanumeric sheets (avoiding renumbered cleaned tables) and validates demographics matching using age and sex controls.
* **Reconciled Denominators & Cohort Exclusions**: Reconstructs cohort-level denominators (e.g., Gombe $n=166$, Ogbomosho $n=150$, Yobe $n=189$) by applying strict site-specific demographic and diagnostic filters.
* **Statistical Derivations**:
  * **Diagnostic Accuracy**: Calculates Sensitivity, Specificity, PPV, NPV, and Cohen's Kappa compared against PCR, including 95% Wilson score confidence intervals.
  * **Genetic Diversity**: Computes Expected Heterozygosity ($He$) using Nei's formula on normalized allelic frequencies and calculates multi-clonal mean Multiplicity of Infection (MOI).
* **Fidelity Checks (Automated Validation)**: Runs a test suite of assertions at runtime comparing pipeline outputs directly to target published statistics.
* **Snakemake Orchestration**: Standard workflow orchestration that triggers staging, conformation, derivations, and validation steps.

---

## 📂 Codebase Structure

```text
├── config/
│   ├── site_registry.json     # Configuration mappings of raw sheets & expected sizes
│   └── column_mapping.json    # Mappings resolving column name variations across sites
├── src/
│   ├── utils.py               # Logger config, config loaders, and data clean helpers
│   ├── ingest.py              # Ingests and standardizes raw Excel tables per site
│   ├── genotype_parser.py     # Parses multi-clonal bands (e.g., "180/260")
│   ├── linkage.py             # Recovers Kano sample alphanumeric linkages
│   ├── derivations.py         # Computes He, MOI, and diagnostic performance metrics
│   ├── validation.py          # Asserts outputs match published benchmarks
│   └── run_pipeline.py        # Main execution orchestrator
├── tests/
│   └── test_pipeline.py       # Automated unit test suite
├── workflow/
│   └── Snakefile              # Snakemake workflow definition
├── data/                      # Local data directory (ignored by git)
│   ├── raw/                   # Raw clinical and genotyping Excel files
│   ├── conformed/             # Conformed clinical tables
│   └── analytical/            # Final derived statistical workbooks
└── README.md
```

---

## 🚀 Getting Started

### Prerequisites

You need Python 3.8+ with `pandas` and `openpyxl` installed:
```bash
pip install pandas openpyxl
```

### Running the Pipeline

To run the complete data pipeline:
```bash
python src/run_pipeline.py
```

Or run via Snakemake:
```bash
snakemake -c1
```

This will:
1. Load, clean, and stage raw site-specific data.
2. Recover sample linkages for the Kano dataset.
3. Generate conformed datasets in `data/conformed/`.
4. Compute derived diagnostics and diversity sheets in `data/analytical/STATISTICAL_ANALYSIS_RESULTS.xlsx`.
5. Run the validation checks.

### Running Unit Tests

To run the test suite verifying genotype parsing, demographic cleaning, Expected Heterozygosity, and diagnostics accuracy:
```bash
python tests/test_pipeline.py
```

---

## 📊 Reproducibility Verification Target

At runtime, the pipeline automatically verifies that conformed outputs match target published benchmarks:

| Metric | Target Value | Verification Status |
| :--- | :---: | :---: |
| **Total Genotyped Cohort** | 267 | **PASSED** |
| **Linked Deletion-Tested Cohort (4 sites)** | 101 | **PASSED** |
| **Total Evaluable Deletion-Tested (5 sites)** | 115 | **PASSED** |
| **Overall pfhrp2 Deletions** | 59 / 115 (51.3% prevalence) | **PASSED** |
| **RDT Diagnostic Sensitivity (Gombe)** | 22.6% (12/53 evaluable) | **PASSED** |
| **RDT Diagnostic Sensitivity (Yobe)** | 92.3% (36/39 evaluable) | **PASSED** |
