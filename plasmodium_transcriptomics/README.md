# FalciTranscriptoPipe: P. falciparum RNA-seq Analysis Pipeline

FalciTranscriptoPipe is a robust, reproducible, and fully-automated bioinformatics pipeline designed to analyze *Plasmodium falciparum* transcriptomics (RNA-seq) data. The pipeline automates the downloading of public SRA data, quality control checking, decoy-aware transcriptome indexing, quantification, and differential gene expression analysis using DESeq2 in R.

For a concrete biological demonstration, this pipeline is configured to compare **Ring stage** vs. **Trophozoite stage** asexual parasites from public dataset **GSE273005** (Illumina paired-end reads).

---

## 📂 Repository Structure

```text
plasmodium_transcriptomics/
├── README.md                     # This documentation file
├── environment.yml               # Conda environment definition for tools and R packages
├── config/
│   ├── config.yaml               # Snakemake configuration (SRA runs, references, thresholds)
│   └── samples.tsv               # Metadata mapping run IDs to biological stages (Ring vs Trophozoite)
├── workflow/
│   ├── Snakefile                 # Snakemake workflow manager orchestrating the execution steps
│   └── scripts/
│       ├── download_references.py # Helper script to fetch references from PlasmoDB
│       └── create_decoy_list.py   # Prep script to build Salmon decoy-aware indexes
├── src/
│   └── deseq2_analysis.R         # R script executing DESeq2, generating tables & plots
├── notebooks/
│   └── rna_seq_report.Rmd        # R Markdown report template compiling results
└── tests/
    └── test_pipeline.py           # Unit tests validating file exists and formatting
```

---

## 🚀 Getting Started

### 1. Prerequisites (Conda / WSL)
Most command-line bioinformatics tools (`salmon`, `fastqc`, `fasterq-dump`) are native to Linux/macOS. 
* **Windows Users:** It is highly recommended to run this pipeline within **WSL2** (Windows Subsystem for Linux) or a Linux Docker container.
* **Conda/Mamba:** Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge) to manage dependencies.

### 2. Environment Installation
Create the conda environment containing all required tools (Snakemake, Salmon, FastQC, MultiQC, R, and DESeq2):

```bash
# Create the environment from the yml configuration
conda env create -f environment.yml

# Activate the environment
conda activate falcitranscripto
```

---

## 🏃 Running the Pipeline

### Step 1: Run Structure Validation Tests
Verify that all scripts are in place and that the configuration files are formatted correctly:
```bash
python tests/test_pipeline.py
```

### Step 2: Run Snakemake Dry-run
Verify the execution DAG (directed acyclic graph) without running any commands:
```bash
snakemake -s workflow/Snakefile -n
```

### Step 3: Execute the Entire Pipeline
Run the pipeline using 4 CPU cores:
```bash
snakemake -s workflow/Snakefile --cores 4
```
This single command will:
1. Download reference FASTA genomes and GFF3 files from PlasmoDB.
2. Build a decoy-aware transcript index using Salmon.
3. Download the raw FASTQ reads from NCBI SRA.
4. Run Quality Control (FastQC) and aggregate results into a MultiQC report.
5. Quantify transcript expression levels using Salmon.
6. Load abundances into R using `tximport` and run differential expression (DESeq2).
7. Generate volcano plots, PCA plots, and heatmaps of top DEGs.

---

## 📊 R Analysis & HTML Reporting

After running the pipeline, you can generate an interactive HTML report compile of all results and plots by rendering the R Markdown template.

### Run R Analysis Directly (Optional)
If you want to modify parameters or rerun the R analysis without Snakemake:
```bash
Rscript src/deseq2_analysis.R
```

### Render the HTML Report
To render the report, run the following in R:
```R
rmarkdown::render("notebooks/rna_seq_report.Rmd")
```
This compiles the report into `notebooks/rna_seq_report.html`, containing interactive tables of differentially expressed genes and high-resolution visualizations.

---

## 🧬 Biological Insights: Ring vs. Trophozoite Stage

* **Ring Stage (0–16 hours post-invasion):** Early stage characterized by a low metabolic rate. The parasite is remodeling the host red blood cell and establishing export pathways.
* **Trophozoite Stage (16–36 hours post-invasion):** The most metabolically active phase. The parasite digests host hemoglobin and rapidly synthesizes nucleic acids and proteins.
* **Expected DEGs:** You will observe strong upregulation in the Trophozoite stage of:
  - **Invasion factors:** Merozoite Surface Proteins (*msp1*, *msp3*, *msp4*), Apical Membrane Antigen 1 (*ama1*), Rhoptry-associated proteins.
  - **Hemoglobin catabolism:** Plasmepsins and falcipains.
  - **Replication/Transcription machinery:** Ribosomal proteins, helicases, and polymerases.

---

## ⚙️ Customization (config.yaml)
You can configure:
- **PlasmoDB release version:** change `plasmodb_release` in `config/config.yaml` to fetch older/newer references.
- **DESeq2 significance levels:** adjust `padj_cutoff` (default: `0.05`) and `lfc_cutoff` (default: `1.0`).
- **SRA sample runs:** edit `config/samples.tsv` to add or remove runs or change biological conditions.
