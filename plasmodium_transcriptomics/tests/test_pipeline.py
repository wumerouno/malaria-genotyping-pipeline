import os
import unittest
import yaml
import pandas as pd

class TestBioinformaticsPipeline(unittest.TestCase):
    def setUp(self):
        # Base directory of the transcriptomics project
        self.base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.config_path = os.path.join(self.base_dir, "config", "config.yaml")
        self.samples_path = os.path.join(self.base_dir, "config", "samples.tsv")

    def test_config_file_exists(self):
        """Check if config.yaml exists and is valid YAML."""
        self.assertTrue(os.path.exists(self.config_path), f"config.yaml not found at {self.config_path}")
        with open(self.config_path, "r") as f:
            try:
                config = yaml.safe_load(f)
                self.assertIsNotNone(config, "config.yaml is empty")
                # Check for critical keys
                required_keys = ["plasmodb_release", "samples_tsv", "reference_dir", "qc_dir", "quant_dir", "de_dir", "plot_dir"]
                for key in required_keys:
                    self.assertIn(key, config, f"Missing required key '{key}' in config.yaml")
            except yaml.YAMLError as exc:
                self.fail(f"Error parsing config.yaml: {exc}")

    def test_samples_file_exists_and_format(self):
        """Check if samples.tsv exists and has correct columns and values."""
        self.assertTrue(os.path.exists(self.samples_path), f"samples.tsv not found at {self.samples_path}")
        try:
            df = pd.read_csv(self.samples_path, sep="\t")
            # Check for required columns
            required_cols = ["sample_id", "run_id", "stage", "replicate"]
            for col in required_cols:
                self.assertIn(col, df.columns, f"Missing column '{col}' in samples.tsv")
            
            # Check SRA ID format (starts with SRR and followed by digits)
            for sra_id in df["run_id"]:
                self.assertTrue(sra_id.startswith("SRR"), f"Invalid SRA Run ID '{sra_id}', must start with 'SRR'")
                self.assertTrue(sra_id[3:].isdigit(), f"Invalid SRA Run ID suffix '{sra_id}'")
                
            # Check that stage has valid values
            valid_stages = ["Ring", "Trophozoite", "Schizont", "Gametocyte"]
            for stage in df["stage"]:
                self.assertIn(stage, valid_stages, f"Stage '{stage}' is invalid. Must be one of {valid_stages}")
        except Exception as e:
            self.fail(f"Error reading or parsing samples.tsv: {e}")

    def test_script_files_exist(self):
        """Verify that all pipeline scripts are present in the repository."""
        scripts = [
            os.path.join(self.base_dir, "workflow", "scripts", "download_references.py"),
            os.path.join(self.base_dir, "workflow", "scripts", "create_decoy_list.py"),
            os.path.join(self.base_dir, "src", "deseq2_analysis.R"),
            os.path.join(self.base_dir, "notebooks", "rna_seq_report.Rmd"),
            os.path.join(self.base_dir, "workflow", "Snakefile")
        ]
        for script_path in scripts:
            self.assertTrue(os.path.exists(script_path), f"Required file not found at: {script_path}")

if __name__ == "__main__":
    unittest.main()
