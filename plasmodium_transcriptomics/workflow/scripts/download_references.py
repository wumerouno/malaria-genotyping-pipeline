import os
import sys
import urllib.request
import argparse

def download_file(url, dest):
    print(f"Downloading {url} to {dest}...")
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    try:
        urllib.request.urlretrieve(url, dest)
        print(f"Successfully downloaded to {dest}")
    except Exception as e:
        print(f"Error downloading {url}: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Download P. falciparum references from PlasmoDB.")
    parser.add_argument("--release", type=int, default=65, help="PlasmoDB release version (default: 65)")
    parser.add_argument("--genome-out", required=True, help="Path to save genome FASTA")
    parser.add_argument("--transcriptome-out", required=True, help="Path to save transcriptome FASTA")
    parser.add_argument("--gff-out", required=True, help="Path to save GFF3 annotation")
    
    args = parser.parse_args()
    
    base_url = f"https://plasmodb.org/common/downloads/release-{args.release}/Pfalciparum3D7"
    
    genome_url = f"{base_url}/fasta/data/PlasmoDB-{args.release}_Pfalciparum3D7_Genome.fasta"
    transcriptome_url = f"{base_url}/fasta/data/PlasmoDB-{args.release}_Pfalciparum3D7_AnnotatedTranscripts.fasta"
    gff_url = f"{base_url}/gff/data/PlasmoDB-{args.release}_Pfalciparum3D7.gff"
    
    download_file(genome_url, args.genome-out)
    download_file(transcriptome_url, args.transcriptome_out)
    download_file(gff_url, args.gff_out)

if __name__ == "__main__":
    # If run within Snakemake, we can parse snakemake.output
    if 'snakemake' in globals():
        release = snakemake.config.get("plasmodb_release", 65)
        # In Snakemake, output is a named list
        genome_out = snakemake.output.genome
        transcriptome_out = snakemake.output.transcriptome
        gff_out = snakemake.output.gff
        
        base_url = f"https://plasmodb.org/common/downloads/release-{release}/Pfalciparum3D7"
        genome_url = f"{base_url}/fasta/data/PlasmoDB-{release}_Pfalciparum3D7_Genome.fasta"
        transcriptome_url = f"{base_url}/fasta/data/PlasmoDB-{release}_Pfalciparum3D7_AnnotatedTranscripts.fasta"
        gff_url = f"{base_url}/gff/data/PlasmoDB-{release}_Pfalciparum3D7.gff"
        
        download_file(genome_url, genome_out)
        download_file(transcriptome_url, transcriptome_out)
        download_file(gff_url, gff_out)
    else:
        # Command line parsing (args.genome-out contains a dash, let's fix the argument parser destination)
        parser = argparse.ArgumentParser(description="Download P. falciparum references from PlasmoDB.")
        parser.add_argument("--release", type=int, default=65, help="PlasmoDB release version (default: 65)")
        parser.add_argument("--genome-out", dest="genome_out", required=True, help="Path to save genome FASTA")
        parser.add_argument("--transcriptome-out", dest="transcriptome_out", required=True, help="Path to save transcriptome FASTA")
        parser.add_argument("--gff-out", dest="gff_out", required=True, help="Path to save GFF3 annotation")
        
        args = parser.parse_args()
        
        base_url = f"https://plasmodb.org/common/downloads/release-{args.release}/Pfalciparum3D7"
        
        genome_url = f"{base_url}/fasta/data/PlasmoDB-{args.release}_Pfalciparum3D7_Genome.fasta"
        transcriptome_url = f"{base_url}/fasta/data/PlasmoDB-{args.release}_Pfalciparum3D7_AnnotatedTranscripts.fasta"
        gff_url = f"{base_url}/gff/data/PlasmoDB-{args.release}_Pfalciparum3D7.gff"
        
        download_file(genome_url, args.genome_out)
        download_file(transcriptome_url, args.transcriptome_out)
        download_file(gff_url, args.gff_out)
