import os
import sys

def generate_decoy_files(genome_path, transcriptome_path, decoy_out, gentrome_out):
    print(f"Parsing genome {genome_path} for decoy headers...")
    decoys = []
    
    # 1. Extract chromosome names from genome FASTA
    try:
        with open(genome_path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    # Extract identifier (e.g. ">Pf3D7_01_v3" -> "Pf3D7_01_v3")
                    header = line.strip().split()[0][1:]
                    decoys.append(header)
    except Exception as e:
        print(f"Error reading genome file: {e}")
        sys.exit(1)
        
    print(f"Found {len(decoys)} chromosomes/scaffolds to use as decoys.")
    
    # Write decoys to file
    try:
        with open(decoy_out, "w") as f:
            for decoy in decoys:
                f.write(decoy + "\n")
        print(f"Decoy list written to {decoy_out}")
    except Exception as e:
        print(f"Error writing decoy file: {e}")
        sys.exit(1)
        
    # 2. Concatenate transcriptome and genome into gentrome
    print(f"Creating gentrome file at {gentrome_out}...")
    try:
        with open(gentrome_out, "w") as out_f:
            # First write transcriptome
            with open(transcriptome_path, "r") as tx_f:
                for line in tx_f:
                    out_f.write(line)
            # Then append genome
            with open(genome_path, "r") as g_f:
                for line in g_f:
                    out_f.write(line)
        print(f"Gentrome successfully created at {gentrome_out}")
    except Exception as e:
        print(f"Error creating gentrome file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if 'snakemake' in globals():
        # Running inside Snakemake
        generate_decoy_files(
            genome_path=snakemake.input.genome,
            transcriptome_path=snakemake.input.transcriptome,
            decoy_out=snakemake.output.decoys,
            gentrome_out=snakemake.output.gentrome
        )
    else:
        # Running from CLI
        import argparse
        parser = argparse.ArgumentParser(description="Create Salmon decoy files.")
        parser.add_argument("--genome", required=True, help="Path to genome FASTA")
        parser.add_argument("--transcriptome", required=True, help="Path to transcriptome FASTA")
        parser.add_argument("--decoy-out", required=True, help="Path to save decoy headers list")
        parser.add_argument("--gentrome-out", required=True, help="Path to save gentrome FASTA")
        
        args = parser.parse_args()
        generate_decoy_files(args.genome, args.transcriptome, args.decoy_out, args.gentrome_out)
