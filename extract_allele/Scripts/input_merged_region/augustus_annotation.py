#!/usr/bin/env python3
"""
Batch gene prediction with AUGUSTUS (Python version)

Species: Aspergillus fumigatus
Input: all .fa files in the given directory (recursively)
Output: GFF3 files saved in the same directory as the input .fa
"""

import os
import subprocess

def run_augustus_on_fasta(fa_path,augustus_species, gff3_status = "off", suffix = ""):
    """
    Run AUGUSTUS on a single FASTA file and save the output in GFF3 format.
    """
    base = os.path.splitext(os.path.basename(fa_path))[0]
    dir_path = os.path.dirname(fa_path)

    if gff3_status == "on":
        output_path = os.path.join(dir_path, f"{base}{suffix}.gff3")
    else:
        output_path = os.path.join(dir_path, f"{base}{suffix}.gff")

    print(f"Running AUGUSTUS on {fa_path} ...")

    # Build command
    cmd = [
        "augustus",
        f"--species={augustus_species}",
        f"--gff3={gff3_status}",
        fa_path
    ]

    # Run command and write output
    with open(output_path, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, check=True)
    return output_path

def annotate_file_path(input_dir,augustus_species, gff3_status = "off"):
    """
    Recursively find all .fa files under INPUT_DIR and run AUGUSTUS for each.
    """

    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".fa"):
                fa_path = os.path.join(root, file)
                try:
                    run_augustus_on_fasta(fa_path,augustus_species, gff3_status)
                except subprocess.CalledProcessError as e:
                    print(f"Error running AUGUSTUS on {fa_path}: {e}")
    return True

if __name__ == "__main__":
    # === Configuration ===
    # input_dir = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/check_coverage/test4_extract_seq/extract_allele/XM_741963.1"
    augustus_species = "aspergillus_oryzae"
    input_dir = "/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/genome_assemblies/reference_genome/GCF_000184455.2_genomic.fna"
    output_path = run_augustus_on_fasta(input_dir, augustus_species, gff3_status="off", suffix="_AUGUSTUS")
