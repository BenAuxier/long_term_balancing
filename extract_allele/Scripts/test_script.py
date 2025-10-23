import os
import sys
import glob

def process_fna_file(filepath):
    """
    Process a single .fna file:
    1. Extract prefix from filename by removing the suffix 'genomic.fna'
    2. Replace each FASTA header (lines starting with '>') with '>prefix_originalHeader'
    3. Save the modified content back to the same file
    """
    filename = os.path.basename(filepath)
    if filename.endswith("genomic.fna"):
        prefix = filename[:-len("genomic.fna")]
    else:
        prefix = os.path.splitext(filename)[0]  # fallback: remove extension

    # Read original file
    with open(filepath, "r") as f:
        lines = f.readlines()

    # Modify headers
    new_lines = []
    for line in lines:
        if line.startswith(">"):
            # Insert prefix after ">"
            line = ">" + prefix + line[1:]
        new_lines.append(line)

    # Write back to the same file
    with open(filepath, "w") as f:
        f.writelines(new_lines)

    print(f"Processed: {filename} with prefix '{prefix}'")


def modify_name_path(assembly_dir):
    # Find all .fna files in the directory
    fna_files = glob.glob(os.path.join(assembly_dir, "*.fna"))

    if not fna_files:
        print(f"No .fna files found in {assembly_dir}")
        return

    for fna_file in fna_files:
        process_fna_file(fna_file)

assembly_dir= f"/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/genome_assemblies"
modify_name_path(assembly_dir)