import os
import subprocess

def run_clinker_batch(sequence_path, results_path):
    """
    Batch run Clinker on all subdirectories containing multiple valid GFF3 files.

    Args:
        sequence_path (str): Path to the parent directory containing multiple subdirectories with .gff3 files.

    The function performs the following steps:
        1. Traverse all first-level subdirectories under sequence_path.
        2. For each subdirectory:
           - Identify valid GFF3 files (non-empty and with content beyond headers).
           - If >=2 valid files exist, run Clinker and save the HTML plot to the output folder.
    """

    output_dir = os.path.join(results_path, "clinker_results")
    os.makedirs(output_dir, exist_ok=True)

    # Traverse each subdirectory under sequence_path
    for gff_dir in os.listdir(sequence_path):
        gff_dir_path = os.path.join(sequence_path, gff_dir)
        if not os.path.isdir(gff_dir_path):
            continue

        print(f"\nProcessing directory: {gff_dir_path}")
        gff_files = []

        # Find all .gff3 files (non-recursive)
        for file in os.listdir(gff_dir_path):
            if file.endswith(".gff3"):
                file_path = os.path.join(gff_dir_path, file)

                # Count non-comment lines
                with open(file_path, "r") as f:
                    non_comment_lines = sum(1 for line in f if not line.startswith("#"))

                if non_comment_lines > 0:
                    print(f"Valid GFF3 file: {file_path} ({non_comment_lines} non-comment lines)")
                    gff_files.append(file_path)
                else:
                    print(f"Skipping empty or header-only file: {file_path}")

        # Must have at least 2 valid GFF3 files
        if len(gff_files) < 2:
            print(f"❌ Error: Less than two valid GFF3 files found in {gff_dir_path}")
            continue

        # Run clinker
        gene_name = os.path.basename(gff_dir_path)
        output_html = os.path.join(output_dir, f"{gene_name}_plot.html")

        print(f"🧬 Running Clinker on {len(gff_files)} valid GFF3 files...")
        try:
            subprocess.run(["clinker", "--force", *gff_files, "-p", output_html], check=True)
            print(f"✅ Clinker analysis finished for {gff_dir_path}. Results saved to {output_html}")
        except subprocess.CalledProcessError as e:
            print(f"⚠️ Error running Clinker on {gff_dir_path}: {e}")

    return output_dir

if __name__ == "__main__":
    # Example usage
    results_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/results"
