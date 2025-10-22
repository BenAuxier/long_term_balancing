import subprocess
import os


def run_bedtools_coverage(gff_file: str, bam_file: str, output_file: str):
    """
    Run bedtools coverage with -mean option to calculate average read depth over GFF regions.

    Args:
        gff_file (str): Path to the GFF3 annotation file.
        bam_file (str): Path to the BAM alignment file.
        output_file (str): Path to save the coverage results.
    """

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Build the bedtools command
    cmd = ["bedtools", "coverage", "-a", gff_file, "-b", bam_file, "-mean"]

    print(f"Running bedtools coverage...\nCommand: {' '.join(cmd)}")
    print(f"Output will be saved to: {output_file}")

    # Run the command and redirect output to file
    with open(output_file, "w") as out:
        subprocess.run(cmd, check=True, stdout=out)

    print(f"✅ Coverage result saved to {output_file}")


if __name__ == "__main__":
    # Example usage
    gff_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/reference/GCF_000002655.1_ASM265v1_mRNA.gff"
    bam_path = ("/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/multi_align_2/40_assembly/normal_align/"
                "all51_to_GCF_000002655.1_asm5.sorted.bam")

    output_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/check_coverage/all51_to_GCF_000002655.1_meandepth.primary.txt"

    run_bedtools_coverage(gff_path, bam_path, output_path)