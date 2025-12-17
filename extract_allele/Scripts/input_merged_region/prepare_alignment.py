import os
import subprocess
import zipfile
import shutil
import glob
import sys
from make_outputs import run_augustus_on_fasta

def calculate_genome_number(assembly_list):
    with open(assembly_list, "r") as f:
        genome_num = sum(1 for line in f if line.strip())
    return genome_num

def download_genomes(main_path,assembly_list, assembly_dir):
    """
    Download genome assemblies from NCBI using accession numbers listed in a text file.

    Args:
        assembly_list (str): Path to the text file containing one accession per line.
        assembly_dir (str): Directory to store the downloaded genome FASTA files.
    """
    # Create output directory if it does not exist
    os.makedirs(assembly_dir, exist_ok=True)
    print(f"Output directory: {assembly_dir}")

    # Check if 'datasets' command is available
    if shutil.which("datasets") is None:
        raise EnvironmentError(
            "❌ 'datasets' command not found. Please install it using:\n"
            "   conda install -c bioconda ncbi-datasets-cli"
        )

    # Read accession list
    with open(assembly_list, "r") as f:
        accessions = [line.strip() for line in f if line.strip()]

    if not accessions:
        print("No accession numbers found in the file.")
        return

    # Process each accession
    for acc in accessions:
        print(f"\n🔹 Downloading {acc} ...")
        zip_path = os.path.join(assembly_dir, f"{acc}.zip")

        try:
            # Run datasets download
            subprocess.run(
                ["datasets", "download", "genome", "accession", acc, "--filename", zip_path],
                check=True
            )

            # Extract only genomic .fna files into assembly_dir (flatten structure)
            with zipfile.ZipFile(zip_path, "r") as zip_ref:
                members = [m for m in zip_ref.namelist() if "genomic.fna" in m]
                if not members:
                    print(f"⚠️ No .fna files found in {acc}")
                else:
                    for member in members:
                        filename = os.path.basename(member)
                        target_path = os.path.join(assembly_dir, f"{acc}_{filename}")
                        with zip_ref.open(member) as source, open(target_path, "wb") as target:
                            shutil.copyfileobj(source, target)
                    print(f"✅ Extracted {len(members)} FASTA file(s) for {acc}")

        except subprocess.CalledProcessError:
            print(f"❌ Error downloading {acc}")

        finally:
            # Clean up the temporary zip
            if os.path.exists(zip_path):
                os.remove(zip_path)

    print("\n🎉 All downloads completed successfully.")
    return assembly_dir

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

def download_reference_genome(reference_genome, ref_path, ref_assembly, ref_gff):
    """
    Download genome sequence (.fna) and annotation (.gff) from NCBI using the accession ID.

    Args:
        reference_genome (str): Genome accession ID (e.g. "GCF_000002655.1").
        assembly_dir (str): Directory to save downloaded and extracted files.

    Returns:
        tuple: Paths to the downloaded FASTA (.fna) and GFF (.gff) files.
    """
    # Ensure output directory exists
    os.makedirs(ref_path, exist_ok=True)
    print(f"📂 Output directory: {ref_path}")

    # Check if 'datasets' command is available
    if shutil.which("datasets") is None:
        raise EnvironmentError(
            "❌ 'datasets' command not found. Please install it with:\n"
            "   conda install -c bioconda ncbi-datasets-cli"
        )

    # Define temporary ZIP file path
    zip_path = os.path.join(ref_path, f"{reference_genome}.zip")

    print(f"\n🔹 Downloading {reference_genome} from NCBI...")

    try:
        # Use NCBI datasets CLI to download genome + annotation
        subprocess.run(
            [
                "datasets", "download",
                "genome", "accession", reference_genome,
                "--filename", zip_path,
                "--include", "genome,gff3"
            ],
            check=True
        )

        # Extract sequence (.fna) and annotation (.gff) to ref_path
        ref_assembly, ref_gff = None, None
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            for member in zip_ref.namelist():
                # Extract genomic sequence
                if member.endswith("genomic.fna"):
                    ref_assembly = os.path.join(ref_path, f"{reference_genome}_genomic.fna")
                    with zip_ref.open(member) as src, open(ref_assembly, "wb") as dst:
                        shutil.copyfileobj(src, dst)
                # Extract annotation (GFF)
                elif member.endswith("genomic.gff"):
                    ref_gff = os.path.join(ref_path, f"{reference_genome}_genomic.gff")
                    with zip_ref.open(member) as src, open(ref_gff, "wb") as dst:
                        shutil.copyfileobj(src, dst)

        if ref_assembly:
            print(f"✅ FASTA saved: {ref_assembly}")
        else:
            print("⚠️ No FASTA file found in the NCBI dataset.")

        if ref_gff:
            print(f"✅ GFF saved:   {ref_gff}")
        else:
            print("⚠️ No GFF file found in the NCBI dataset.")

    except subprocess.CalledProcessError as e:
        print(f"❌ Download failed for {reference_genome}: {e}")

    finally:
        # Clean up temp ZIP
        if os.path.exists(zip_path):
            os.remove(zip_path)

    return ref_assembly, ref_gff

def extract_annotations(ref_gff, gff_filtered, type_annotation, key_words):
    """
    Extract lines with a given type_annotation (e.g., 'mRNA') in the 3rd column from a GFF file,
    and optionally filter lines that contain all key_words.

    Args:
        input_gff (str): Path to the input GFF file (e.g., "GCF_000002655.1_ASM265v1_genomic.gff")
        type_annotation (str): The annotation type to extract (e.g., "mRNA")
        key_words (list, optional): List of required keywords that must all appear in the line.
                                    If None or empty, no keyword filtering is applied.

    Returns:
        str: Path to the output GFF file
    """
    count = 0
    with open(ref_gff, "r") as infile, open(gff_filtered, "w") as outfile:
        for line in infile:
            # Skip comment or empty lines
            if line.startswith("#") or not line.strip():
                continue

            # Split columns by tab and check 3rd column
            columns = line.strip().split("\t")
            if len(columns) >= 3 and columns[2] == type_annotation:
                # If key_words are provided, then each keyword must appear in that line.
                if key_words:
                    if all(kw in line for kw in key_words):
                        outfile.write(line)
                        count += 1
                else:
                    # No keyword filtering, just enter
                    outfile.write(line)
                    count += 1

    print(f"✅ Extracted {count} {type_annotation} entries to {gff_filtered}")
    return gff_filtered


def align_assemblies_to_reference(ref_assembly, assembly_dir, bam_path, bam_file):
    """
    Align multiple genome assemblies (.fna) to a reference genome using minimap2 and samtools.

    Args:
        reference_fna (str): Path to the reference genome FASTA file.
        assembly_dir (str): Directory containing assembly .fna files to align.
        output_dir (str): Directory to store BAM and index files.
        output_prefix (str): Prefix for the output BAM file name.
    """
    # Create output directory if needed
    os.makedirs(bam_path, exist_ok=True)

    # Find all .fna files
    assembly_files = sorted(
        glob.glob(os.path.join(assembly_dir, "*.fna")) +
        glob.glob(os.path.join(assembly_dir, "*.fa"))
    )

    if not assembly_files:
        print(f"❌ No .fna or .fa files found in {assembly_dir}")
        return

    print(f"🔹 Found {len(assembly_files)} assemblies for alignment.")
    print(f"   Reference: {ref_assembly}")

    # Build minimap2 command pipeline
    cmd = (
        f"minimap2 -ax asm5 --secondary=no {ref_assembly} {' '.join(assembly_files)} "
        f"| samtools view -bS - "
        f"| samtools sort -o {bam_file}"
    )

    print(f"\n🚀 Running alignment pipeline...\n{cmd}\n")

    # Run pipeline in shell
    subprocess.run(cmd, shell=True, check=True)

    # Index BAM file
    subprocess.run(["samtools", "index", bam_file], check=True)

    print(f"✅ Alignment complete.\n   BAM: {bam_file}\n   BAI: {bam_file}.bai")

    return bam_file

def prepare_analyze_alignment(main_path, assembly_dir, ref_path, ref_assembly, ref_gff, bam_path, bam_file, reference_genome, assembly_list):
    """"""
    # download assemblies
    download_genomes(main_path,assembly_list, assembly_dir)
    print("assemblies downloaded")

    # modify the chromosome name of the assemblies
    modify_name_path(assembly_dir)

    # download reference genome
    ref_assembly, ref_gff = download_reference_genome(reference_genome, ref_path, ref_assembly, ref_gff)

    # alignment
    align_assemblies_to_reference(ref_assembly, assembly_dir, bam_path, bam_file)

    return bam_path

def prepare_augustus_reference(ref_assembly, augustus_species):
    # annotate reference genome with augustus, output gff
    ref_gff_augustus = run_augustus_on_fasta(ref_assembly, augustus_species, gff3_status="off", suffix="_AUGUSTUS")

    # annotate reference genome with augustus, output gff
    ref_gff3_augustus = run_augustus_on_fasta(ref_assembly, augustus_species, gff3_status="on", suffix="_AUGUSTUS")


def run_bedtools_depth(gff_file, main_path, species, reference_genome):
    """
    Run bedtools coverage with -mean option to calculate average read depth over GFF regions.

    Args:
        gff_file (str): Path to the GFF3 annotation file.
        bam_file (str): Path to the BAM alignment file.
        output_file (str): Path to save the coverage results.
    """
    output_file = f"{main_path}/depth_calculation/{species}_{reference_genome}_meandepth_test.txt"
    bam_file = f"{main_path}/alignment/alignment_{species}.sorted.bam"

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Build the bedtools command
    cmd = ["bedtools", "coverage", "-a", gff_file, "-b", bam_file, "-mean"]

    print(f"Running bedtools coverage...\nCommand: {' '.join(cmd)}")
    print(f"Output will be saved to: {output_file}")

    # Run the command and redirect output to file
    with open(output_file, "w") as out:
        subprocess.run(cmd, check=True, stdout=out)

    print(f"✅ Result saved to {output_file}")
    return output_file

if __name__ == "__main__":
    base_path = "/lustre/BIF/nobackup/leng010/test"
    reference_genome = "GCF_000184455.3"  # genome annotation should be GCF version
    species = "aspergillus_oryzae"
    type_annotation = "gene"  # type of annotation used in depth calculation, the third column
    key_words = ["protein_coding"]
    assembly_list = f"{base_path}/genome_accessions/{species}_test.txt"
    main_path = f"{base_path}/{species}"

    assembly_dir, ref_assembly, ref_gff, gff_filtered, bam_path = prepare_analyze_alignment(
        base_path,
        species,
        reference_genome,
        type_annotation,
        assembly_list,
        key_words)

    depth_path = run_bedtools_depth(gff_filtered, main_path, species, reference_genome)


