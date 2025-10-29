import os
import subprocess
import zipfile
import shutil
import glob
import sys

def calculate_genome_number(assembly_list):
    with open(assembly_list, "r") as f:
        genome_num = sum(1 for line in f if line.strip())
    return genome_num

def download_genomes(main_path,assembly_list):
    """
    Download genome assemblies from NCBI using accession numbers listed in a text file.

    Args:
        assembly_list (str): Path to the text file containing one accession per line.
        assembly_dir (str): Directory to store the downloaded genome FASTA files.
    """

    assembly_dir = f"{main_path}/genome_assemblies"

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
    return assembly_dir, assembly_list

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


def modify_name_path(ref_assembly):
    # Find all .fna files in the directory
    fna_files = glob.glob(os.path.join(ref_assembly, "*.fna"))

    if not fna_files:
        print(f"No .fna files found in {ref_assembly}")
        return

    for fna_file in fna_files:
        process_fna_file(fna_file)

def download_reference_genome(reference_genome: str, assembly_dir: str):
    """
    Download genome sequence (.fna) and annotation (.gff) from NCBI using the accession ID.

    Args:
        reference_genome (str): Genome accession ID (e.g. "GCF_000002655.1").
        assembly_dir (str): Directory to save downloaded and extracted files.

    Returns:
        tuple: Paths to the downloaded FASTA (.fna) and GFF (.gff) files.
    """
    output_dir = f"{assembly_dir}/reference_genome"

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    print(f"📂 Output directory: {output_dir}")

    # Check if 'datasets' command is available
    if shutil.which("datasets") is None:
        raise EnvironmentError(
            "❌ 'datasets' command not found. Please install it with:\n"
            "   conda install -c bioconda ncbi-datasets-cli"
        )

    # Define temporary ZIP file path
    zip_path = os.path.join(output_dir, f"{reference_genome}.zip")

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

        # Extract sequence (.fna) and annotation (.gff) to output_dir
        fasta_path, gff_path = None, None
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            for member in zip_ref.namelist():
                # Extract genomic sequence
                if member.endswith("genomic.fna"):
                    fasta_path = os.path.join(output_dir, f"{reference_genome}_genomic.fna")
                    with zip_ref.open(member) as src, open(fasta_path, "wb") as dst:
                        shutil.copyfileobj(src, dst)
                # Extract annotation (GFF)
                elif member.endswith("genomic.gff"):
                    gff_path = os.path.join(output_dir, f"{reference_genome}_genomic.gff")
                    with zip_ref.open(member) as src, open(gff_path, "wb") as dst:
                        shutil.copyfileobj(src, dst)

        if fasta_path:
            print(f"✅ FASTA saved: {fasta_path}")
        else:
            print("⚠️ No FASTA file found in the NCBI dataset.")

        if gff_path:
            print(f"✅ GFF saved:   {gff_path}")
        else:
            print("⚠️ No GFF file found in the NCBI dataset.")

    except subprocess.CalledProcessError as e:
        print(f"❌ Download failed for {reference_genome}: {e}")

    finally:
        # Clean up temp ZIP
        if os.path.exists(zip_path):
            os.remove(zip_path)

    return fasta_path, gff_path

def extract_annotations(input_gff: str, type_annotation: str, key_words: list = None):
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
    output_gff = f"{input_gff[:-4]}_{type_annotation}.gff"
    count = 0
    with open(input_gff, "r") as infile, open(output_gff, "w") as outfile:
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

    print(f"✅ Extracted {count} {type_annotation} entries to {output_gff}")
    return output_gff



def align_assemblies_to_reference(reference_fna: str, assembly_dir: str, main_path: str, species: str):
    """
    Align multiple genome assemblies (.fna) to a reference genome using minimap2 and samtools.

    Args:
        reference_fna (str): Path to the reference genome FASTA file.
        assembly_dir (str): Directory containing assembly .fna files to align.
        output_dir (str): Directory to store BAM and index files.
        output_prefix (str): Prefix for the output BAM file name.
    """
    # Create output directory if needed
    output_dir = f"{main_path}/alignment"
    os.makedirs(output_dir, exist_ok=True)

    # Find all .fna files
    assembly_files = sorted(glob.glob(os.path.join(assembly_dir, "*.fna")))
    if not assembly_files:
        print(f"❌ No .fna files found in {assembly_dir}")
        return

    print(f"🔹 Found {len(assembly_files)} assemblies for alignment.")
    print(f"   Reference: {reference_fna}")

    # Output BAM file
    output_bam = f"{output_dir}/alignment_{species}.sorted.bam"

    # Build minimap2 command pipeline
    cmd = (
        f"minimap2 -ax asm5 --secondary=no {reference_fna} {' '.join(assembly_files)} "
        f"| samtools view -bS - "
        f"| samtools sort -o {output_bam}"
    )

    print(f"\n🚀 Running alignment pipeline...\n{cmd}\n")

    # Run pipeline in shell
    subprocess.run(cmd, shell=True, check=True)

    # Index BAM file
    subprocess.run(["samtools", "index", output_bam], check=True)

    print(f"✅ Alignment complete.\n   BAM: {output_bam}\n   BAI: {output_bam}.bai")

    return output_bam

def prepare_anallyze_alignment(base_path, species, reference_genome, type_annotation,assembly_list, key_words):
    """"""

    # path to specific species
    main_path = f"{base_path}/{species}"

    # download assemblies
    assembly_dir, assembly_list = download_genomes(main_path,assembly_list)

    # download reference genome
    ref_assembly, ref_gff = download_reference_genome(reference_genome, assembly_dir)

    # modify the chromosome name of the assemblies
    modify_name_path(assembly_dir)

    # filter the annotation with type_annotation
    gff_filtered = extract_annotations(ref_gff, type_annotation, key_words)

    # alignment
    bam_path = align_assemblies_to_reference(ref_assembly, assembly_dir, main_path, species)

    return assembly_dir, ref_assembly, ref_gff, gff_filtered, bam_path


def run_bedtools_depth(gff_file, main_path, species, reference_genome):
    """
    Run bedtools coverage with -mean option to calculate average read depth over GFF regions.

    Args:
        gff_file (str): Path to the GFF3 annotation file.
        bam_file (str): Path to the BAM alignment file.
        output_file (str): Path to save the coverage results.
    """
    output_file = f"{main_path}/depth_calculation/{species}_{reference_genome}_meandepth.txt"
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

    print(f"✅ Coverage result saved to {output_file}")
    return bam_file

if __name__ == "__main__":
    base_path = "/lustre/BIF/nobackup/leng010/test"
    reference_genome = "GCF_000184455.3"  # genome annotation should be GCF version
    species = "aspergillus_oryzae"
    type_annotation = "gene"  # type of annotation used in depth calculation, the third column
    key_words = ["protein_coding"]
    assembly_list = f"{base_path}/genome_accessions/{species}_test.txt"
    main_path = f"{base_path}/{species}"

    assembly_dir, ref_assembly, ref_gff, gff_filtered, bam_path = prepare_anallyze_alignment(
        base_path,
        species,
        reference_genome,
        type_annotation,
        assembly_list,
        key_words)

    depth_path = run_bedtools_depth(gff_filtered, main_path, species, reference_genome)


