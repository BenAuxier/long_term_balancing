import os
import subprocess
import zipfile
import shutil
import glob

def download_genomes(assembly_list: str, assembly_dir: str):
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

def download_reference_genome(reference_genome: str, assembly_dir: str):
    """
    Download genome sequence (.fna) and annotation (.gff) from NCBI using the accession ID.

    Args:
        reference_genome (str): Genome accession ID (e.g. "GCF_000002655.1").
        assembly_dir (str): Directory to save downloaded and extracted files.

    Returns:
        tuple: Paths to the downloaded FASTA (.fna) and GFF (.gff) files.
    """
    reference_path = f"{assembly_dir}/reference_genome"

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

def extract_mrna_annotations(input_gff: str, type_annotation: str):
    """
    Extract lines with type_annotation (i.e., 'mRNA') in the 3rd column from a GFF file.

    Args:
        input_gff (str): Path to the input GFF file (e.g., "GCF_000002655.1_ASM265v1_genomic.gff")
        output_gff (str): Path to the output file (e.g., "GCF_000002655.1_ASM265v1_type_annotation.gff")
    """
    output_gff = f"{gff_path[:-4]}_{type_annotation}.gff"
    count = 0
    with open(input_gff, "r") as infile, open(output_gff, "w") as outfile:
        for line in infile:
            # Skip comment or empty lines
            if line.startswith("#") or not line.strip():
                continue

            # Split columns by tab and check 3rd column
            columns = line.strip().split("\t")
            if len(columns) >= 3 and columns[2] == type_annotation:
                outfile.write(line)
                count += 1

    print(f"✅ Extracted {count} {type_annotation} entries to {output_gff}")
    return output_gff

def align_assemblies_to_reference(reference_fna: str, assembly_dir: str, output_dir: str, output_prefix: str):
    """
    Align multiple genome assemblies (.fna) to a reference genome using minimap2 and samtools.

    Args:
        reference_fna (str): Path to the reference genome FASTA file.
        assembly_dir (str): Directory containing assembly .fna files to align.
        output_dir (str): Directory to store BAM and index files.
        output_prefix (str): Prefix for the output BAM file name.
    """

    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)

    # Find all .fna files
    assembly_files = sorted(glob.glob(os.path.join(assembly_dir, "*.fna")))
    if not assembly_files:
        print(f"❌ No .fna files found in {assembly_dir}")
        return

    print(f"🔹 Found {len(assembly_files)} assemblies for alignment.")
    print(f"   Reference: {reference_fna}")

    # Output BAM file
    output_bam = os.path.join(output_dir, f"{output_prefix}.sorted.bam")

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


if __name__ == "__main__":

    # download assemblies
    assembly_dir = "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/genome_assemblies"
    assembly_list = os.path.join(assembly_dir, "genome_accessions.txt")
    #download_genomes(assembly_list, assembly_dir)

    # download reference genome
    # Example usage
    reference_genome = "GCF_000002655.1"
    ref_assembly, ref_gff = download_reference_genome(reference_genome, assembly_dir)

    # filter the annotation with type_annotation
    type_annotation = "mRNA"
    gff_filtered = extract_mrna_annotations(ref_gff,type_annotation)

    # alignment
    output_dir = "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/alignment"
    output_prefix = f"Asper_{reference_genome}"

    bam_path = align_assemblies_to_reference(ref_assembly, assembly_dir, output_dir, output_prefix)

