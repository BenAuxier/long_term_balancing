import subprocess
import os

from visualization_clinker import run_clinker_batch

# modify annotation and rerun clinker to find where is the problem

def run_augustus_on_fasta(fa_path: str,augustus_species):
    """
    Run AUGUSTUS on a single FASTA file and save the output in GFF3 format.
    """
    base = os.path.splitext(os.path.basename(fa_path))[0]
    dir_path = os.path.dirname(fa_path)
    output_gff3 = os.path.join(dir_path, f"{base}.gff3")

    print(f"Running AUGUSTUS on {fa_path} ...")

    # Build command
    cmd = [
        "augustus",
        f"--species={augustus_species}",
        "--gff3=on",
        fa_path
    ]

    # Run command and write output
    with open(output_gff3, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, check=True)

def annotate_file_path(input_dir,augustus_species,keyword):
    """
    Recursively find all .fa files under INPUT_DIR and run AUGUSTUS for each.
    """

    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(keyword):
                fa_path = os.path.join(root, file)
                try:

                    run_augustus_on_fasta(fa_path,augustus_species)
                except subprocess.CalledProcessError as e:
                    print(f"Error running AUGUSTUS on {fa_path}: {e}")
    return True


genomic_region = "XM_003708869.1"

species = "magnaporthe_grisea"
oringin_path = f"/lustre/BIF/nobackup/leng010/test/magnaporthe_grisea/extract_sequences/{genomic_region}"
path = "/lustre/BIF/nobackup/leng010/test/magnaporthe_grisea/test_clinker"
subprocess.run(["cp", "-r", oringin_path, path], check=True)
test_path = f"/lustre/BIF/nobackup/leng010/test/magnaporthe_grisea/test_clinker/{genomic_region}"

# copy a file (such as reference) to check annotation in clinker plot
def copy_file(test_path):
    for root, _, files in os.walk(test_path):
        for file in files:
            if file.startswith(f"{genomic_region}_reference_genome"):
                file_path = f"{test_path}/{file}"
                file_copy_path = f"{test_path}/copy_{file}"
                print(["cp", file_path, file_copy_path])
                subprocess.run(["cp", file_path, file_copy_path], check=True)

copy_file(test_path)

keyword = ".fa"
#annotate_file_path(test_path,species, keyword)

# run clinker
clinker_output_dir = run_clinker_batch(path, test_path)




