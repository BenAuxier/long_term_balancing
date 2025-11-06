import os
from prepare_alignment import align_assemblies_to_reference

import os
import glob
import subprocess


def find_reference_sequence(sequence_path: str) -> str:
    """
    Extract the genomic_region from the specified path and locate the corresponding reference genome file (.fa).
    Return the full path of the matched file; return None if not found.
    """
    # Extract genomic_region
    genomic_region = os.path.basename(sequence_path)

    # Constructing Matching Patterns
    sequence_pattern = f"{sequence_path}/{genomic_region}_reference_genome*.fa"
    annotation_pattern = f"{sequence_path}/{genomic_region}_reference_genome*.gff3"

    # Find matching files
    ref_sequence = glob.glob(sequence_pattern)[0]
    ref_annotation = glob.glob(annotation_pattern)[0]

    # return results
    if sequence_pattern and annotation_pattern:
        return genomic_region, ref_sequence,ref_annotation
    else:
        return None


base_path = "/lustre/BIF/nobackup/leng010/test"
species = "magnaporthe_grisea"
main_path = f"{base_path}/{species}"

sequence_path = "/lustre/BIF/nobackup/leng010/test/magnaporthe_grisea/extract_sequences/XM_003708866.1"

# find reference annotation
genomic_region, ref_sequence, ref_annotation = find_reference_sequence(sequence_path)

# generate output path
output_path = f"{main_path}/results/sequence_alignments/{genomic_region}"
os.makedirs(output_path, exist_ok=True)
bam_file = f"{output_path}/{genomic_region}_alignment.sorted.bam"

#make alignment
#align_assemblies_to_reference(ref_sequence, sequence_path, bam_file)
#subprocess.run(["cp", ref_sequence, output_path], check=True)
#subprocess.run(["cp", ref_annotation, output_path], check=True)


def annotate_file_path(input_dir):
    """
    Recursively find all .fa files under INPUT_DIR and run alignment for each.
    """
    for root, dirs, files in os.walk(input_dir):
        for dir in dirs:
            genomic_region = dir
            sequence_path = f"{input_dir}/{dir}"
            print(genomic_region,sequence_path)


input_dir = "/lustre/BIF/nobackup/leng010/test/magnaporthe_grisea/extract_sequences"
annotate_file_path(input_dir)