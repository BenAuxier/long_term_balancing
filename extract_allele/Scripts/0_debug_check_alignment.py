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
        return ref_sequence,ref_annotation
    else:
        return None

def annotate_file_path(input_dir, output_main_path):
    """
    Recursively find all .fa files under INPUT_DIR and run alignment for each.
    """
    for root, dirs, files in os.walk(input_dir):
        for dir in dirs:
            genomic_region = dir
            sequence_path = f"{input_dir}/{genomic_region}"

            print("==========================================",
                  f"Start alignment for genomic region {genomic_region}",
                  f" based files in {sequence_path}.",
                  "==========================================")

            # find reference annotation
            ref_sequence, ref_annotation = find_reference_sequence(sequence_path)

            # generate output path
            output_path = f"{output_main_path}/{genomic_region}"
            os.makedirs(output_path, exist_ok=True)
            bam_file = f"{output_path}/{genomic_region}_alignment.sorted.bam"

            # make alignment
            align_assemblies_to_reference(ref_sequence, sequence_path, bam_file)

            # copy the reference information to the output_path
            subprocess.run(["cp", ref_sequence, output_path], check=True)
            subprocess.run(["cp", ref_annotation, output_path], check=True)

            print(f"Finished alignment for {genomic_region}! Files are saved to {output_path}")
    return True


if __name__ == "__main__":
    base_path = "/lustre/BIF/nobackup/leng010/test"
    species = "magnaporthe_grisea"
    main_path = f"{base_path}/{species}"
    input_dir = f"{main_path}/extract_sequences"
    output_main_path = f"{main_path}/results/sequence_alignments"
    annotate_file_path(input_dir,output_main_path)
