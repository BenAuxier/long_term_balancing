import os
from prepare_alignment import align_assemblies_to_reference


base_path = "/lustre/BIF/nobackup/leng010/test"
reference_genome = "GCF_000002495.2"
species = "magnaporthe_grisea"
main_path = f"{base_path}/{species}"

ref_assembly= f"{main_path}/genome_assemblies/reference_genome/{reference_genome}_genomic.fna"
assembly_dir= f"{main_path}/genome_assemblies"
output_dir = f"{main_path}/alignment_test"
os.makedirs(output_dir, exist_ok=True)  # Create output directory if needed
bam_path = f"{output_dir}/alignment_{species}.sorted.bam"

bam_path = align_assemblies_to_reference(ref_assembly, assembly_dir, output_dir, bam_path)

















