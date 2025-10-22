"""
Extract the allele of each gene in multiple genomes

"""
from calculate_coverage import run_bedtools_coverage
from mix_region import process_results
from mix_region import process_data
from load_reference import load_annotation
from analyze_position import analyze_all_candidate_position
from make_outputs import extract_outputs
from visualization_clinker import run_clinker_batch




# settings
up_num = 5
down_num = 5
assembly_num = 7
lower_limit = 10
upper_limit = 40
extend = 5000
reference_genome = "GCF_000002655.1"
species = "aspergillus_fumigatus"

##########################################################################
# file paths, including all input files
# The modified csv file that contains the mRNA with depth of interests.
main_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/check_coverage"
depth_path = f"{main_path}/depth_calculation/{species}_{reference_genome}_meandepth.txt"

# Reference genome annotation of the BAM file
reference_path = f"{main_path}/reference_genome"
gff_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/reference/GCF_000002655.1_ASM265v1_genomic.gff"

# Reference genome assembly
ref_assembly = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/reference/GCF_000002655.1_ASM265v1_genomic.fna"

# Bam file, used multiple genome assemblies align to the reference genome
bam_path = ("/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/multi_align_2/40_assembly/normal_align/"
            "all51_to_GCF_000002655.1_asm5.sorted.bam")
# The txt file that includes the genome assemblies used in alignment.
assembly_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/multi_align_2/40_assembly/accessions.txt"
# path to all genome assemblies
assembly_dir = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/multi_align_2/40_assembly/genome_assemblies"
# output base path of the fasta files
output_path = f"{main_path}/test4_extract_seq/extract_allele"

#########################################################################
# download genome assemblies

# download reference genome and related annotation

# edit the names of the contigs

# alignment










##########################################################################
# analyze the depth of the genomic regions of
run_bedtools_coverage(gff_path, bam_path, depth_path)

# load annotation data from gff annotation
annotation_sorted, annotation_sorted_dict = load_annotation(gff_path,"mRNA")

# processes the input candidate mRNAs
# Input and output file paths
candidate_data = process_data(depth_path)
candidate_merge = process_results(depth_path,lower_limit, upper_limit,annotation_sorted_dict)


# test the main code
candidate_data_test = dict(list(candidate_merge.items())[0:])
# print(candidate_data_test)
candidate_data_summary = analyze_all_candidate_position(candidate_data_test, annotation_sorted, candidate_data,
                                                        bam_path, assembly_path, up_num, down_num, lower_limit)

# extract sequences
outputs = extract_outputs(candidate_data_summary, reference_genome, gff_path, output_path, extend, ref_assembly,
                          assembly_dir,assembly_num,candidate_data,species)

#run clinker
run_clinker_batch(output_path)

print("finished")