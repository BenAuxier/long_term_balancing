#penicillium_chrysogenum.py
"""
Extract the allele of each gene in multiple genomes

"""
from prepare_alignment import prepare_anallyze_alignment
from prepare_alignment import calculate_genome_number
from depth_calculation import calculate_depth_all
from merge_region import process_results
from merge_region import process_data
from load_reference import load_annotation
from load_reference import count_gff_features
from analyze_position import analyze_all_candidate_position
from make_outputs import extract_outputs
from visualization_clinker import run_clinker_batch
from doublecheck_alignment import annotate_file_path
import os

# information, settings
reference_genome = "GCF_013085055.1" # genome annotation should be GCF version (RefSeq)
species = "fusarium_oxysporum" # the species
augustus_species = "fusarium_graminearum" # the reference species used in AUGUSTUS

type_annotation = "mRNA" # type of annotation used in depth calculation, the third column
ID_label = "transcript_id" # this is the key that the gene/mRNA id follows

key_words = None # the keywords that have to be included in the annotation

#########################################################
# file paths, including all input files
base_path = "/lustre/BIF/nobackup/leng010/test"

#assembly_list, this file need to create manually
assembly_list = f"{base_path}/genome_accessions/{species}.txt"
##########################################################################
# path to specific species
main_path = f"{base_path}/{species}"
os.makedirs(main_path, exist_ok=True)

assembly_dir, ref_assembly, ref_gff, gff_filtered, bam_path = prepare_anallyze_alignment(base_path, species, reference_genome, type_annotation,assembly_list, key_words)

#########################################################################
""""""
assembly_dir= f"{main_path}/genome_assemblies"
ref_assembly= f"{assembly_dir}/reference_genome/{reference_genome}_genomic.fna"
ref_gff= f"{assembly_dir}/reference_genome/{reference_genome}_genomic.gff"
gff_filtered= f"{assembly_dir}/reference_genome/{reference_genome}_genomic_{type_annotation}.gff"
bam_path = f"{main_path}/alignment/alignment_{species}.sorted.bam"

##########################################################################################
# verify some basic details
# check reference annotation .gff file
feature_counts = count_gff_features(ref_gff)

# calculate number of assembly used in alignment
genome_num = calculate_genome_number(assembly_list)

# settings
up_num = 5
down_num = 5
assembly_num = 7
lower_limit = genome_num * 0.2
upper_limit = genome_num * 0.8
minimal_alignment = genome_num * 0.3
extend = 5000

##########################################################################
# analyze the depth of the genomic regions of
depth_path = calculate_depth_all(bam_path, main_path, gff_filtered)
#depth_path = f"{main_path}/depth_calculation/mean_depth.txt"

# load annotation data from gff annotation
annotation_sorted, annotation_sorted_dict = load_annotation(gff_filtered, ID_label, type_annotation)

# processes the input candidate mRNAs
# Input and output file paths
candidate_data = process_data(depth_path, ID_label)
#print(candidate_data)
candidate_merge = process_results(depth_path,lower_limit, upper_limit,annotation_sorted_dict, ID_label)

# test the main code
#candidate_merge = dict(list(candidate_merge.items())[0:5])

candidate_data_summary = analyze_all_candidate_position(candidate_merge, annotation_sorted, candidate_data,
                        bam_path, assembly_list, up_num, down_num, lower_limit, minimal_alignment,type_annotation)

# extract sequences
results_path,sequence_path = extract_outputs(candidate_data_summary, reference_genome, ref_gff, main_path, extend, ref_assembly,
                          assembly_dir, assembly_num, candidate_data, augustus_species)

#run clinker
clinker_output_dir = run_clinker_batch(sequence_path, results_path)

# align sequence onto reference sequence to doublecheck and debug
output_main_path = f"{results_path}/sequence_alignments"
annotate_file_path(sequence_path,output_main_path)

print("finished")