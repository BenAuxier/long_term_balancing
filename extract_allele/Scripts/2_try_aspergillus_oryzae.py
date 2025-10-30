"""
Extract the allele of each gene in multiple genomes

"""
from prepare_alignment import prepare_anallyze_alignment
from prepare_alignment import run_bedtools_depth
from prepare_alignment import calculate_genome_number
from merge_region import process_results
from merge_region import process_data
from load_reference import load_annotation
from load_reference import count_gff_features
from analyze_position import analyze_all_candidate_position
from make_outputs import extract_outputs
from visualization_clinker import run_clinker_batch
import os

# information
reference_genome = "GCF_000184455.2" # genome annotation should be GCF version (RefSeq)
species = "aspergillus_oryzae" # the species
augustus_species = "aspergillus_oryzae" # the reference species used in AUGUSTUS
type_annotation = "mRNA" # type of annotation used in depth calculation, the third column
ID_label = "transcript_id" # this is the key that the gene/mRNA id follows

key_words = None # the keywords that have to be included in the annotation

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
"""
assembly_dir = f"{main_path}/genome_assemblies" #path to all genome assemblies

# The txt file that includes the genome assemblies used in alignment.
# need to be created manually!
assembly_list = f"{main_path}/genome_accessions.txt 

ref_assembly # Reference genome assembly
ref_gff # Reference genome annotation
bam_path # Path to bam file

assembly_dir= f"/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/genome_assemblies"
assembly_list= f"/lustre/BIF/nobackup/leng010/test/genome_accessions/aspergillus_oryzae.txt"
ref_assembly= "/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/genome_assemblies/reference_genome/GCF_000184455.2_genomic.fna"
ref_gff= "/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/genome_assemblies/reference_genome/GCF_000184455.2_genomic.gff"
gff_filtered= "/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/genome_assemblies/reference_genome/GCF_000184455.2_genomic_mRNA.gff"
bam_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/alignment/alignment_aspergillus_oryzae.sorted.bam"
"""
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
depth_path = run_bedtools_depth(gff_filtered, main_path, species, reference_genome)
#depth_path = f"{main_path}/depth_calculation/{species}_{reference_genome}_meandepth.txt"

# load annotation data from gff annotation
annotation_sorted, annotation_sorted_dict = load_annotation(gff_filtered, ID_label, type_annotation)

# processes the input candidate mRNAs
# Input and output file paths
candidate_data = process_data(depth_path, ID_label)
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

print("finished")

