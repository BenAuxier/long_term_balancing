"""
Extract the allele of each gene in multiple genomes

"""
from prepare_alignment import prepare_anallyze_alignment
from calculate_depth import run_bedtools_coverage
from mix_region import process_results
from mix_region import process_data
from load_reference_annotation import load_annotation
from analyze_position import analyze_all_candidate_position
from make_outputs import extract_outputs
from visualization_clinker import run_clinker_batch

# information
reference_genome = "GCF_000002655.1"
species = "aspergillus_fumigatus"
AUG_species = "aspergillus_fumigatus"
#type_annotation = "mRNA" # type of annotation used in depth calculation
#annotation_name = "transcript_id"
type_annotation = "gene"
annotation_name = "locus_tag"

##########################################################################
# file paths, including all input files
# The modified csv file that contains the mRNA with depth of interests.
base_path = "/lustre/BIF/nobackup/leng010/test"
# path to specific species
main_path = f"{base_path}/{species}"

#assembly_list, this file need to create manually
assembly_list = f"{base_path}/genome_accessions/{species}.txt"
assembly_list = f"{base_path}/genome_accessions/{species}_test.txt"
#########################################################################
#assembly_dir, ref_assembly, ref_gff, gff_filtered, bam_path = prepare_anallyze_alignment(base_path, species, reference_genome, type_annotation,assembly_list)


"""
assembly_dir = f"{main_path}/genome_assemblies" #path to all genome assemblies

# The txt file that includes the genome assemblies used in alignment.
# need to be created manually!
assembly_list = f"{main_path}/genome_accessions.txt 

ref_assembly # Reference genome assembly
ref_gff # Reference genome annotation
bam_path # Path to bam file
"""

assembly_dir= f"/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/genome_assemblies"
ref_assembly= "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/genome_assemblies/reference_genome/GCF_000002655.1_genomic.fna"
ref_gff= "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/genome_assemblies/reference_genome/GCF_000002655.1_genomic.gff"
gff_filtered= "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/genome_assemblies/reference_genome/GCF_000002655.1_genomic_gene.gff"
bam_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/alignment/alignment_aspergillus_fumigatus.sorted.bam"

##########################################################################################
# calculate number of assemblt used
with open(assembly_list, "r") as f:
    genome_num = sum(1 for line in f if line.strip())

# settings
up_num = 5
down_num = 5
assembly_num = 7
lower_limit = genome_num * 0.2
upper_limit = genome_num * 0.8
minimal_alignment = genome_num * 0.3
extend = 5000

##########################################################################
# path to store the depth file
depth_path = f"{main_path}/depth_calculation/{species}_{reference_genome}_meandepth.txt"
# analyze the depth of the genomic regions of
#run_bedtools_coverage(gff_filtered, bam_path, depth_path)

# load annotation data from gff annotation
annotation_sorted, annotation_sorted_dict = load_annotation(gff_filtered, type_annotation, annotation_name)

# processes the input candidate mRNAs
# Input and output file paths
candidate_data = process_data(depth_path)
candidate_merge = process_results(depth_path,lower_limit, upper_limit,annotation_sorted_dict)

# test the main code
candidate_data_test = dict(list(candidate_merge.items())[0:])

# print(candidate_data_test)
candidate_data_summary = analyze_all_candidate_position(candidate_data_test, annotation_sorted, candidate_data,
                        bam_path, assembly_list, up_num, down_num, lower_limit, minimal_alignment,annotation_name, type_annotation)

# extract sequences
candidates_path,sequence_path = extract_outputs(candidate_data_summary, reference_genome, ref_gff, main_path, extend, ref_assembly,
                          assembly_dir,assembly_num,candidate_data,species)

print(sequence_path)
#run clinker
clinker_output_dir = run_clinker_batch(sequence_path, main_path)

print("finished")
