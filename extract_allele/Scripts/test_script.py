"""
Extract the allele of each gene in multiple genomes

"""
from prepare_alignment import prepare_anallyze_alignment
from calculate_depth import run_bedtools_coverage
from merge_region import process_results
from merge_region import process_data
from load_reference_annotation import load_annotation
from analyze_position import analyze_all_candidate_position
from make_outputs import extract_outputs
from visualization_clinker import run_clinker_batch
import os

# information
reference_genome = "GCA_000184455.3"
species = "aspergillus_oryzae"
type_annotation = "gene" # type of annotation used in depth calculation
annotation_name = "locus_tag"

##########################################################################
# file paths, including all input files
# The modified csv file that contains the mRNA with depth of interests.
base_path = "/lustre/BIF/nobackup/leng010/test"
# path to specific species
main_path = f"{base_path}/{species}"
os.makedirs(main_path, exist_ok=True)


#assembly_list, this file need to create manually
assembly_list = f"{base_path}/genome_accessions/{species}.txt"

#assembly_dir, ref_assembly, ref_gff, gff_filtered, bam_path = prepare_anallyze_alignment(base_path, species, reference_genome, type_annotation,assembly_list)

#########################################################################
"""
assembly_dir = f"{main_path}/genome_assemblies" #path to all genome assemblies

# The txt file that includes the genome assemblies used in alignment.
# need to be created manually!
assembly_list = f"{main_path}/genome_accessions.txt 

ref_assembly # Reference genome assembly
ref_gff # Reference genome annotation
bam_path # Path to bam file
"""
assembly_dir= f"/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/genome_assemblies"
assembly_list= f"/lustre/BIF/nobackup/leng010/test/genome_accessions/aspergillus_oryzae.txt"
ref_assembly= "/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/genome_assemblies/reference_genome/GCA_000184455.3_genomic.fna"
ref_gff= "/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/genome_assemblies/reference_genome/GCA_000184455.3_genomic.gff"
gff_filtered= "/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/genome_assemblies/reference_genome/GCA_000184455.3_genomic_gene.gff"
bam_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/alignment/alignment_aspergillus_oryzae.sorted.bam"

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


print("minimal_alignment",minimal_alignment)



