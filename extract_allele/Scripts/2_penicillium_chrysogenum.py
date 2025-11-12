#penicillium_chrysogenum.py
"""
Extract the allele of each gene in multiple genomes

"""
from run_analysis import run_whole_analysis


# information, settings
reference_genome = "GCF_028827035.1" # genome annotation should be GCF version (RefSeq)
species = "penicillium_chrysogenum" # the species
augustus_species = "aspergillus_fumigatus" # the reference species used in AUGUSTUS

type_annotation = "mRNA" # type of annotation used in depth calculation, the third column
ID_label = "transcript_id" # this is the key that the gene/mRNA id follows
key_words = None # the keywords that have to be included in the annotation

# file paths, including all input files
base_path = "/lustre/BIF/nobackup/leng010/test"

run_whole_analysis(reference_genome, species, augustus_species, type_annotation, ID_label, key_words, base_path)