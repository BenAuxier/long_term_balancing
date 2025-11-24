"""
Extract the allele of each gene in multiple genomes

"""
from run_analysis import run_whole_analysis


# information
reference_genome = "GCF_000002655.1" # genome annotation should be GCF version
species = "aspergillus_fumigatus"
augustus_species = "aspergillus_fumigatus"

type_annotation_ref = "mRNA" # type of annotation used in depth calculation, the third column
type_annotation_augustus = "transcript" # type of annotation used in depth calculation, the third column

key_words = None # the keywords that have to be included in the annotation
ID_ref_label = "locus_tag" #
ID_augustus_label = "gene_id" # this is the key that the gene/mRNA id follows in gff file


# file paths, including all input files
base_path = "/lustre/BIF/nobackup/leng010/test"

run_whole_analysis(reference_genome, species, augustus_species, type_annotation_ref, type_annotation_augustus, ID_ref_label, ID_augustus_label, key_words, base_path)

