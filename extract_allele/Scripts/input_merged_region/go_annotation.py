import subprocess
import os
from load_clinker_csv import interpro_annotation
def go_annotation(gff_refseq, ref_assembly, proteins_file, annotation_file):

    # extract protein sequences for all genes in the genome
    cmd = ["gffread", gff_refseq, "-g", ref_assembly, "-y", proteins_file]
    subprocess.run(cmd)

    #make annotation with interpro
    #annotation_output = interpro_annotation(proteins_file, annotation_file)

    # make annotation with eggNOG-mapper
    cmd_emapper = ["emapper.py", "-i", proteins_file,"--output", annotation_file,
                   "--cpu", "8", "--go_evidence", "all"]
    subprocess.run(cmd_emapper)

    return annotation_file

def extract_go_terms(annotation_file):





    return True



if __name__ == "__main__":
    # path to specific species
    # information
    reference_genome = "GCF_000002655.1"  # genome annotation should be GCF version
    species = "aspergillus_fumigatus"

    # file paths, including all input files
    base_path = "/lustre/BIF/nobackup/leng010/test"
    main_path = f"{base_path}/{species}"
    ref_path = f"{main_path}/reference_genome"
    ref_assembly = f"{ref_path}/{reference_genome}_genomic.fna"
    gff_refseq = f"{ref_path}/{reference_genome}_genomic.gff"

    results_path = f"{main_path}/results"
    region_output_file = f"{results_path}/genomic_region_genes.csv"
    refseq_candidate_file = f"{results_path}/all_candidate_genes.txt"

    go_path = f"{main_path}/go_analysis"
    os.makedirs(go_path, exist_ok=True)
    proteins_file = f"{main_path}/go_analysis/{reference_genome}_all_proteins.fasta"
    annotation_file = f"{main_path}/go_analysis/{reference_genome}_annotation"
    annotation_file = go_annotation(gff_refseq, ref_assembly, proteins_file, annotation_file)
    print(annotation_file)


