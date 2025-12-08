import subprocess
import os
import csv
from collections import defaultdict

from load_reference import create_ID_dictionary
from load_reference import csv_to_dict
from load_clinker_csv import interpro_annotation

####################################################3

def rename_fasta_headers(fasta_file, output_file, id_to_gene):
    """
    Replace the ID of each sequence in the FASTA file with the corresponding gene name.
    :param fasta_file: the FASTA file path
    :param output_file: output FASTA file path
    :param id_to_gene: Dictionary, {mRNA ID: Gene ID}
    """
    with open(fasta_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                seq_id = line.split("-")[1].strip()
                new_header = id_to_gene.get(seq_id, seq_id)  # If the ID is not found in the dictionary, retain the original ID.
                f_out.write(f">{new_header}\n")
            else:
                f_out.write(line)
    print(f"Renamed gene names in {fasta_file}")


def go_annotation(gff_refseq, ref_assembly, proteins_file_raw, proteins_file, annotation_file,mrna_dict):

    # extract protein sequences for all genes in the genome
    cmd = ["gffread", gff_refseq, "-g", ref_assembly, "-y", proteins_file_raw]
    #subprocess.run(cmd)
    print(f"Proteins are extracted from the genome {gff_refseq}")

    # rename the headers with gene_id
    #rename_fasta_headers(proteins_file_raw, proteins_file, mrna_dict)

    #make annotation with interpro
    #annotation_interpro = f"{annotation_file}_interpro"
    #interpro_annotation(proteins_file, annotation_interpro)

    # run eggnog to annotate files
    #subprocess.run(["conda", "activate", "eggnog"])
    #have to run conda activate eggnog before running cmd
    eggnog_path = "emapper.py"

    # make annotation with eggNOG-mapper
    cmd_emapper = [eggnog_path, "-i", proteins_file,"--output", annotation_file,
                   "--override",
                   "--cpu", "8", "--go_evidence", "all"]
    subprocess.run(cmd_emapper)

    return annotation_file

def extract_go_terms(eggnog_annotation, basic_name, go_basic_obo):

    # read GO -> namespace (BP/MF/CC)
    go_namespace = {}
    with open(go_basic_obo) as f:
        current_id = None
        for line in f:
            if line.startswith("id: GO:"):
                current_id = line.strip().split(" ")[1]
            elif line.startswith("namespace:"):
                ns = line.strip().split(": ")[1]
                go_namespace[current_id] = ns

    gene2go_all = defaultdict(list)
    gene2go_BP = defaultdict(list)
    gene2go_MF = defaultdict(list)
    gene2go_CC = defaultdict(list)

    with open(eggnog_annotation) as f:
        for line in f:

            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            gene = parts[0]
            gos = parts[9]

            if gos == "-" or gos == "":
                continue

            #print(gene, gos)

            go_list = gos.split(",")

            for go in go_list:
                gene2go_all[gene].append(go)
                if go in go_namespace:

                    if go_namespace[go] == "biological_process":
                        gene2go_BP[gene].append(go)
                    elif go_namespace[go] == "molecular_function":
                        gene2go_MF[gene].append(go)
                    elif go_namespace[go] == "cellular_component":
                        gene2go_CC[gene].append(go)


    # save files
    with open(f"{basic_name}.tsv", "w") as out:
        for gene, gos in gene2go_all.items():
            for go in gos:
                out.write(f"{gene}\t{go}\n")

    with open(f"{basic_name}_BP.tsv", "w") as out:
        for gene, gos in gene2go_BP.items():
            for go in gos:
                out.write(f"{gene}\t{go}\n")

    with open(f"{basic_name}_MF.tsv", "w") as out:
        for gene, gos in gene2go_MF.items():
            for go in gos:
                out.write(f"{gene}\t{go}\n")

    with open(f"{basic_name}_CC.tsv", "w") as out:
        for gene, gos in gene2go_CC.items():
            for go in gos:
                out.write(f"{gene}\t{go}\n")

    return True

def find_background_genes(gff_path, type_annotation_ref, ID_ref_label, background_output):
    background_gene = []
    with open(gff_path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t", fieldnames=[
            "seq_ID", "source", "type", "start", "end", "score", "strand", "phase", "attributes"
        ])
        for row in reader:
            if row["seq_ID"].startswith("#"):
                continue  # skip header/comment lines
            if row["type"] != type_annotation_ref:
                continue

            # parse attributes (GFF3 format: key=value;key=value;...)
            attr_dict = {}
            for attr in row["attributes"].split(";"):
                if attr.strip() == "":
                    continue
                if "=" in attr:
                    key, value = attr.strip().split("=", 1)
                    row[key] = value

            gene_id = row[ID_ref_label]
            background_gene.append(gene_id)

    background_gene = list(set(background_gene))
    with open(background_output, "w") as out:
        for gene in background_gene:
            out.write(f"{gene}\n")

    return background_gene


if __name__ == "__main__":

    #############################################################
    # path to specific species
    # basic information
    #reference_genome = "GCF_000002655.1"  # genome annotation should be GCF version
    #species = "aspergillus_fumigatus"

    reference_genome = "GCF_000184455.2"  # genome annotation should be GCF version (RefSeq)
    species = "aspergillus_oryzae"  # the species

    ID_ref_label = "locus_tag"
    type_annotation_ref = "mRNA"

    # file paths, including all input files
    base_path = "/lustre/BIF/nobackup/leng010/test"
    main_path = f"{base_path}/{species}"
    ref_path = f"{main_path}/reference_genome"
    ref_assembly = f"{ref_path}/{reference_genome}_genomic.fna"
    gff_refseq = f"{ref_path}/{reference_genome}_genomic.gff"

    results_path = f"{main_path}/results"
    region_output_file = f"{results_path}/genomic_region_genes.csv"
    refseq_candidate_file = f"{results_path}/all_candidate_genes.txt"

    #####################################################################

    # dictionary between locus_tag and CDS (XM) ID
    id_dict_file = f"{ref_path}/{reference_genome}_id_dict.csv"
    create_ID_dictionary(gff_refseq, id_dict_file, ID_ref_label, "CDS")
    # optional: "protein_id", "gene_id", "mrna_id"
    mrna_dict = csv_to_dict(id_dict_file, "mrna_id", "gene_id")

    go_path = f"{main_path}/go_analysis"
    os.makedirs(go_path, exist_ok=True)
    proteins_file_raw = f"{go_path}/{reference_genome}_raw.fasta"
    proteins_file = f"{go_path}/{reference_genome}_all_proteins.fasta"
    annotation_file = f"{go_path}/{reference_genome}_annotation"

    go_annotation(gff_refseq, ref_assembly, proteins_file_raw, proteins_file, annotation_file, mrna_dict)

    # get the result of eggnog annotation
    eggnog_annotation = f"{go_path}/{reference_genome}_annotation.emapper.annotations"
    genetogo_path = f"{go_path}/gene2go"
    os.makedirs(genetogo_path, exist_ok=True)
    basic_name = f"{genetogo_path}/{reference_genome}_gene2go"
    # download from http://purl.obolibrary.org/obo/go.obo
    go_basic_obo = "/lustre/BIF/nobackup/leng010/dataset/go_term/go.obo"

    extract_go_terms(eggnog_annotation, basic_name, go_basic_obo)

    background_output = f"{go_path}/background_genes.txt"
    find_background_genes(gff_refseq, type_annotation_ref, ID_ref_label, background_output)



