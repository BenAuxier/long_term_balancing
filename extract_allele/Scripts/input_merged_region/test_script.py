"""
Extract the allele of each gene in multiple genomes

"""
import os
from prepare_alignment import prepare_analyze_alignment
from prepare_alignment import calculate_genome_number
from prepare_alignment import extract_annotations
from depth_calculation import calculate_depth_all
from merge_region import process_merging
from merge_region import process_data_augustus
from load_reference import count_gff_features
from load_reference import create_ID_dictionary
from load_reference import load_annotation_refseq
from load_reference import load_annotation_augustus
from analyze_position import analyze_all_candidate_position
from analyze_position import save_json
from analyze_position import load_json
from make_outputs import extract_candidates
from make_outputs import extract_sequences
from make_outputs import extract_sequences_interpro
from visualization_clinker import run_clinker_visualization
from visualization_clinker import run_clinker_data
from doublecheck_alignment import annotate_file_path


def run_whole_analysis(reference_genome, species, augustus_species, type_annotation_ref, type_annotation_augustus,
                       ID_ref_label, ID_augustus_label, key_words, base_path):
    # assembly_list, this file need to create manually
    assembly_list = f"{base_path}/genome_accessions/{species}.txt"
    # assembly_list = f"{base_path}/genome_accessions/{species}_test.txt"

    ##########################################################################
    # path to specific species
    main_path = f"{base_path}/{species}"
    assembly_dir = f"{main_path}/genome_assemblies"
    ref_path = f"{main_path}/reference_genome"
    ref_assembly = f"{ref_path}/{reference_genome}_genomic.fna"
    gff_refseq = f"{ref_path}/{reference_genome}_genomic.gff"
    gff_refseq_filtered = f"{gff_refseq[:-4]}_{type_annotation_ref}.gff"
    gff_augustus = f"{ref_path}/{reference_genome}_genomic_AUGUSTUS.gff"
    gff3_augustus = f"{ref_path}/{reference_genome}_genomic_AUGUSTUS.gff3"
    gff_augustus_filtered = f"{gff_augustus[:-4]}_{type_annotation_augustus}.gff"
    bam_path = f"{main_path}/alignment"
    bam_file = f"{main_path}/alignment/alignment_{species}.sorted.bam"
    filtered_region_nucleotides = f"{main_path}/depth_calculation/filtered_region_nucleotides.txt"

    ##########################################################################################

    # verify some basic details
    # check reference annotation .gff file
    count_gff_features(gff_augustus)

    # calculate number of assembly used in alignment
    genome_num = calculate_genome_number(assembly_list)

    # settings
    up_num = 5
    down_num = 5
    assembly_num = 7
    assembly_num_interpro = 15

    lower_limit = genome_num * 0.2
    upper_limit = genome_num * 0.8
    print(lower_limit, upper_limit)

    minimal_alignment = genome_num * 0.3
    extend = 5000
    min_length_gene = 300  # the minimal length of candidate gene
    transfer_id = True  # whether transfer genomic region name to CDS ID

    min_length_region = 200  # the minimal length of the genomic region merged in "calculate_depth_all"
    base_interval = "2"  # the interval between merged base
    min_overlap = 100  # the minimal length a candidate gene overlap with balancing selection region

    similarity_visualization = "0.3"

    ##########################################################################

    # analyze the depth of the genomic regions
    gene_depth = f"{main_path}/depth_calculation/mean_depth_gene.txt"
    gene_region_depth = f"{main_path}/depth_calculation/mean_depth_region.txt"

    # load annotation data from gff annotation
    annotation_refseq, annotation_dict_refseq = load_annotation_refseq(gff_refseq_filtered, ID_ref_label,
                                                                       type_annotation_ref)
    annotation_augustus, annotation_dict_augustus = load_annotation_augustus(gff_augustus_filtered, ID_augustus_label,
                                                                             type_annotation_augustus)

    output_json = f"{main_path}/temp/candidate_data_summary.json"
    candidate_data_summary = load_json(output_json)

    return candidate_data_summary, annotation_refseq, annotation_augustus

import csv
def find_region_gene(candidate_data_summary, annotation_refseq, annotation_augustus, region_output_file, refseq_candidate_file):
    """

    :param candidate_data_summary:

    [{
        "region_name": name,
        "candidate_data": candidate region, # information of this candidate, {seq: [merged_info1, ...], ...}
        "position_info": up_down_locations, # find upstream and downstream positions, dict
        "up_down_loci": all_up_down_loci, # the selected locations
    }, ......]

    [{"region_name": candidate["region_name"],
     "position_info": up_down_locations,
     "up_down_loci": all_up_down_loci,
    },...]
    :param annotation_refseq: {seq_ID: [row_dict, ...]}
    :param annotation_augustus: {seq_ID: [row_dict, ...]}
    :param output_path:
    :return:
    """
    summary_genes = []
    all_refseq_genes = []
    for summary in candidate_data_summary:
        region_name = summary["region_name"]
        region_info = summary["position_info"]["position"]
        region_seq = region_info["seq_ID"]
        region_start = region_info["start"]
        region_end = region_info["end"]

        annotation_refseq_seq = annotation_refseq[region_seq]
        annotation_augustus_seq = annotation_augustus[region_seq]

        refseq_gene = []
        augustus_gene = []

        for row in annotation_refseq_seq:
            candidate_start = row["start"]
            candidate_end = row["end"]
            candidate_id = row["id"]

            if (region_start <= candidate_start <= region_end) or (region_start <= candidate_end <= region_end):
                refseq_gene.append(candidate_id)
                all_refseq_genes.extend(refseq_gene)

        for row in annotation_augustus_seq:
            candidate_start = row["start"]
            candidate_end = row["end"]
            candidate_id = row["id"]

            if (region_start <= candidate_start <= region_end) or (region_start <= candidate_end <= region_end):
                augustus_gene.append(candidate_id)

        summary_genes.append([region_name, refseq_gene, augustus_gene])

    with open(region_output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["genomic_region", "RefSeq_genes", "Augustus_genes"])
        writer.writerows(summary_genes)
    print(f"Genomic region file written to {region_output_file}")

    all_refseq_genes = list(set(all_refseq_genes))
    with open(refseq_candidate_file, "w", newline="", encoding="utf-8") as f:
        for item in all_refseq_genes:
            f.write(str(item) + "\n")

    print(f"Refseq candidate genes written to {refseq_candidate_file}")

    return summary_genes


if __name__ == "__main__":
    # information
    reference_genome = "GCF_000002655.1"  # genome annotation should be GCF version
    species = "aspergillus_fumigatus"
    augustus_species = "aspergillus_fumigatus"

    type_annotation_ref = "mRNA"  # type of annotation used in depth calculation, the third column
    type_annotation_augustus = "transcript"  # type of annotation used in depth calculation, the third column

    key_words = None  # the keywords that have to be included in the annotation
    ID_ref_label = "locus_tag"  #
    ID_augustus_label = "gene_id"  # this is the key that the gene/mRNA id follows in gff file

    # file paths, including all input files
    base_path = "/lustre/BIF/nobackup/leng010/test"
    main_path = f"{base_path}/{species}"
    results_path = f"{main_path}/results"

    candidate_data_summary, annotation_refseq, annotation_augustus = run_whole_analysis(reference_genome, species, augustus_species, type_annotation_ref, type_annotation_augustus,
                       ID_ref_label, ID_augustus_label, key_words, base_path)

    region_output_file = f"{results_path}/genomic_region_genes.csv"
    refseq_candidate_file = f"{results_path}/all_candidate_genes.txt"

    summary_genes = find_region_gene(candidate_data_summary, annotation_refseq, annotation_augustus, region_output_file, refseq_candidate_file)







