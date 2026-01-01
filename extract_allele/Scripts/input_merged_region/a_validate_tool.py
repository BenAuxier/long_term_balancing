"""
Extract the allele of each gene in multiple genomes

"""
import os
import argparse
import subprocess
import csv

from prepare_alignment import prepare_analyze_alignment
from prepare_alignment import prepare_augustus_reference
from prepare_alignment import calculate_genome_number
from prepare_alignment import extract_annotations
from prepare_alignment import align_assemblies_to_reference
from depth_calculation import calculate_depth_all
from merge_region import process_merging
from merge_region import process_data_augustus
from load_reference import count_gff_features
from load_reference import create_ID_dictionary
from load_reference import csv_to_dict
from load_reference import load_annotation_refseq
from load_reference import load_annotation_augustus
from analyze_position import analyze_all_candidate_position
from analyze_position import save_json
from analyze_position import load_json
from make_outputs import find_region_gene
from make_outputs import extract_candidates
from make_outputs import extract_sequences
from make_outputs import extract_sequences_interpro
from visualization_clinker import run_clinker_visualization
from visualization_clinker import run_clinker_data
from doublecheck_alignment import annotate_file_path
from load_clinker_csv import analysis_interpro
from sta_summary import prepare_statistic_data
from concurrent.futures import ThreadPoolExecutor, as_completed


import random
import warnings

def read_non_empty_lines(txt_path):
    """
    Read non-empty lines from a txt file and return a list.
    """
    lines = []
    with open(txt_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line:
                lines.append(line)
    return lines


def sample_nonempty_lines(
    input_txt,
    output_txt,
    n,
    seed=None
):
    """
    Randomly sample n non-empty lines from a txt file.

    Parameters
    ----------
    input_txt : str
        Input txt file path
    output_txt : str
        Output txt file path
    n : int
        Number of lines to sample
    seed : int or None
        Random seed for reproducibility
    """

    if seed is not None:
        random.seed(seed)

    # read the lines
    with open(input_txt) as f:
        lines = [line.strip() for line in f if line.strip()]

    total = len(lines)

    if total == 0:
        raise ValueError("Input file contains no non-empty lines.")

    if n >= total:
        warnings.warn(
            f"Requested n={n} exceeds available lines ({total}); exporting all lines."
        )
        selected = lines
    else:
        selected = random.sample(lines, n)

    # write output
    with open(output_txt, "w") as out:
        for line in selected:
            out.write(line + "\n")

def loop_whole_analysis_old(genome_number, replication, reference_genome, species, augustus_species, type_annotation_ref, type_annotation_augustus,
                       ID_ref_label, ID_augustus_label, key_words, base_path):
    # with all assemblies
    assembly_list = f"{base_path}/genome_accessions/{species}.txt"

    test_path = f"/lustre/BIF/nobackup/leng010/test/{species}_test"
    os.makedirs(test_path, exist_ok=True)

    for number_test in genome_number:
        # create new path
        test_path_case = f"{test_path}/{number_test}_genomes"
        os.makedirs(test_path_case, exist_ok=True)

        for i in range(replication):
            replication_number = i + 1

            print(f"Analyzing [{number_test}] genome assemblies with replication of [{replication}].")

            path_replication = f"{test_path_case}/replication_{replication_number}"
            os.makedirs(path_replication, exist_ok=True)

            assembly_random = f"{path_replication}/{species}_{number_test}_{replication_number}.txt"

            sample_nonempty_lines(assembly_list,assembly_random,number_test,None)

            run_whole_analysis(assembly_random, reference_genome, species,
                               augustus_species, type_annotation_ref,
                               type_annotation_augustus,
                               ID_ref_label, ID_augustus_label,
                               key_words, base_path, path_replication)

            alignment_path = f"{path_replication}/alignment"

            cmd = ["rm", "-r", alignment_path]
            subprocess.run(cmd)


def run_one_case(number_test, replication_number,
                 test_path, replication, species,
                 assembly_list, reference_genome,
                 augustus_species, type_annotation_ref,
                 type_annotation_augustus,
                 ID_ref_label, ID_augustus_label,
                 key_words, base_path):

    print(f"Analyzing [{number_test}] genome assemblies "
          f"with replication [{replication_number}/{replication}].")

    # create paths
    test_path_case = f"{test_path}/{number_test}_genomes"
    os.makedirs(test_path_case, exist_ok=True)

    path_replication = f"{test_path_case}/replication_{replication_number}"
    os.makedirs(path_replication, exist_ok=True)

    assembly_random = f"{path_replication}/{species}_{number_test}_{replication_number}.txt"

    sample_nonempty_lines(
        assembly_list, assembly_random, number_test, None
    )

    run_whole_analysis(
        assembly_random, reference_genome, species,
        augustus_species, type_annotation_ref,
        type_annotation_augustus,
        ID_ref_label, ID_augustus_label,
        key_words, base_path, path_replication
    )

    alignment_path = f"{path_replication}/alignment"
    subprocess.run(["rm", "-r", alignment_path], check=False)

def loop_whole_analysis(genome_number, replication, reference_genome, species, augustus_species, type_annotation_ref, type_annotation_augustus,
                       ID_ref_label, ID_augustus_label, key_words, base_path):
    # with all assemblies
    assembly_list = f"{base_path}/genome_accessions/{species}.txt"

    test_path = f"/lustre/BIF/nobackup/leng010/test/{species}_test"
    os.makedirs(test_path, exist_ok=True)

    max_threads = 4

    tasks = []

    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        for number_test in genome_number:
            for i in range(replication):
                replication_number = i + 1

                future = executor.submit(
                    run_one_case,
                    number_test, replication_number,
                    test_path, replication, species,
                    assembly_list, reference_genome,
                    augustus_species, type_annotation_ref,
                    type_annotation_augustus,
                    ID_ref_label, ID_augustus_label,
                    key_words, base_path
                )
                tasks.append(future)

        for future in as_completed(tasks):
            try:
                future.result()
            except Exception as e:
                print("Task failed:", e)

def run_whole_analysis(assembly_random, reference_genome, species, augustus_species,
                       type_annotation_ref, type_annotation_augustus,
                       ID_ref_label, ID_augustus_label, key_words, base_path, working_path):

    ##########################################################################
    # path to specific species
    main_path_old = f"{base_path}/{species}"
    assembly_dir = f"{main_path_old}/genome_assemblies"
    ref_path = f"{main_path_old}/reference_genome"
    ref_assembly = f"{ref_path}/{reference_genome}_genomic.fna"
    ref_fai = f"{ref_assembly}.fai"
    gff_refseq = f"{ref_path}/{reference_genome}_genomic.gff"
    gff_refseq_filtered = f"{gff_refseq[:-4]}_{type_annotation_ref}.gff"
    gff_augustus = f"{ref_path}/{reference_genome}_genomic_AUGUSTUS.gff"
    gff3_augustus = f"{ref_path}/{reference_genome}_genomic_AUGUSTUS.gff3"
    gff_augustus_filtered = f"{gff_augustus[:-4]}_{type_annotation_augustus}.gff"

    bam_path = f"{working_path}/alignment"
    bam_file = f"{working_path}/alignment/alignment_{species}.sorted.bam"
    filtered_region_nucleotides = f"{working_path}/depth_calculation/filtered_region_nucleotides.txt"

    ##########################################################################################

    # no need to download the assemblies again
    #alignment
    align_assemblies_to_reference(assembly_random, ref_assembly, assembly_dir, bam_path, bam_file)
    """
    prepare_augustus_reference(ref_assembly, augustus_species)

    # filter the annotation with type_annotation
    gff_refseq_filtered = extract_annotations(gff_refseq, gff_refseq_filtered, type_annotation_ref, key_words)
    gff_augustus_filtered = extract_annotations(gff_augustus, gff_augustus_filtered, type_annotation_augustus,
                                                key_words)
    """

    # calculate number of assembly used in alignment
    genome_num = calculate_genome_number(assembly_random)

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
    base_interval = "2"  # the max interval between merged base
    min_overlap = 100  # the minimal length a candidate gene overlap with balancing selection region

    identity_visualization = "0.3"  # the default identity threshold in visualization

    ##########################################################################
    # analyze the depth of the genomic regions
    gene_depth = f"{working_path}/depth_calculation/mean_depth_gene.txt"
    gene_region_depth = f"{working_path}/depth_calculation/mean_depth_region.txt"
    base_depth = f"{working_path}/depth_calculation/tmp/single_nucleotides_depth.txt"
    """"""
    calculate_depth_all(gene_depth, gene_region_depth, bam_file, working_path, gff_augustus_filtered,
                        lower_limit, upper_limit, base_interval, min_length_region)

    # load annotation data from gff annotation
    annotation_refseq, annotation_dict_refseq = load_annotation_refseq(gff_refseq_filtered, ID_ref_label,
                                                                       type_annotation_ref)
    annotation_augustus, annotation_dict_augustus = load_annotation_augustus(gff_augustus_filtered, ID_augustus_label,
                                                                             type_annotation_augustus)

    # print(annotation_refseq)
    CDS_dict = False
    # dictionary between locus_tag and CDS (XM) ID
    id_dict_file = f"{ref_path}/{reference_genome}_id_dict.csv"
    create_ID_dictionary(gff_refseq, id_dict_file, ID_ref_label, "CDS")
    # optional: "protein_id", "gene_id", "mrna_id"
    CDS_dict = csv_to_dict(id_dict_file, "gene_id", "protein_id")

    # optional: "protein_id", "gene_id", "mrna_id"
    mrna_dict = csv_to_dict(id_dict_file, "mrna_id", "gene_id")

    # processes the input candidate mRNAs
    print("loading candidate data")
    gene_depth_data = process_data_augustus(gene_depth, ID_augustus_label)

    print("merging candidate data")
    candidate_data, candidate_merge = process_merging(gene_region_depth, ID_augustus_label, lower_limit, upper_limit,
                                                      annotation_dict_augustus, min_length_gene, min_overlap)

    # test
    # candidate_merge = {"NC_007196.1": candidate_merge["NC_007196.1"]}

    print("analyzing candidate data")
    output_json = f"{working_path}/temp/candidate_data_summary.json"
    """"""
    candidate_data_summary = analyze_all_candidate_position(candidate_merge, annotation_augustus, gene_depth_data,
                                                            bam_file, assembly_random, up_num, down_num, lower_limit,
                                                            minimal_alignment, type_annotation_ref)
    # save tmp file for candidate_data_summary
    save_json(candidate_data_summary, output_json)

    # reload candidate_data_summary
    candidate_data_summary = load_json(output_json)

    # find genes within the genomic regions
    results_path = f"{working_path}/results"
    os.makedirs(results_path, exist_ok=True)
    region_output_file = f"{results_path}/genomic_region_genes.csv"
    # candidate genes output
    go_path = f"{working_path}/go_analysis"
    os.makedirs(go_path, exist_ok=True)
    refseq_candidate_file = f"{results_path}/all_candidate_genes.txt"

    summary_genes = find_region_gene(candidate_data_summary, annotation_refseq, annotation_augustus, region_output_file,
                                     refseq_candidate_file)

    ## save final candidate genes to an Excel file
    print("saving candidate genes")
    result_file = f"{results_path}/{species}_final_candidates.xlsx"
    extract_candidates(candidate_data_summary, result_file, candidate_data,
                       genome_num, annotation_refseq, CDS_dict, ref_fai)

    #######################################################################
    statistics_path = f"{working_path}/statistics_data"
    os.makedirs(statistics_path, exist_ok=True)

    prepare_statistic_data(base_depth, statistics_path, gene_depth, gene_region_depth, result_file,
                           region_output_file, assembly_random, refseq_candidate_file)

def collect_gene_data(gene_num, genome_number, replication, test_path, complete_gene_list, output_file):


    BS_genes = read_non_empty_lines(complete_gene_list)
    num_BS_genes = len(BS_genes)
    all_statis = []

    for number_test in genome_number:
        test_path_case = f"{test_path}/{number_test}_genomes"

        for i in range(replication):
            replication_number = i + 1
            path_replication = f"{test_path_case}/replication_{replication_number}"
            gene_test = f"{path_replication}/statistics_data/all_candidate_genes.txt"
            test_list = read_non_empty_lines(gene_test)

            # number of genes found in the replication
            detected_genes = len(test_list)
            intersection_true = list(set(BS_genes) & set(test_list))

            # True results
            true_positive = len(intersection_true)

            # False negative
            only_in_all = list(set(BS_genes) - set(test_list))
            false_negative = len(only_in_all)

            # False positive
            only_in_case = list(set(test_list) - set(BS_genes))
            false_positive = len(only_in_case)

            #True negative
            true_negative = gene_num - num_BS_genes

            # sensitivity
            sensitivity = true_positive/num_BS_genes

            # False Negative Rate (FNR)
            fnr = false_negative/num_BS_genes

            # False discovery rate (FDR)
            fdr = false_positive / num_BS_genes

            # False Positive Rate (FPR)
            fpr = false_positive / (false_positive + true_negative)

            all_statis.append([number_test, replication_number,
                               true_positive, false_positive, true_negative, false_negative,
                               sensitivity, fnr, fdr, fpr])

    header = ["Assembly_number", "Replication",
              "TP", "FP", "TN", "FN",
              "Sensitivity", "FNR", "FDR", "FPR"]

    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(all_statis)

    print(f"Summary are saved to {output_file}")


if __name__ == "__main__":
    """"""
    parser = argparse.ArgumentParser(
        description="Run whole genome annotation analysis"
    )

    # information
    parser.add_argument(
        "--reference_genome",
        required=True,
        help="Reference genome accession (GCF version)"
    )
    parser.add_argument(
        "--species",
        required=True,
        help="Species name"
    )
    parser.add_argument(
        "--augustus_species",
        required=True,
        help="Augustus species name"
    )

    parser.add_argument(
        "--type_annotation_ref",
        default="mRNA",
        help="Annotation type for reference (3rd column in GFF)"
    )
    parser.add_argument(
        "--type_annotation_augustus",
        default="transcript",
        help="Annotation type for Augustus (3rd column in GFF)"
    )

    parser.add_argument(
        "--ID_ref_label",
        default="locus_tag",
        help="ID label key in reference GFF"
    )
    parser.add_argument(
        "--ID_augustus_label",
        default="gene_id",
        help="ID label key in Augustus GFF"
    )

    parser.add_argument(
        "--key_words",
        nargs="*",
        default=None,
        help="Keywords that must be included in annotation (optional)"
    )

    # paths
    parser.add_argument(
        "--base_path",
        required=True,
        help="Base directory containing all input files"
    )

    args = parser.parse_args()



    genome_number = [10, 15, 20, 30, 50, 75, 100, 150, 200, 250, 300]
    # test
    #genome_number = [5, 10, 20, 40, 80, 160]
    replication = 12

    """
    loop_whole_analysis(
        genome_number, replication,
        args.reference_genome,
        args.species,
        args.augustus_species,
        args.type_annotation_ref,
        args.type_annotation_augustus,
        args.ID_ref_label,
        args.ID_augustus_label,
        args.key_words,
        args.base_path
    )
    """

    species = args.species
    reference_genome = args.reference_genome

    test_path = f"/lustre/BIF/nobackup/leng010/test/{species}_test"
    complete_gene_list = f"/lustre/BIF/nobackup/leng010/test/{species}/statistics_data/all_candidate_genes.txt"

    test_summary = f"{test_path}/{species}_test_summary.csv"

    main_path = f"/lustre/BIF/nobackup/leng010/test/{species}"

    ref_path = f"{main_path}/reference_genome"
    ref_assembly = f"{ref_path}/{reference_genome}_genomic.fna"
    ref_fai = f"{ref_assembly}.fai"
    gff_refseq = f"{ref_path}/{reference_genome}__AUGUSTUS_transcript.gff"

    gff_refseq_filtered = f"{main_path}/reference_genome/{reference_genome}_genomic_{args.type_annotation_ref}.gff"

    refseq_genes = read_non_empty_lines(gff_refseq_filtered)
    num_refseq_genes = len(refseq_genes)
    print(num_refseq_genes)


    collect_gene_data(num_refseq_genes, genome_number, replication, test_path, complete_gene_list, test_summary)




