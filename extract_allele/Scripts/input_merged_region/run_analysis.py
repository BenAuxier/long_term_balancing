"""
Extract the allele of each gene in multiple genomes

"""
import os
import argparse
import subprocess

from prepare_alignment import prepare_analyze_alignment
from prepare_alignment import prepare_augustus_reference
from prepare_alignment import calculate_genome_number
from prepare_alignment import extract_annotations
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

def run_whole_analysis(reference_genome, species, augustus_species, type_annotation_ref, type_annotation_augustus, ID_ref_label, ID_augustus_label, key_words, base_path):
    # assembly_list, this file need to create manually
    assembly_list = f"{base_path}/genome_accessions/{species}.txt"
    #assembly_list = f"{base_path}/genome_accessions/{species}_test.txt"

    ##########################################################################
    # path to specific species
    main_path = f"{base_path}/{species}"
    assembly_dir = f"{main_path}/genome_assemblies"
    ref_path = f"{main_path}/reference_genome"
    ref_assembly = f"{ref_path}/{reference_genome}_genomic.fna"
    ref_fai = f"{ref_assembly}.fai"
    gff_refseq = f"{ref_path}/{reference_genome}_genomic.gff"
    gff_refseq_filtered = f"{gff_refseq[:-4]}_{type_annotation_ref}.gff"
    gff_augustus = f"{ref_path}/{reference_genome}_genomic_AUGUSTUS.gff"
    gff3_augustus = f"{ref_path}/{reference_genome}_genomic_AUGUSTUS.gff3"
    gff_augustus_filtered = f"{gff_augustus[:-4]}_{type_annotation_augustus}.gff"
    bam_path = f"{main_path}/alignment"
    bam_file = f"{main_path}/alignment/alignment_{species}.sorted.bam"
    filtered_region_nucleotides = f"{main_path}/depth_calculation/filtered_region_nucleotides.txt"

    ##########################################################################################
    """
    prepare_analyze_alignment(main_path, assembly_dir, ref_path, ref_assembly, gff_refseq, bam_path,
                              bam_file, reference_genome, assembly_list)
    """
    prepare_augustus_reference(ref_assembly, augustus_species)

    subprocess.run(["samtools", "faidx", ref_assembly])


    # filter the annotation with type_annotation
    gff_refseq_filtered = extract_annotations(gff_refseq, gff_refseq_filtered, type_annotation_ref, key_words)
    gff_augustus_filtered = extract_annotations(gff_augustus, gff_augustus_filtered, type_annotation_augustus, key_words)


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

    min_length_region = 200 # the minimal length of the genomic region merged in "calculate_depth_all"
    base_interval = "2" # the max interval between merged base
    min_overlap = 100 # the minimal length a candidate gene overlap with balancing selection region

    identity_visualization = "0.3" # the default identity threshold in visualization

    ##########################################################################
    # analyze the depth of the genomic regions
    gene_depth = f"{main_path}/depth_calculation/mean_depth_gene.txt"
    gene_region_depth = f"{main_path}/depth_calculation/mean_depth_region.txt"
    base_depth = f"{main_path}/depth_calculation/tmp/single_nucleotides_depth.txt"
    """"""
    calculate_depth_all(gene_depth, gene_region_depth, bam_file, main_path, gff_augustus_filtered,
                            lower_limit, upper_limit, base_interval, min_length_region)


    # load annotation data from gff annotation
    annotation_refseq, annotation_dict_refseq = load_annotation_refseq(gff_refseq_filtered, ID_ref_label, type_annotation_ref)
    annotation_augustus, annotation_dict_augustus = load_annotation_augustus(gff_augustus_filtered, ID_augustus_label, type_annotation_augustus)

    #print(annotation_refseq)
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

    #test
    #candidate_merge = {"NC_007196.1": candidate_merge["NC_007196.1"]}

    print("analyzing candidate data")
    output_json = f"{main_path}/temp/candidate_data_summary.json"
    """"""
    candidate_data_summary = analyze_all_candidate_position(candidate_merge, annotation_augustus, gene_depth_data,
                                                            bam_file, assembly_list, up_num, down_num, lower_limit,
                                                            minimal_alignment, type_annotation_ref)
    # save tmp file for candidate_data_summary
    save_json(candidate_data_summary, output_json)


    # reload candidate_data_summary
    candidate_data_summary = load_json(output_json)

    # find genes within the genomic regions
    results_path = f"{main_path}/results"
    os.makedirs(results_path, exist_ok=True)
    region_output_file = f"{results_path}/genomic_region_genes.csv"
    # candidate genes output
    go_path = f"{main_path}/go_analysis"
    os.makedirs(go_path, exist_ok=True)
    refseq_candidate_file = f"{go_path}/all_candidate_genes.txt"

    summary_genes = find_region_gene(candidate_data_summary, annotation_refseq, annotation_augustus, region_output_file,
                                     refseq_candidate_file)

    ## save final candidate genes to an Excel file
    print("saving candidate genes")
    result_file = f"{results_path}/{species}_final_candidates.xlsx"

    extract_candidates(candidate_data_summary, result_file, candidate_data,
                       genome_num,annotation_refseq, CDS_dict, ref_fai)

    #######################################################################
    statistics_path = f"{main_path}/statistics_data"
    os.makedirs(statistics_path, exist_ok=True)

    prepare_statistic_data(base_depth, statistics_path, gene_depth, gene_region_depth, result_file,
                           region_output_file, assembly_list, refseq_candidate_file)

    #######################################################################
    ## extract sequences for clinker visualization
    print("saving candidate genes")
    sequence_visualization = f"{main_path}/extract_sequences/clinker_visualization"
    clinker_output_path = f"{results_path}/clinker_results"
    sequence_interpro = f"{main_path}/extract_sequences/clinker_interpro"

    """"""
    extract_sequences(candidate_data_summary, reference_genome, gff_refseq, gff3_augustus,
                    sequence_visualization, extend, ref_assembly, assembly_dir, assembly_num, augustus_species)
    
    # run clinker for visualization
    print("running clinker")
    run_clinker_visualization(sequence_visualization, clinker_output_path, identity_visualization)

    ######################################################################################
    # interpro
    interpro_analysis_path = f"{main_path}/interpro_analysis"
    clinker_data_path = f"{interpro_analysis_path}/clinker_comparison"
    transformed_data_path = f"{interpro_analysis_path}/transformed_clinker_data"
    protein_path = f"{sequence_visualization}/protein_extraction"
    annotation_path = f"{protein_path}/interpro_annotation"
    summary_out_file = f"{main_path}/results/{species}_final_candidates.xlsx"
    final_output = f"{results_path}/interpro_annotation.xlsx"

    # similarity
    high_threshold = 0.85
    basic_threshold = 0.4
    low_threshold = 0.6
    # query the upstream reference gene.
    ratio_threhold = 0.6
    #####################################################################################

    additional_assembly = assembly_num_interpro - assembly_num
    """"""
    ## extract sequences for interpro analysis
    extract_sequences_interpro(sequence_interpro, candidate_data_summary,
                               sequence_visualization, extend, assembly_dir,
                               additional_assembly, augustus_species)
    
    # run clinker for comparison data
    run_clinker_data(sequence_interpro, clinker_data_path, "0.01")

    analysis_interpro(clinker_data_path, transformed_data_path, sequence_interpro, protein_path, high_threshold,
                      basic_threshold, ratio_threhold, annotation_path, summary_out_file, final_output)
    """"""

if __name__ == "__main__":
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

    run_whole_analysis(
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






