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
                              bam_file, reference_genome, assembly_list, augustus_species)
    
    # filter the annotation with type_annotation
    gff_refseq_filtered = extract_annotations(gff_refseq, gff_refseq_filtered, type_annotation_ref, key_words)
    gff_augustus_filtered = extract_annotations(gff_augustus, gff_augustus_filtered, type_annotation_augustus, key_words)
    """
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
    base_interval = "2" # the interval between merged base
    min_overlap = 100 # the minimal length a candidate gene overlap with balancing selection region

    similarity_visualization = "0.3"

    ##########################################################################

    # analyze the depth of the genomic regions
    gene_depth = f"{main_path}/depth_calculation/mean_depth_gene.txt"
    gene_region_depth = f"{main_path}/depth_calculation/mean_depth_region.txt"
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
    CDS_dict = csv_to_dict(id_dict_file, "protein_id", "gene_id")

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
                       genome_num,annotation_refseq, CDS_dict)

    #######################################################################
    ## extract sequences for clinker visualization
    print("saving candidate genes")
    sequence_visualization = f"{main_path}/extract_sequences/clinker_visualization"
    clinker_output_path = f"{results_path}/clinker_results"
    sequence_interpro = f"{main_path}/extract_sequences/clinker_interpro"
    clinker_data_path = f"{results_path}/clinker_comparison"

    """"""
    extract_sequences(candidate_data_summary, reference_genome, gff_refseq, gff3_augustus,
                    sequence_visualization, extend, ref_assembly, assembly_dir, assembly_num, augustus_species)
    
    # run clinker for visualization
    print("running clinker")
    run_clinker_visualization(sequence_visualization, clinker_output_path, similarity_visualization)

    #####################################################################################
    additional_assembly = assembly_num_interpro - assembly_num
    ## extract sequences for interpro analysis
    extract_sequences_interpro(sequence_interpro, candidate_data_summary,
                               sequence_visualization, extend, assembly_dir,
                               additional_assembly, augustus_species)

    # run clinker for comparison data
    run_clinker_data(sequence_interpro, clinker_data_path, "0.01")

    ######################################################################################
    #interpro
    from load_clinker_csv import analysis_interpro
    interpro_analysis_path = f"{main_path}/interpro_analysis"

    comparison_data_path = f"{interpro_analysis_path}/clinker_comparison"
    transformed_data_path = f"{interpro_analysis_path}/transformed_clinker_data"
    protein_path = f"{sequence_visualization}/protein_extraction"
    annotation_path = f"{protein_path}/interpro_annotation"
    excel_file = f"{main_path}/results/{species}_final_candidates.xlsx"
    final_output = f"{results_path}/interpro_annotation.xlsx"

    # similarity
    high_threshold = 0.85
    basic_threshold = 0.4
    low_threshold = 0.6
    # query the upstream reference gene.
    ratio_threhold = 0.6

    analysis_interpro(comparison_data_path, transformed_data_path, sequence_visualization, protein_path, high_threshold,
                      basic_threshold, ratio_threhold, annotation_path, excel_file, final_output)

    """
    # align sequence onto reference sequence to doublecheck and debug
    realign_output_path = f"{results_path}/sequence_alignments"
    annotate_file_path(sequence_path, realign_output_path)

    print("finished")"""

if __name__ == "__main__":
    # information
    reference_genome = "GCF_000184455.2"  # genome annotation should be GCF version (RefSeq)
    species = "aspergillus_oryzae"  # the species
    augustus_species = "aspergillus_oryzae"  # the reference species used in AUGUSTUS

    type_annotation_ref = "mRNA"  # type of annotation used in depth calculation, the third column
    ID_label = "transcript_id"  # this is the key that the gene/mRNA id follows
    key_words = None  # the keywords that have to be included in the annotation

    # file paths, including all input files
    base_path = "/lustre/BIF/nobackup/leng010/test"




