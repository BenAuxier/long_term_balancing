



if __name__ == "__main__":
    # Convert to DataFrame
    # data of MAT1-2-4
    species = "aspergillus_fumigatus"
    species = "aspergillus_oryzae"
    main_path = f"/lustre/BIF/nobackup/leng010/test/{species}"
    sequence_path = f"{main_path}/extract_sequences/clinker_interpro"
    result_path = f"{main_path}/results"

    comparison_data_path = f"{result_path}/clinker_comparison"
    transformed_data_path = f"{result_path}/clinker_comparison/transformed_data"
    protein_path = f"{sequence_path}/protein_extraction"
    annotation_path = f"{protein_path}/interpro_annotation"
    excel_file = f"{main_path}/results/{species}_final_candidates.xlsx"
    final_output = f"{result_path}/interpro_annotation.xlsx"

    #similarity
    high_threshold = 0.85
    basic_threshold = 0.4
    low_threshold = 0.6
    # query the upstream reference gene.
    ratio_threhold = 0.6


