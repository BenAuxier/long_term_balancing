import pandas as pd


def extract_genes_by_keywords(
    excel_file: str,
    keyword_file: str,
    output_excel: str,
    region_col: str = "genomic_region"):
    """
    Extract all rows belonging to genomic regions that contain at least one keyword.
    :param excel_file:
    :param keyword_file:
    :param output_excel:
    :param region_col:
    :return:
    """

    # -----------------------------
    # Read Excel file
    # -----------------------------
    df = pd.read_excel(excel_file)

    # -----------------------------
    # Build dictionary:
    # {genomic_region: [row_indices]}
    # -----------------------------
    region_dict = {}
    for idx, region in df[region_col].items():
        region_dict.setdefault(region, []).append(idx)

    # -----------------------------
    # Read keywords from text file
    # -----------------------------
    with open(keyword_file, "r") as f:
        keywords = [line.strip() for line in f if line.strip()]

    # -----------------------------
    # Identify genomic regions to extract
    # -----------------------------
    selected_regions = set()

    for idx, row in df.iterrows():
        # Convert the entire row to a single string for keyword matching
        row_text = " ".join(map(str, row.values))

        for kw in keywords:
            if kw in row_text:
                selected_regions.add(row[region_col])
                break

    # -----------------------------
    # Collect all rows for selected regions
    # -----------------------------
    selected_indices = []
    for region in selected_regions:
        selected_indices.extend(region_dict.get(region, []))

    # Remove duplicates and keep original order
    selected_indices = sorted(set(selected_indices))

    # Write to new Excel file
    df_selected = df.loc[selected_indices]
    df_selected.to_excel(output_excel, index=False)

    print(f"Extraction completed: {len(df_selected)} rows written to {output_excel}")


if __name__ == "__main__":
    # Convert to DataFrame
    # data of MAT1-2-4
    species = "aspergillus_fumigatus"
    #species = "aspergillus_oryzae"
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

    excel_file = final_output
    domain_keyword_file = "/lustre/BIF/nobackup/leng010/test/het_domain.txt"
    output_excel = f"{result_path}/interpro_annotation_het_domain.xlsx"
    region_col = "genomic_region"

    extract_genes_by_keywords(excel_file, domain_keyword_file, output_excel, region_col)

    print("saved")