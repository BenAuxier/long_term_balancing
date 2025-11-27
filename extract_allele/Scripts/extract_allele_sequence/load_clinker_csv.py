import pandas as pd
import os

def load_clinker_csv(file_path, output_path):
    data = []
    current_title = None
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()

            # Title row (identified by “vs”)
            if " vs " in line:
                current_title = line.split("vs") # g3347.t1-g3348.t1__GCA_020501325.1__JAIBUL010000135.1__1-31369__ref_allele

                query_sequence = current_title[0].strip().split("__")
                region_name = query_sequence[0]
                region_gene = region_name.split("-")
                if len(region_gene) == 1:
                    upstream_gene = region_gene[0]
                    downstream_gene = region_gene[0]
                elif len(region_gene) == 2:
                    upstream_gene = region_gene[0]
                    downstream_gene = region_gene[1]

                query_genome = query_sequence[1]
                query_seq = query_sequence[2]
                query_pos = query_sequence[3]
                query_label = query_sequence[4]

                target_sequence = current_title[1].strip().split("__")
                print(target_sequence)
                target_genome = target_sequence[1]
                target_seq = target_sequence[2]
                target_pos = target_sequence[3]
                target_label = target_sequence[4]

                continue

            # Skip separator lines and blank lines
            if not line or line.startswith("-") or line.startswith("Query"):
                continue

            # Parsing table rows
            parts = line.split()
            if len(parts) == 4:
                query, target, identity, similarity = parts
                info = {
                    "Comparison": current_title,
                    "Genomic_region": region_name,
                    "Upstream_gene": upstream_gene,
                    "Downstream_gene": downstream_gene,
                    "Query_genome": query_genome,
                    "Query_label": query_label,
                    "Target_genome": target_genome,
                    "Target_label": target_label,
                    "Query_gene": query,
                    "Target_gene": target,
                    "Identity": float(identity),
                    "Similarity": float(similarity)
                }

                data.append(info)

                # reverse the order of query and target data
                info_reverse = {
                    "Comparison": current_title,
                    "Genomic_region": region_name,
                    "Upstream_gene": upstream_gene,
                    "Downstream_gene": downstream_gene,
                    "Query_genome": target_genome,
                    "Query_label": target_label,
                    "Target_genome": query_genome,
                    "Target_label": query_label,
                    "Query_gene": target,
                    "Target_gene": query,
                    "Identity": float(identity),
                    "Similarity": float(similarity)
                }

                data.append(info_reverse)
    #print(data)
    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)
    print(f"csv saved to {output_path}")
    return df

def transform_clinker_results(input_path, output_path):

    os.makedirs(output_path, exist_ok=True)

    for txt_file in os.listdir(input_path):
        if txt_file.endswith(".txt"):
            txt_file_path = os.path.join(input_path, txt_file)
            output_file_name = f"{txt_file[:-4]}_transformed.csv"
            output_file_path = os.path.join(output_path, output_file_name)
            #print(txt_file_path, output_file_path)

            gene_data = load_clinker_csv(txt_file_path, output_file_path)


    print(f"transform complete, files are saved to {output_path}")

def read_transform_data(input_data):
    df = pd.read_csv(input_data)
    return df


def count_similar_genes(df, query_gene, high_similarity_threshold):
    # query a gene in df, count the number of assemblies that have a high similarity gene

    df_gene = df[df["Query_gene"].str.contains(query_gene, na=False)]

    for index, row in df_gene.iterrows():
        #print(index, row)
        continue

    assemblies = df_gene["Target_genome"].unique()

    similar_genes = []
    for target_assembly in assemblies:
        df_target = df_gene[df_gene["Target_genome"] == target_assembly]

        max_row = df_target.loc[df_target["Identity"].idxmax()]

        max_identity = max_row["Identity"]
        max_gene = max_row["Target_gene"]
        #print(max_identity, max_gene)

        if max_identity >= high_similarity_threshold:
            similar_genes.append({target_assembly: max_gene})
    print(similar_genes)

    return similar_genes


def analyze_clinker_results(input_data):
    df = read_transform_data(input_data)

    # read the up and down-stream candidate genes
    upstream_gene = df.loc[0, "Upstream_gene"]
    downstream_gene = df.loc[0, "Downstream_gene"]

    # all data between reference genome vs. reference and divergent allele
    # only select the data when query sequence is ref_AUGUSTUS
    df_ref_data = df[
        (df["Query_label"] == "ref_AUGUSTUS") &
        (df["Target_label"].isin(["ref_allele", "diff_allele"]))
        ]

    #extract the data for ref and divergent allele separately
    df_ref_ref = df_ref_data[df_ref_data["Target_label"] == "ref_allele"]
    df_ref_diver = df_ref_data[df_ref_data["Target_label"] == "diff_allele"]

    # number of ref and diff assemblies
    num_ref_assemblies = len(df_ref_ref["Target_genome"].unique())
    num_diver_assemblies = len(df_ref_diver["Target_genome"].unique())
    num_all_assemblies = num_ref_assemblies + num_diver_assemblies

    high_similarity_threshold = 0.85
    low_similarity_threshold = 0.6

    query_gene = upstream_gene.split(".")[0]

    # check whether the candidate gene have high similarity with ref genes and low similarity with diff
    """similar_genes_ref = count_similar_genes(df_ref_ref, query_gene, high_similarity_threshold)
    num_ref = len(similar_genes_ref)
    ratio_ref = num_ref / num_ref_assemblies

    similar_genes_diff = count_similar_genes(df_ref_diver, query_gene, high_similarity_threshold)
    num_diff = len(similar_genes_diff)
    ratio_diff = num_diff / num_ref_assemblies"""

    similar_genes_ref = count_similar_genes(df_ref_ref, query_gene, high_similarity_threshold)






    return True



if __name__ == "__main__":
    # Convert to DataFrame
    # data of MAT1-2-4
    result_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/results"
    file_path = f"{result_path}/clinker_results/data/g3347.t1-g3348.t1_data.csv"
    output_path = f"{result_path}/clinker_results/data/g3347.t1-g3348.t1_reload.csv"
    #df = load_clinker_csv(file_path, output_path)

    input_path = f"{result_path}/clinker_comparasion"
    output_path = f"{result_path}/clinker_comparasion/transformed_data"
    #transform_clinker_results(input_path, output_path)

    input_data = f"{result_path}/clinker_comparasion/transformed_data/g3347.t1-g3348.t1_data_transformed.csv"
    analyze_clinker_results(input_data)




