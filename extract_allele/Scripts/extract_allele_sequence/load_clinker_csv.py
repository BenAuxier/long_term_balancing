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
                data.append({
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
                })
    print(data)
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
            print(txt_file_path, output_file_path)

            gene_data = load_clinker_csv(txt_file_path, output_file_path)


    print(f"transform complete, files are saved to {output_path}")


def analyze_clinker_results(input_path, output_path):



if __name__ == "__main__":
    # Convert to DataFrame
    # data of MAT1-2-4
    result_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/results"
    file_path = f"{result_path}/clinker_results/data/g3347.t1-g3348.t1_data.csv"
    output_path = f"{result_path}/clinker_results/data/g3347.t1-g3348.t1_reload.csv"
    df = load_clinker_csv(file_path, output_path)

    input_path = f"{result_path}/clinker_comparasion"
    output_path = f"{result_path}/clinker_comparasion/transformed_data"
    transform_clinker_results(input_path, output_path)

