import pandas as pd



def load_clinker_csv(file_path):
    data = []
    current_title = None
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()

            # Title row (identified by “vs”)
            if " vs " in line:
                current_title = line.split("vs")
                quary_sequence = current_title[0].strip().split("__")
                region_name = quary_sequence[0]
                region_gene = region_name.split("-")
                if len(region_gene) == 1:
                    upstream_gene = region_gene[0]
                    downstream_gene = region_gene[0]
                elif len(region_gene) == 2:
                    upstream_gene = region_gene[0]
                    downstream_gene = region_gene[1]

                query_genome = quary_sequence[1]
                query_seq = quary_sequence[2]
                query_start = quary_sequence[3]
                query_end = quary_sequence[4]
                query_label = quary_sequence[5]

                target_sequence = current_title[1].strip()

                print(quary_sequence)

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
                    "Query": query,
                    "Target": target,
                    "Identity": float(identity),
                    "Similarity": float(similarity)
                })

    df = pd.DataFrame(data)
    #print(df)

if __name__ == "__main__":
    # Convert to DataFrame
    # data of MAT1-2-4
    result_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/results"
    file_path = f"{result_path}/clinker_results/data/g3347.t1-g3348.t1_data.csv"
    load_clinker_csv(file_path)

