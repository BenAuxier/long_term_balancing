import pandas as pd
import os
import glob
import csv
import subprocess

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

    similar_genes = {}
    for target_assembly in assemblies:
        df_target = df_gene[df_gene["Target_genome"] == target_assembly]

        max_row = df_target.loc[df_target["Identity"].idxmax()]

        max_identity = max_row["Identity"]
        max_gene = max_row["Target_gene"]
        #print(max_identity, max_gene)

        if max_identity >= high_similarity_threshold:
            similar_genes[target_assembly] = max_gene

    return similar_genes

def search_ref_gene(query_gene, df_ref_ref, df_ref_diver, orientation, ratio_threhold, high_similarity_threshold):
    # Finding the reference gene that have high similarity with ref and diff assemblies.
    local_genes = df_ref_ref["Query_gene"].unique()

    if orientation == "upstream":
        query_genes_sorted = sorted(local_genes, reverse=True)
    elif orientation == "downstream":
        query_genes_sorted = sorted(local_genes, reverse=False)

    for i in range(len(query_genes_sorted)):
        gene = query_genes_sorted[i]
        if query_gene in gene:
            break
    candidate_ref_genes = query_genes_sorted[i+1:]

    num_ref_assemblies = len(df_ref_ref["Target_genome"].unique())
    num_diver_assemblies = len(df_ref_diver["Target_genome"].unique())

    round = 0
    for candidate_gene in candidate_ref_genes:
        round += 1
        if round > 4:
            break
        similar_genes_ref = count_similar_genes(df_ref_ref, candidate_gene, high_similarity_threshold)
        num_ref_genes = len(similar_genes_ref)
        ref_ratio = num_ref_genes / num_ref_assemblies

        similar_genes_diver = count_similar_genes(df_ref_diver, candidate_gene, high_similarity_threshold)
        num_diver_genes = len(similar_genes_diver)
        diver_ratio = num_diver_genes / num_diver_assemblies

        if ref_ratio >= ratio_threhold and diver_ratio >= ratio_threhold:
            # return similar up and down stream genes
            return similar_genes_ref, similar_genes_diver
        else:
            continue

    return False, False

def find_intermediate(df, assembly, up_ref_gene, down_ref_gene):

    # select the data of this genome assembly
    df_selected = df[df["Target_genome"] == assembly]
    local_genes = df_selected["Target_gene"].unique()
    query_genes_sorted = sorted(local_genes, reverse=False)
    inter_gene = []
    for i in range(len(query_genes_sorted)):
        gene = query_genes_sorted[i]
        if gene == up_ref_gene:
            start_pos = i+1
        if gene == down_ref_gene:
            end_pos = i
    if start_pos and end_pos:
        inter_gene = query_genes_sorted[start_pos:end_pos]

    return inter_gene


def analyze_ref_gene(genomic_region, upstream_gene, downstream_gene, df_ref_ref, df_ref_diver, ratio_threhold, high_similarity_threshold):
    # finding the assmeblies and genes to extract genes in between

    up_genes_ref, up_genes_diver = search_ref_gene(upstream_gene, df_ref_ref, df_ref_diver, "upstream",
                                                                ratio_threhold, high_similarity_threshold)
    down_genes_ref, down_genes_diver = search_ref_gene(downstream_gene, df_ref_ref, df_ref_diver,"downstream",
                                                                      ratio_threhold, high_similarity_threshold)
    extraction_info = []
    for up_assembly, up_gene in up_genes_ref.items():
        if up_assembly not in down_genes_ref.keys():
            continue
        assembly = up_assembly
        down_gene = down_genes_ref[assembly]
        inter_gene_ref = find_intermediate(df_ref_ref, assembly, up_gene, down_gene)

        case_info_ref = {
            "genomic_region": genomic_region,
            "assembly_ID": assembly,
            "assembly_label": "ref_allele",
            #"up_ref_gene": up_gene,
            #"down_ref_gene": down_gene,
            "included_genes": inter_gene_ref
        }
        extraction_info.append(case_info_ref)

    for up_assembly, up_gene in up_genes_diver.items():
        if up_assembly not in down_genes_diver.keys():
            continue
        assembly = up_assembly
        down_gene = down_genes_diver[assembly]
        inter_gene_diver = find_intermediate(df_ref_diver, assembly, up_gene, down_gene)

        case_info_diver = {
            "genomic_region": genomic_region,
            "assembly_ID": assembly,
            "assembly_label": "diver_allele",
            #"up_diver_gene": up_gene,
            #"down_diver_gene": down_genes_diver[assembly],
            "included_genes": inter_gene_diver
        }
        extraction_info.append(case_info_diver)

    return extraction_info


def find_genome_assembly_path(file_dir, genome, suffix):
    # search for the file including genome accession in its name
    assembly_pattern = os.path.join(file_dir, f"*{genome}*{suffix}")
    matches = glob.glob(assembly_pattern)
    if not matches:
        warning = f"Warning: No genome assembly found for {genome}"
        print(warning)
        # return warning
    genome_assembly_path = matches[0]
    return genome_assembly_path


def extract_aminoacid(annotation_dir, gene_info, output_file):
    """

    :param gene_info: example, {'assembly_ID': 'GCA_018804105.1', 'assembly_label': 'diver_allele',
    'inter_gene_diver': ['g7.t1.cds', 'g8.t1.cds']}
    :return:
    """
    genomic_region = gene_info["genomic_region"]
    assembly_id = gene_info["assembly_ID"]
    assembly_label = gene_info["assembly_label"]
    inter_gene_diver = gene_info["included_genes"]

    candidate_genes = []
    for gene_cds in inter_gene_diver:
        gene = gene_cds.split(".")[0]
        candidate_genes.extend([gene])

    gff_path = find_genome_assembly_path(annotation_dir, assembly_id, ".gff3")

    for gene in candidate_genes:
        start_header = f"start gene {gene}"
        end_header = f"end gene {gene}"
        with open(gff_path, "r", encoding="utf-8-sig") as f:
            in_target_gene = False
            in_protein_seq = False
            seq_lines = []

            for line in f:
                line = line.strip()
                if not line.startswith("#"):
                    continue  # skip annotation
                # Identifying the start of the target gene
                if start_header in line:
                    in_target_gene = True
                    continue
                # Identifying the end of the target gene
                if end_header in line:
                    in_target_gene = False
                    break  # finish

                if in_target_gene:
                    # Begin reading the protein sequence
                    line = line.lstrip("# ").strip()
                    if line.startswith("protein sequence = ["):
                        in_protein_seq = True
                        seq_line = line.split("[", 1)[1]  # remove
                        seq_lines.append(seq_line)
                        continue
                    elif in_protein_seq:
                        if line.endswith("]"):
                            seq_lines.append(line.rstrip("]"))  # remove
                            in_protein_seq = False
                        else:
                            seq_lines.append(line)

            protein_seq = "".join(seq_lines).replace(" ", "").replace("\n", "")

            # Write to output file, FASTA format
            with open(output_file, "a", encoding="utf-8") as out:
                out.write(f">{genomic_region}__{assembly_id}__{assembly_label}__{gene}\n")
                # Line breaks every 60 characters for easy viewing.
                for i in range(0, len(protein_seq), 60):
                    out.write(protein_seq[i:i + 60] + "\n")


def analyze_clinker_results(input_data, annotation_dir, sequence_file):
    """


    :param input_data: the transformed clinker data
    :param annotation_dir: path to the annotation directory
    :param sequence_file: output file
    :return:
    """
    df = read_transform_data(input_data)

    # read the up and down-stream candidate genes
    genomic_region = df.loc[0, "Genomic_region"]
    upstream_gene = df.loc[0, "Upstream_gene"].split(".")[0]
    downstream_gene = df.loc[0, "Downstream_gene"].split(".")[0]

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

    # check whether the candidate gene have high similarity with ref alleles
    # and low similarity with divergent alleles
    """similar_genes_ref = count_similar_genes(df_ref_ref, query_gene, high_similarity_threshold)
    num_ref = len(similar_genes_ref.keys())
    ratio_ref = num_ref / num_ref_assemblies"""

    # query the upstream reference gene.
    ratio_threhold = 0.6

    # find the gene that could be used as reference gene in reference genome
    info_extraction = analyze_ref_gene(genomic_region, upstream_gene, downstream_gene, df_ref_ref,
                     df_ref_diver, ratio_threhold, high_similarity_threshold)
    # [{'genomic_region': 'g3347.t1-g3348.t1', 'assembly_ID': 'GCA_018804105.1',
    # 'assembly_label': 'diver_allele',
    # 'included_genes': ['g7.t1.cds', 'g8.t1.cds']}, ...]

    # If the output file already exists, delete it
    if os.path.exists(sequence_file):
        os.remove(sequence_file)
    # extract amino acid sequences
    for gene_info in info_extraction:
        extract_aminoacid(annotation_dir, gene_info, sequence_file)
    print("Protein sequences of candidate genomic region", genomic_region, "are saved to", sequence_file)

    return True

def interpro_annotation(sequence_file, annotation_file):

    cmd_interpro = ["interproscan.sh", "-f", "tsv", "-t", "p", "-i", sequence_file,
                    "--goterms", "-appl", "Pfam", "-b", annotation_file]
    subprocess.run(cmd_interpro)

    annotation_file = f"{annotation_file}.tsv"
    print("Interpro annotation results are saved to", annotation_file)
    print("==========================================================")
    return annotation_file


def load_interpro(interpro_file):
    columns = [
        "seq_id", "md5", "length", "source", "domain_info", "domain_name",
        "start", "end", "evalue", "status", "date",
        "ipr_accession", "ipr_name", "go_terms", "extra"
    ]

    interpro_data = []

    with open(interpro_file, "r", encoding="utf-8-sig") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")  # parse TSV
            # make dictionary
            row_dict = dict(zip(columns, fields))
            row_dict["go_terms"] = [t.strip("|") for t in row_dict["go_terms"].split("(InterPro)") if t.strip("|")]

            interpro_data.append(row_dict)

    return interpro_data

def analyze_annotation(interpro_data):
    go_info = []
    for row in interpro_data:
        # only keep Pfam annotation
        if row["source"] != "Pfam":
            continue
        gene_info = row["seq_id"].split("__")
        genomic_region = gene_info[0]
        assembly_id = gene_info[1]
        assembly_label = gene_info[2]
        gene = gene_info[3]
        domain_name = row["domain_name"]
        ipr_name = row["ipr_name"]
        go_terms = row["go_terms"]

        all_info = {
            "genomic_region": genomic_region,
            "assembly_id": assembly_id,
            "assembly_label": assembly_label,
            "gene": gene,
            "domain_name": domain_name,
            "ipr_name": ipr_name,
            "go_terms": go_terms
        }
        go_info.append(all_info)

    return go_info


def analyze_all_region(transformed_data_path,sequence_path, protein_path, annotation_path):
    files = os.listdir(transformed_data_path)
    annotation_files = []
    for transformed_data in files:
        genomic_region = transformed_data.split("_")[0]

        # extract amino acid sequences for each assembly of this genomic region
        input_data = f"{transformed_data_path}/{genomic_region}_data_transformed.csv"
        annotation_dir = f"{sequence_path}/{genomic_region}"
        sequence_file = f"{protein_path}/{genomic_region}_protein.fasta"
        analyze_clinker_results(input_data, annotation_dir, sequence_file)

        # make interpro annotation
        annotation_file = f"{annotation_path}/{genomic_region}_protein"
        annotation_file = interpro_annotation(sequence_file, annotation_file)
        #annotation_file = f"{annotation_file}.tsv"

def output_annotation_results(annotation_path):
    files = os.listdir(annotation_path)
    annotation_files = []
    for annotation_data in files:
        # load annotation
        genomic_region = annotation_data.split("_")[0]
        annotation_file = f"{annotation_path}/{genomic_region}_protein"
        annotation_data = load_interpro(annotation_file)
        go_info = analyze_annotation(annotation_data)
        annotation_files.extend(go_info)



if __name__ == "__main__":
    # Convert to DataFrame
    # data of MAT1-2-4
    main_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus"
    sequence_path = f"{main_path}/extract_sequences"
    result_path = f"{main_path}/results"

    comparison_data_path = f"{result_path}/clinker_comparison"
    transformed_data_path = f"{result_path}/clinker_comparison/transformed_data"
    os.makedirs(transformed_data_path, exist_ok=True)
    #transform_clinker_results(comparison_data_path, transformed_data_path)

    protein_path = f"{result_path}/clinker_comparison/protein_extraction"
    os.makedirs(protein_path, exist_ok=True)
    annotation_path = f"{result_path}/interpro_annotation"
    os.makedirs(annotation_path, exist_ok=True)

    # extract sequences and prepare annotation using interpro
    analyze_all_region(transformed_data_path, sequence_path, protein_path, annotation_path)














