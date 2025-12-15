import pandas as pd
import os
import glob
import csv
import subprocess
import ast
from concurrent.futures import ProcessPoolExecutor, as_completed

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


def count_similar_genes(df, query_gene, identity_threshold):
    """

    :param df:
    :param query_gene:
    :param identity_threshold:
    :return:
    """
    # query a gene in df, count the number of assemblies that have a high similarity gene

    df_gene = df[df["Query_gene"].str.contains(query_gene, na=False)]

    for index, row in df_gene.iterrows():
        #print(index, row)
        continue

    assemblies = df_gene["Target_genome"].unique().tolist()

    similar_genes = {}

    # analyze each assembly
    for target_assembly in assemblies:
        df_target = df_gene[df_gene["Target_genome"] == target_assembly]

        # find the row (comparison) with the highest Identity
        max_row = df_target.loc[df_target["Identity"].idxmax()]

        # the value of identity and the related target gene
        max_identity = max_row["Identity"]
        max_gene = max_row["Target_gene"]
        #print(max_identity, max_gene)

        if max_identity >= identity_threshold:
            similar_genes[target_assembly] = max_gene

    return similar_genes

def search_ref_gene(query_gene, df_ref_ref, df_ref_diver, orientation, ratio_threhold, high_threshold, basic_threshold):
    """

    :param query_gene:
    :param df_ref_ref:
    :param df_ref_diver:
    :param orientation:
    :param ratio_threhold:
    :param high_threshold:
    :return:
    """

    # Finding the reference gene that have high similarity with ref and diver assemblies.
    local_genes = df_ref_ref["Query_gene"].unique().tolist()

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

    if num_ref_assemblies * num_diver_assemblies == 0:
        return False, False, False

    round = 0
    for candidate_gene in candidate_ref_genes:
        round += 1
        if round > 4:
            break

        #find the high similarity gene of this ref gene in ref and diver assemblies
        high_similar_ref = count_similar_genes(df_ref_ref, candidate_gene, high_threshold)
        num_ref_genes = len(high_similar_ref)
        ref_ratio = num_ref_genes / num_ref_assemblies

        high_similar_diver = count_similar_genes(df_ref_diver, candidate_gene, high_threshold)
        num_diver_genes = len(high_similar_diver)
        diver_ratio = num_diver_genes / num_diver_assemblies

        if ref_ratio < ratio_threhold or diver_ratio < ratio_threhold:
            # return similar up and down stream genes
            continue

        elif ref_ratio >= ratio_threhold and diver_ratio >= ratio_threhold:
            # select data with basic threshold (less similar)
            similar_genes_ref = count_similar_genes(df_ref_ref, candidate_gene, basic_threshold)
            similar_genes_diver = count_similar_genes(df_ref_diver, candidate_gene, basic_threshold)

            return similar_genes_ref, similar_genes_diver, round
            continue





    return False, False, False

def find_intermediate(df, assembly, up_ref_gene, down_ref_gene):

    # select the data of this genome assembly
    df_selected = df[df["Target_genome"] == assembly]
    local_genes = df_selected["Target_gene"].unique().tolist()

    # have not sorted, defaultly should be ok
    query_genes_sorted = local_genes

    # find intermediate genes
    inter_gene = []
    for i in range(len(query_genes_sorted)):
        gene = query_genes_sorted[i]
        if gene == up_ref_gene:
            start_pos = i+1
        if gene == down_ref_gene:
            end_pos = i

    if start_pos is not None and end_pos is not None and start_pos < end_pos:
        inter_gene = query_genes_sorted[start_pos:end_pos]
    else:
        inter_gene = ["no_gene"]

    return inter_gene


def analyze_ref_gene(genomic_region, upstream_gene, downstream_gene, df_ref_ref, df_ref_diver, ratio_threhold, high_threshold, basic_threshold):

    # finding the assmeblies and genes to extract genes in between

    up_genes_ref, up_genes_diver, up_rank = search_ref_gene(upstream_gene, df_ref_ref, df_ref_diver, "upstream",
                                                                ratio_threhold, high_threshold, basic_threshold)
    down_genes_ref, down_genes_diver, down_rank = search_ref_gene(downstream_gene, df_ref_ref, df_ref_diver,"downstream",
                                                                      ratio_threhold, high_threshold, basic_threshold)

    # if certain ref gene is not found:
    if not all([up_genes_ref, up_genes_diver, down_genes_ref, down_genes_diver]):
        return False

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
            "ref_included": f"up {up_rank-1}, down {down_rank-1}",
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
            "ref_included": f"up {up_rank-1}, down {down_rank-1}",
            #"up_diver_gene": up_gene,
            #"down_diver_gene": down_genes_diver[assembly],
            "included_genes": inter_gene_diver
        }
        extraction_info.append(case_info_diver)

    return extraction_info


def find_genome_assembly_path(file_dir, genome, suffix):
    """

    :param file_dir:
    :param genome:
    :param suffix:
    :return:
    """
    # search for the file including genome accession in its name
    assembly_pattern = os.path.join(file_dir, f"*{genome}*{suffix}")
    matches = glob.glob(assembly_pattern)
    if not matches:
        warning = f"Warning: No genome assembly found for {genome}"
        print(warning)
        return False
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

    if not gff_path:
        return False

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


def analyze_clinker_results(input_data, annotation_dir, sequence_file, high_threshold, basic_threshold, ratio_threhold):
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
        (df["Target_label"].isin(["ref_allele", "diver_allele"]))
        ]

    #extract the data for ref and divergent allele separately
    df_ref_ref = df_ref_data[df_ref_data["Target_label"] == "ref_allele"]
    df_ref_diver = df_ref_data[df_ref_data["Target_label"] == "diver_allele"]

    # number of ref and divergent assemblies
    num_ref_assemblies = len(df_ref_ref["Target_genome"].unique())
    num_diver_assemblies = len(df_ref_diver["Target_genome"].unique())
    num_all_assemblies = num_ref_assemblies + num_diver_assemblies

    query_gene = upstream_gene.split(".")[0]

    # check whether the candidate gene have high similarity with ref alleles
    # and low similarity with divergent alleles
    """similar_genes_ref = count_similar_genes(df_ref_ref, query_gene, high_similarity_threshold)
    num_ref = len(similar_genes_ref.keys())
    ratio_ref = num_ref / num_ref_assemblies"""

    # find the gene that could be used as reference gene in reference genome
    info_extraction = analyze_ref_gene(genomic_region, upstream_gene, downstream_gene, df_ref_ref,
                     df_ref_diver, ratio_threhold, high_threshold, basic_threshold)

    # [{'genomic_region': 'g3347.t1-g3348.t1', 'assembly_ID': 'GCA_018804105.1',
    # 'assembly_label': 'diver_allele',
    # 'included_genes': ['g7.t1.cds', 'g8.t1.cds']}, ...]

    if info_extraction == False:
        return False

    # If the output file already exists, delete it
    if os.path.exists(sequence_file):
        os.remove(sequence_file)

    # extract amino acid sequences
    for gene_info in info_extraction:
        extract_aminoacid(annotation_dir, gene_info, sequence_file)

    print("Protein sequences of candidate genomic region", genomic_region, "are saved to", sequence_file)

    return info_extraction

def interpro_annotation(sequence_file, annotation_file):
    """

    :param sequence_file:
    :param annotation_file:
    :return:
    """

    final_output = f"{annotation_file}.tsv"
    # if file exist, do not run interpro
    if os.path.exists(final_output):
        print(f"Annotation file already exists: {final_output}")
        print("Skip running InterProScan.")
        print("==========================================================")
        return final_output

    cmd_interpro = ["interproscan.sh", "-f", "tsv", "-t", "p", "-i", sequence_file,
                    "--goterms", "-appl", "Pfam", "-b", annotation_file]
    subprocess.run(cmd_interpro)

    print("Interpro annotation results are saved to", final_output)
    print("==========================================================")
    return final_output


def extract_all_sequences(transformed_data_path,sequence_path, protein_path, high_threshold, basic_threshold, ratio_threhold):
    """
    extract sequences for all candidates, and make dictionary including all extraction info
    :param transformed_data_path:
    :param sequence_path:
    :param protein_path:
    :param annotation_path:
    :param high_threshold:
    :param basic_threshold:
    :param ratio_threhold:
    :return:
    """
    os.makedirs(protein_path, exist_ok=True)
    files = os.listdir(transformed_data_path)
    candidate_info = []
    for transformed_data in files:
        genomic_region = transformed_data.split("_")[0]

        # extract amino acid sequences for each assembly of this genomic region
        input_data = f"{transformed_data_path}/{genomic_region}_data_transformed.csv"
        annotation_dir = f"{sequence_path}/{genomic_region}"
        sequence_file = f"{protein_path}/{genomic_region}_protein.fasta"
        info_extraction = analyze_clinker_results(input_data, annotation_dir, sequence_file, high_threshold, basic_threshold, ratio_threhold)

        # if there is no ideal reference gene, skip this gene
        if info_extraction == False:
            print(f"{genomic_region} did not find reference gene in clinker data, skipping.")
            print(f"==========================================================")
            continue

        # save candidate data information in dict
        candidate_info.extend(info_extraction)

        candidate_dict = {}
        for info in candidate_info:
            region = info["genomic_region"]
            assembly = info["assembly_ID"]

            # if genomic_region not exist
            if region not in candidate_dict:
                candidate_dict[region] = {}
            # save information
            candidate_dict[region][assembly] = info

    return candidate_dict

def interpro_all_sequences(protein_path, annotation_path):
    """

    :param protein_path:
    :param annotation_path:
    :return:
    """
    # perform interpro annotation for each candidate sequence
    os.makedirs(annotation_path, exist_ok=True)

    sequence_files = os.listdir(protein_path)
    tasks = []
    with ProcessPoolExecutor(max_workers=4) as executor:
        for sequence_file in sequence_files:
            sequence_path = os.path.join(protein_path, sequence_file)
            # Output file base name
            annotation_file = os.path.join(annotation_path, sequence_file.strip(".fasta"))

            # Submit the task (simply replace it with existing function).
            # interpro_annotation(sequence_path, annotation_file)
            future = executor.submit(interpro_annotation, sequence_path, annotation_file)
            tasks.append(future)

        # Output in order of completion
        for f in as_completed(tasks):
            result = f.result()
            print(f"[InterPro Done] {result}")

    # old version
    """for sequence_file in sequence_files:
        sequence_path = f"{protein_path}/{sequence_file}"
        annotation_file = f"{annotation_path}/{genomic_region}_protein"
        annotation_file = interpro_annotation(sequence_path, annotation_file)
        # annotation_file = f"{annotation_file}.tsv" """

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

            # only keep Pfam annotation
            if row_dict["source"] != "Pfam":
                continue

            interpro_data.append(row_dict)

    return interpro_data

def interpro_data_dict(interpro_data):

    # first extract data and make a dictionary
    interpro_dict = {}

    for row in interpro_data:
        gene_info = row["seq_id"].split("__")
        assembly_id = gene_info[1]
        gene = gene_info[3]

        if assembly_id not in interpro_dict:
            interpro_dict[assembly_id] = {}

        interpro_dict[assembly_id][gene] = {
            "domain_name": "unpredicted" if row["domain_name"] == "-" else row["domain_name"],
            "ipr_name": "unpredicted" if row["ipr_name"] == "-" else row["ipr_name"],
            "go_terms": "unpredicted" if row["go_terms"] == ["-"] else row["go_terms"]
        }

    return interpro_dict

def list_to_pipe_string(lst):
    """
    Convert a list of items to a string separated by ' | '.
    If the list is empty or None, return an empty string.
    """
    if len(lst) == 0:
        return ""
    list_string = " | ".join(str(item) for item in lst)
    return list_string

def analyze_annotation(assembly, included_genes, interpro_dict):
    """

    :param assembly:
    :param included_genes:
    :param interpro_dict:
    :return:
    """

    domain_info = []
    go_info = []
    ipr_info = []
    for gene in included_genes:
        gene_id = gene.split(".")[0]

        if gene_id not in interpro_dict[assembly].keys():
            domain_info.append("unpredicted")
            go_info.append("unpredicted")
            ipr_info.append("unpredicted")
            continue

        gene_dict = interpro_dict[assembly][gene_id]
        domain_name = gene_dict["domain_name"]
        ipr_name = gene_dict["ipr_name"]
        go_terms = gene_dict["go_terms"]

        domain_info.append(domain_name)
        ipr_info.append(ipr_name)
        go_info.append(go_terms)

    return domain_info, ipr_info, go_info

def build_region_to_genes(excel_file):
    """
    return: dictionary
    {
        genomic_region: [RefSeq_genes included in this region, in gene ID]
    }
    """
    df = pd.read_excel(excel_file)

    result = {}

    for _, row in df.iterrows():
        region = str(row["genomic_region"]).strip()

        # Convert the RefSeq_genes field
        genes = row["RefSeq_genes"]
        if isinstance(genes, str):
            try:
                genes = ast.literal_eval(genes)
            except:
                genes = []
        elif pd.isna(genes):
            genes = []

        # create key
        if region not in result:
            result[region] = []

        # add gene
        result[region].extend(genes)

    return result

def output_annotation_results(candidate_dict, annotation_path, id_dict, output_file):
    """

    :param candidate_dict:
    :param annotation_path:
    :param output_file:
    :return:
    """
    files = os.listdir(annotation_path)

    annotation_files = []
    for annotation_data in files:
        if not annotation_data.endswith(".tsv"):
            continue
        annotation_file = f"{annotation_path}/{annotation_data}"
        # load annotation
        genomic_region = annotation_data.split("_")[0]

        # test
        #if not genomic_region == "g5348.t1-g5349.t1":
            #continue

        region_info = candidate_dict[genomic_region]

        # load analysis results
        interpro_data = load_interpro(annotation_file)

        # if no prediction for this genomic region
        if not interpro_data:
            all_info = {
                "genomic_region": genomic_region,
                "RefSeq_genes": list_to_pipe_string(id_dict[genomic_region]),
                "ref_included": "",
                "assembly_id": "",
                "assembly_label": "",
                "gene": "",
                "domain_name": "unpredicted",
                "ipr_name": "unpredicted",
                "go_terms": "unpredicted"
            }
            annotation_files.append(all_info)
            continue

        #transform interpro info to dict
        interpro_dict = interpro_data_dict(interpro_data)

        for assembly, info in region_info.items():
            #Here can see which assemblies are present and which are absent.
            assembly_label = info["assembly_label"]
            included_genes = info["included_genes"]
            ref_included = info["ref_included"]

            # if no prediction for this assembly
            if assembly not in interpro_dict.keys():
                all_info = {
                    "genomic_region": genomic_region,
                    "RefSeq_genes": list_to_pipe_string(id_dict[genomic_region]),
                    "ref_included": ref_included,
                    "assembly_id": assembly,
                    "assembly_label": assembly_label,
                    "gene": list_to_pipe_string(included_genes),
                    "domain_name": "unpredicted",
                    "ipr_name": "unpredicted",
                    "go_terms": "unpredicted"
                }
                annotation_files.append(all_info)
                continue

            # load annotation info for the gene of this assembly
            domain_info, ipr_info, go_info = analyze_annotation(assembly, included_genes, interpro_dict)

            all_info = {
                "genomic_region": genomic_region,
                "RefSeq_genes": list_to_pipe_string(id_dict[genomic_region]),
                "ref_included": ref_included,
                "assembly_id": assembly,
                "assembly_label": assembly_label,
                "gene": list_to_pipe_string(included_genes),
                "domain_name": list_to_pipe_string(domain_info),
                "ipr_name": list_to_pipe_string(ipr_info),
                "go_terms": list_to_pipe_string(go_info)
            }
            annotation_files.append(all_info)

    # Create a DataFrame (keys automatically become column headers)
    df = pd.DataFrame(annotation_files)

    # Columns used as grouping keys
    keys = ["genomic_region", "RefSeq_genes", "ref_included", "assembly_label",
            "domain_name", "ipr_name", "go_terms"]

    # Combine groups and assembly
    df2 = (
        df.groupby(keys, dropna=False)
        .agg({
            "assembly_id": list,
            "gene": list
        })
        .reset_index()
    )

    # Newly added statistical column: Number of assemblies after merging
    df2["assembly_number"] = df2["assembly_id"].apply(len)

    # Save to Excel
    df2.to_excel(output_file, index=False)

    print(f"Final candidate data (genes) saved to: {output_file}")

def analysis_interpro(comparison_data_path, transformed_data_path, sequence_path, protein_path, high_threshold,
                                           basic_threshold, ratio_threhold, annotation_path, excel_file, final_output):

    # transform the data of clinker results
    transform_clinker_results(comparison_data_path, transformed_data_path)

    # extract sequences and prepare annotation using interpro
    candidate_dict = extract_all_sequences(transformed_data_path, sequence_path, protein_path, high_threshold,
                                           basic_threshold, ratio_threhold)

    interpro_all_sequences(protein_path, annotation_path)

    # read the dictionary of gene id
    id_dict = build_region_to_genes(excel_file)

    # load annotation results and generate output
    output_annotation_results(candidate_dict, annotation_path, id_dict, final_output)


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

    analysis_interpro(comparison_data_path, transformed_data_path, sequence_path, protein_path, high_threshold,
                      basic_threshold, ratio_threhold, annotation_path, excel_file, final_output)













