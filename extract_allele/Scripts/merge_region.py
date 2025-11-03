# 1. process depth result data
def process_data(input_file,ID_label = "locus_tag"):
    """
    Process the depth data from a txt file (input_file).
    :param input_file: path to the analysis result (depth of alignment at each position).
    :return: candidate_data: all candidate_data of the txt file. Example:
    [{'seq_ID': 'NC_007194.1',
     'type': 'mRNA',
     'start': 216,
     'end': 836,
     'strand': '+',
     'depth': '50.0000000',
     'Parent': 'gene-AFUA_1G00100',
     'locus_tag': 'AFUA_1G00100',
     'product': 'putative MFS monocarboxylate transporter',
     'transcript_id': 'XM_001481640.1'}, ...]
     """

    # Read and process data
    candidate_data = []
    all_keys = set()

    with open(input_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")

            # Basic columns (first 8 columns + attributes + depth)
            base_cols = [
                "seq_ID", "source", "type", "start", "end",
                "score", "strand", "phase", "attributes", "depth"
            ]
            row = dict(zip(base_cols, parts))
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])

            # Parse the 'attributes' column
            attributes = row["attributes"].split(";")
            for attr in attributes:
                attr = attr.strip()
                if not attr:
                    continue
                if "=" in attr:
                    key, value = attr.split("=", 1)
                else:
                    continue
                row[key] = value
                all_keys.add(key)

            # select column
            keys_to_keep = ['seq_ID', 'type', 'start', 'end', 'strand', 'depth','locus_tag', 'product']
            row_new = {}
            for key in keys_to_keep:
                row_new[key] = row.get(key, ".")

            row_new["id"] = row.get(ID_label, ".")

            candidate_data.append(row_new)
    return candidate_data


def select_depth(candidate_data, lower_limit, upper_limit):
    """
    Only keeps the candidate_data with depth between lower_limit and upper_limit.
    :param candidate_data:
    :param lower_limit:
    :param upper_limit:
    :return: filtered_candidate_data, same as candidate_data, a list including
    dictionaries for each locus position
    """

    filtered_candidate_data = []
    for row in candidate_data:
        depth = row["depth"]
        depth = float(depth)
        if lower_limit < depth < upper_limit:
            filtered_candidate_data.append(row)

    return filtered_candidate_data


def dict_candidate_data_transfer(candidate_data):
    """
    Transfer a list into dictionary
    :param candidate_data: a list including dictionaries for each mRNA position
    :return: a dictionary including dictionaries for each mRNA position
    """
    candidate_data_dict = {}
    for candidate in candidate_data:
        id = candidate["id"]
        candidate_data_dict[id] = candidate
    return candidate_data_dict


def extract_candidate_position_list(filtered_candidate_data):
    """
    Transfer the list to a dictionary, only keep the necessary information
    :param filtered_candidate_data: A list including each candidate_data as a dictionary.
    :return: candidate_data_seq, a dictionary, with key:value =
    sequence_id : [candiate_info1, candiate_info2, ...]

    {'NC_007194.1':
    [{'seq_ID': 'NC_007194.1',
    'start': 29598,
    'end': 29981,
    'depth': '12.0000000',
    'id': 'XM_744675.1',
    'locus_tag': 'AFUA_1G00160'}, ...]
    xxxx: ......
    }

    """
    candidate_data_seq = {}
    for row in filtered_candidate_data:
        if not row["seq_ID"] in candidate_data_seq.keys():
            candidate_data_seq[row["seq_ID"]] = []

        row_data = {
            "seq_ID": row["seq_ID"],
            "start": row["start"],
            "end": row["end"],
            "depth": row["depth"],
            #"locus_tag": row["locus_tag"]
            "id": row["id"]
        }
        candidate_data_seq[row["seq_ID"]].append(row_data)

    return candidate_data_seq

def merge_candidate_position(candidate_data_seq, annotation_dict):
    """
    Merge adjacent or only one-position-separated candidate regions into a
    larger candidate region.
    :param candidate_data_seq:
    {'NC_007194.1':
    [{'seq_ID': 'NC_007194.1',
    'start': 29598,
    'end': 29981,
    'depth': '12.0000000',
    'id': 'XM_744675.1',
    'locus_tag': 'AFUA_1G00160'}, ...]
    xxxx: ......
    }
    :param annotation_dict: dict, genome annotation
    {seq_ID: {id: annotation, ...}}
    :return: candidate_merge, dict
    {
    NC_007197.1:
    [{'seq_ID': 'NC_007197.1',
    'region_name': 'XM_741280.1-XM_741278.1',
    'start_gene': 'XM_741280.1',
    'end_gene': 'XM_741278.1',
    'start': 223151,
    'end': 228584,
    'rank': 72,
    'gene_included': ['XM_741280.1', 'XM_741279.1', 'XM_741278.1'],
    'gene_number': 3
    }, ...],
    seq_ID_xxx: xxx,
    ......
    }
    """
    candidate_merge = {}

    for seq, candidates in candidate_data_seq.items():
        annotation_seq = annotation_dict[seq]
        candidate_merge[seq] = []
        merged_data = {}

        for i in range(0, len(candidates)):
            candidate_i = candidates[i]
            gene_id_i = candidate_i["id"]
            annotation_candidate_i = annotation_seq[gene_id_i]

            if merged_data == {}:  # if nothing in merged data, select the current annotation

                # if there is not a merged data, import this annotation as the merged data
                merged_data = {
                    "seq_ID": seq,
                    "region_name": gene_id_i,
                    "start_gene": gene_id_i,
                    "end_gene": gene_id_i,
                    "start": annotation_candidate_i["start"],
                    "end": annotation_candidate_i["end"],
                    "rank": annotation_candidate_i["rank"],
                    "gene_included": [gene_id_i],
                    "gene_number": 1
                }

            # Then there is already a merged data from i-n to i
            start_i = merged_data["start"]
            rank_i = merged_data["rank"]

            # if there is only one candidate in the genomic sequence
            if len(candidates) == 1:
                candidate_merge[seq].append(merged_data)
                break

            if i == len(candidates) - 1:
                candidate_merge[seq].append(merged_data)
                break

            elif i <= len(candidates) - 2:
                # candidate i + 1
                candidate_i_1 = candidates[i + 1]
                id_i_1 = candidate_i_1["id"]
                #print(candidates)
                annotation_candidate_i_1 = annotation_seq[id_i_1]
                rank_i_1 = annotation_candidate_i_1["rank"]

                # compare rank of i+1 to i
                if abs(rank_i_1 - rank_i) <= 2:  # mix
                    merged_data_new = {
                        "seq_ID": merged_data["seq_ID"],
                        "region_name": f"{merged_data["start_gene"]}-{id_i_1}",
                        "start_gene": merged_data["start_gene"],
                        "end_gene": id_i_1,
                        "start": start_i,
                        "end": annotation_candidate_i_1["end"],
                        "rank": rank_i_1,
                        "gene_included": merged_data["gene_included"]
                    }
                    merged_data_new["gene_included"].append(id_i_1)
                    merged_data_new["gene_number"] = len(merged_data_new["gene_included"])

                    merged_data = merged_data_new

                else:
                    candidate_merge[seq].append(merged_data)
                    merged_data = {}

    return candidate_merge

def process_results(depth_path,lower_limit, upper_limit,annotation_sorted_dict, ID_label):
    candidate_data = process_data(depth_path, ID_label)
    filtered_candidate_data = select_depth(candidate_data, lower_limit, upper_limit)  # the candidate mRNA data
    print(lower_limit, upper_limit)
    candidate_data_seq = extract_candidate_position_list(filtered_candidate_data)
    candidate_merge = merge_candidate_position(candidate_data_seq, annotation_sorted_dict)
    return candidate_merge









