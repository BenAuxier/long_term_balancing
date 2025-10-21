import csv
import pysam
import pandas as pd
from collections import defaultdict
import random
import os
import glob
from Bio import SeqIO

# 1. process depth result data
def process_data(input_file):
    """
    Process the depth data from a txt file (input_file).
    :param input_file: path to the analysis result (depth of alignment at each position).
    :return: rows: all rows of the txt file. Example:
    [{'seq_ID': 'NC_007194.1',
     'type': 'mRNA',
     'start': 216,
     'end': 836,
     'strand': '+',
     'depth': '50.0000000',
     'Parent': 'gene-AFUA_1G00100',
     'Name': 'XM_001481640.1',
     'locus_tag': 'AFUA_1G00100',
     'product': 'putative MFS monocarboxylate transporter',
     'transcript_id': 'XM_001481640.1'}, ...]
     """

    # Read and process data
    rows = []
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
            keys_to_keep = ['seq_ID', 'type', 'start', 'end', 'strand', 'depth', 'Parent', 'Name',
                            'locus_tag', 'product', 'transcript_id']
            row_new = {}
            for key in keys_to_keep:
                row_new[key] = row[key]

            rows.append(row_new)
    return rows


def select_depth(rows, lower_limit, upper_limit):
    """
    Only keeps the rows with depth between lower_limit and upper_limit.
    :param rows:
    :param lower_limit:
    :param upper_limit:
    :return: filtered_rows, same as rows, a list including
    dictionaries for each mRNA position
    """

    filtered_rows = []
    for row in rows:
        depth = row["depth"]
        depth = float(depth)
        if lower_limit < depth < upper_limit:
            # print(row)
            filtered_rows.append(row)

    return filtered_rows


def dict_rows_transfer(rows):
    """
    Transfer a list into dictionary
    :param rows: a list including dictionaries for each mRNA position
    :return: a dictionary including dictionaries for each mRNA position
    """
    rows_dict = {}
    for row in rows:
        name = row["Name"]
        rows_dict[name] = row
    return rows_dict


def extract_candidate_position_list(filtered_rows):
    """
    Transfer the list to a dictionary, only keep the necessary information
    :param filtered_rows: A list including each rows as a dictionary.
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
    for row in filtered_rows:
        if not row["seq_ID"] in candidate_data_seq.keys():
            candidate_data_seq[row["seq_ID"]] = []

        row_data = {
            "seq_ID": row["seq_ID"],
            "start": row["start"],
            "end": row["end"],
            "depth": row["depth"],
            "id": row["Name"],
            "locus_tag": row["locus_tag"]
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
    {seq_ID: {transcript_id: annotation, ...}}
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
        mixed_data = {}

        for i in range(0, len(candidates)):
            candidate_i = candidates[i]
            gene_id_i = candidate_i["id"]
            annotation_candidate_i = annotation_seq[gene_id_i]

            if mixed_data == {}:  # if nothing in mixed data, select the current annotation

                # if there is not a mixed data, import this annotation as the mixed data
                mixed_data = {
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

            # Then there is already a mixed data from i-n to i
            start_i = mixed_data["start"]
            rank_i = mixed_data["rank"]

            # if there is only one candidate in the genomic sequence
            if len(candidates) == 1:
                candidate_merge[seq].append(mixed_data)
                break

            if i == len(candidates) - 1:
                candidate_merge[seq].append(mixed_data)
                break

            elif i <= len(candidates) - 2:
                # candidate i + 1
                candidate_i_1 = candidates[i + 1]
                id_i_1 = candidate_i_1["id"]
                annotation_candidate_i_1 = annotation_seq[id_i_1]
                rank_i_1 = annotation_candidate_i_1["rank"]

                # compare rank of i+1 to i
                if abs(rank_i_1 - rank_i) <= 2:  # mix
                    mixed_data_new = {
                        "seq_ID": mixed_data["seq_ID"],
                        "region_name": f"{mixed_data["start_gene"]}-{id_i_1}",
                        "start_gene": mixed_data["start_gene"],
                        "end_gene": id_i_1,
                        "start": start_i,
                        "end": annotation_candidate_i_1["end"],
                        "rank": rank_i_1,
                        "gene_included": mixed_data["gene_included"]
                    }
                    mixed_data_new["gene_included"].append(id_i_1)
                    mixed_data_new["gene_number"] = len(mixed_data_new["gene_included"])

                    mixed_data = mixed_data_new

                else:
                    candidate_merge[seq].append(mixed_data)
                    # print(mixed_data["region_name"])
                    mixed_data = {}

    return candidate_merge











