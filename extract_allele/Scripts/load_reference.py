import csv
from collections import defaultdict

def count_gff_features(gff_file):
    """
    Count the occurrences of each feature type (3rd column) in a GFF file.

    Parameters
    ----------
    gff_file : str
        Path to the input GFF file.

    Returns
    -------
    dict
        A dictionary where keys are feature types (3rd column)
        and values are their counts.
    """
    feature_counts = {}

    with open(gff_file, 'r', encoding='utf-8') as file:
        for line in file:
            # Skip comments or empty lines
            if line.startswith('#') or not line.strip():
                continue

            # Split the line into columns by tab
            columns = line.strip().split('\t')

            # Ensure the line has at least 3 columns
            if len(columns) < 3:
                continue

            # Extract the feature type (3rd column)
            feature_type = columns[2]

            # Count occurrences
            feature_counts[feature_type] = feature_counts.get(feature_type, 0) + 1

    for feature, count in feature_counts.items():
        print(f"{feature}\t{count}")

    return feature_counts

# load annotation
def read_gff(gff_path, ID_label = "locus_tag", keep_type=None):
    """
    Read a GFF3 file and organize it by seq_ID.

    :param gff_path: path to the GFF3 file
    :param keep_type: optional, only keep candidate_data with this feature type (e.g., "exon")
    :return: a dictionary, {seq_ID: [row_dict, ...]} sorted by start position
    """
    gff_dict = defaultdict(list)

    with open(gff_path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t", fieldnames=[
            "seq_ID", "source", "type", "start", "end", "score", "strand", "phase", "attributes"
        ])

        for row in reader:
            if row["seq_ID"].startswith("#"):
                continue  # skip header/comment lines
            if keep_type and row["type"] != keep_type:
                continue

            # parse attributes (GFF3 format: key=value;key=value;...)
            attr_dict = {}
            for attr in row["attributes"].split(";"):
                if attr.strip() == "":
                    continue
                if "=" in attr:
                    key, value = attr.strip().split("=", 1)
                    attr_dict[key] = value

            row_data = {
                "seq_ID": row["seq_ID"],
                "source": row.get("source", "."),
                "type": row.get("type", "."),
                "start": int(row["start"]),
                "end": int(row["end"]),
                "score": row.get("score", "."),
                "strand": row.get("strand", "."),
                "phase": row.get("phase", "."),
                "id": attr_dict.get(ID_label, "."),
                "locus_tag": attr_dict.get("locus_tag", "."),
                #"transcript_id": attr_dict.get("transcript_id", "."),
                "attributes": row.get("attributes", ".")
            }

            gff_dict[row["seq_ID"]].append(row_data)

    # sort each list by start position
    for seq in gff_dict:
        gff_dict[seq].sort(key=lambda x: x["start"])

    annotation_sorted = annotation_rank(dict(gff_dict))

    return annotation_sorted


def annotation_rank(annotation_sorted):
    for chromosome, annotation_info in annotation_sorted.items():
        i = 1
        for annotation in annotation_info:
            annotation["rank"] = i
            i += 1
    return annotation_sorted


def read_gff_dict(annotation_sorted):
    """
    Convert annotation storage format.
    :param annotation_sorted: a dictionary, {seq_ID: [row_dict, ...]} sorted by start position
    :return: dict, {seq_ID: {id: annotation, ...}}
    """
    annotation_dict = {}
    for seq, annotation_list in annotation_sorted.items():
        annotation_dict[seq] = {}
        for annotation in annotation_list:
            annotation_id = annotation["id"]
            annotation_dict[seq][annotation_id] = annotation

    return annotation_dict

def create_ID_dictionary(gff_path,ID_label_1,ID_label_2,keep_type=None):
    #all_annotation = read_gff(gff_path, ID_label_1, "CDS")
    gff_dict = defaultdict(list)

    with open(gff_path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t", fieldnames=[
            "seq_ID", "source", "type", "start", "end", "score", "strand", "phase", "attributes"
        ])

        for row in reader:
            if row["seq_ID"].startswith("#"):
                continue  # skip header/comment lines
            if keep_type and row["type"] != keep_type:
                continue

            # parse attributes (GFF3 format: key=value;key=value;...)
            attr_dict = {}
            for attr in row["attributes"].split(";"):
                if attr.strip() == "":
                    continue
                if "=" in attr:
                    key, value = attr.strip().split("=", 1)
                    row[key] = value


def load_annotation(gff_path, ID_label, type_annotation):
    annotation_sorted = read_gff(gff_path, ID_label, type_annotation)
    annotation_sorted_dict = read_gff_dict(annotation_sorted)
    return annotation_sorted, annotation_sorted_dict

if __name__ == "__main__":
    gff_path = ""
    ID_label_1 = ""
    ID_label_2 = ""
