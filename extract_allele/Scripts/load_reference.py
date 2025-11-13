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
        Read a GFF file and organize it by seq_ID.

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

            row_data = {}

            # parse attributes (GFF3 format: key=value;key=value;...)
            attr_dict = {}
            for attr in row["attributes"].split(";"):
                if attr.strip() == "":
                    continue
                if "=" in attr:
                    key, value = attr.strip().split("=", 1)
                    row[key] = value

            for key, value in row.items():
                if key == "attributes":
                    continue
                if key == "start" or key == "end":
                    row_data[key] = int(value)
                    continue
                row_data[key] = value

            row_data["attributes"] = row["attributes"]
            row_data["id"] = row_data.get(ID_label, ".")

            gff_dict[row["seq_ID"]].append(row_data)

    # sort each list by start position
    for seq in gff_dict:
        gff_dict[seq].sort(key=lambda x: x["start"])

    annotation_sorted = annotation_rank(dict(gff_dict))

    return annotation_sorted


def read_gff_old(gff_path, ID_label = "locus_tag", keep_type=None):
    """
    not used version!
    Read a GFF file and organize it by seq_ID.

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

def create_ID_dictionary(gff_path, ID_label = "locus_tag"):
    """
    Create ID dictionary between XP_ and other such as XM_ ID of CDS.
    Used to transfer the ID to XP, to label them in clinker
    :param gff_path:
    :return:
    """
    CDS_dict = {}

    CDS_annotation = read_gff(gff_path, ID_label, "CDS")

    for seq, annotation_list in CDS_annotation.items():
        for annotation in annotation_list:
            try:
                protein_id = annotation["protein_id"]
            except KeyError:
                continue

            if ID_label == "locus_tag":
                using_ID = annotation["locus_tag"]
            elif ID_label == "transcript_id":
                using_ID = annotation["Parent"][4:].split('-')[-1]

            CDS_dict[using_ID] = protein_id

    return CDS_dict

def load_annotation(gff_path, ID_label, type_annotation):
    annotation_sorted = read_gff(gff_path, ID_label, type_annotation)
    annotation_sorted_dict = read_gff_dict(annotation_sorted)
    return annotation_sorted, annotation_sorted_dict

if __name__ == "__main__":
    gff_path = "/lustre/BIF/nobackup/leng010/test/magnaporthe_grisea/genome_assemblies/reference_genome/GCF_000002495.2_genomic.gff"
    ID_label = "transcript_id"
    CDS_dict = create_ID_dictionary(gff_path,ID_label)
    print(CDS_dict)
