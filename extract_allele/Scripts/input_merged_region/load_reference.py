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



def read_gff_augustus(gff_path, ID_label = "locus_tag", keep_type=None):
    """
        Read an augustus GFF file and organize it by seq_ID.

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

            row_data["id"] = row["attributes"] # whatever gene or mRNA, Augustus annotates only one id in attributes

            for key, value in row.items():
                if key == "start" or key == "end":
                    row_data[key] = int(value)
                    continue
                row_data[key] = value

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

def create_ID_dictionary(gff_path, output_csv, ID_label = "locus_tag", keep_type = "CDS"):
    """
    Create ID dictionary between XP_ and other such as XM_ ID of CDS.
    Used to transfer the ID to XP, to label them in clinker
    :param gff_path:
    :return:
    """
    ID_list = []

    ID_annotation = read_gff(gff_path, ID_label, keep_type)

    for seq, annotation_list in ID_annotation.items():
        for annotation in annotation_list:
            try:
                protein_id = annotation["protein_id"]
            except KeyError:
                continue

            gene_id = annotation[ID_label]
            mrna_id = annotation["Parent"][4:].split('-')[-1]

            ID_list.append([protein_id, gene_id, mrna_id])

    header = ["protein_id", "gene_id", "mrna_id"]
    # save as csv
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)  # 写入表头
        writer.writerows(ID_list)  # 写入数据

    return True

def csv_to_dict(file_path, key_header, value_header):
    """
    Convert two columns from a CSV file into a dictionary.
    :param file_path: CSV file path
    :param key_header: optional: "protein_id", "gene_id", "mrna_id"
    :param value_header: optional: "protein_id", "gene_id", "mrna_id"
    :return: dict
    """
    result = {}
    with open(file_path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = row[key_header]
            value = row[value_header]
            result[key] = value

    return result

def load_annotation_refseq(gff_path, ID_label, type_annotation):
    annotation_sorted = read_gff(gff_path, ID_label, type_annotation)
    annotation_sorted_dict = read_gff_dict(annotation_sorted)
    return annotation_sorted, annotation_sorted_dict

def load_annotation_augustus(gff_path_augustus, ID_label, type_annotation):
    annotation_sorted_augustus = read_gff_augustus(gff_path_augustus, ID_label, type_annotation)
    annotation_sorted_dict_augustus = read_gff_dict(annotation_sorted_augustus)
    return annotation_sorted_augustus, annotation_sorted_dict_augustus


if __name__ == "__main__":
    gff_path = "/lustre/BIF/nobackup/leng010/test/magnaporthe_grisea/genome_assemblies/reference_genome/GCF_000002495.2_genomic.gff"
    ID_label = "transcript_id"
    CDS_dict = create_ID_dictionary(gff_path,ID_label)
    print(CDS_dict)
