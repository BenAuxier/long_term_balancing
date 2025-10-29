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

    return feature_counts


# Example usage:
if __name__ == "__main__":
    gff_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/genome_assemblies/reference_genome/GCA_000184455.3_genomic.gff"
    result = count_gff_features(gff_path)

    # Print results
    for feature, count in result.items():
        print(f"{feature}\t{count}")
