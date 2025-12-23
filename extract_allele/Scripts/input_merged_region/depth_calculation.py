import pandas as pd
import subprocess
import os

def calculate_base_depth(bam_path,main_path):
    """
    Run bedtools genomecov to calculate genome depth
    """
    output_path = f"{main_path}/depth_calculation/tmp"
    # Ensure output directory exists
    os.makedirs(output_path, exist_ok=True)

    base_depth = f"{output_path}/single_nucleotides_depth.txt"

    # Define command as a list for subprocess
    cmd = [
        "bedtools", "genomecov",
        "-ibam", bam_path,
        "-d"
    ]

    # Open output file for writing
    with open(base_depth, "w") as outfile:
        # Run the command and redirect stdout to the file
        subprocess.run(cmd, stdout=outfile, check=True)

    print(f"depth for each nucleotide are saved to {base_depth}")
    return base_depth

def calculate_average_depth(depth_file, gff_file, output_file):
    """
    calculate the average depth of for each region in gff_file using selected nucleotides in depth_file
    :param depth_file:
    :param gff_file:
    :param output_file:
    :return:
    """
    # Read Depth Files
    depth_cols = ["seqid", "pos", "depth"]
    depth_df = pd.read_csv(depth_file, sep="\t", names=depth_cols)

    # To accelerate query speed, grouped by chromosome
    depth_groups = {seqid: df for seqid, df in depth_df.groupby("seqid")}

    # Read GFF files (skipping comment lines)
    gff_cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gff_df = pd.read_csv(
        gff_file, sep="\t", comment="#", names=gff_cols, dtype={"start": int, "end": int}
    )

    avg_depths = []
    bases = []
    gene_length = []
    for _, row in gff_df.iterrows():
        seqid = row["seqid"]
        start, end = row["start"], row["end"]

        # Extract intervals from depth data corresponding to chromosomes
        avg_depth = 0
        included_bases = 0

        if seqid in depth_groups:
            df = depth_groups[seqid]
            region_depths = df[(df["pos"] >= start) & (df["pos"] <= end)]
            included_bases = region_depths["pos"].nunique()

            if included_bases > 0:
                avg_depth = region_depths["depth"].mean()
            #else:
                #avg_depth = 0
        #else:
            #avg_depth = 0
        avg_depths.append(avg_depth)
        bases.append(included_bases)

    gff_df["average_depth"] = avg_depths
    gff_df["included_bases"] = bases

    # save result
    gff_df.to_csv(output_file, sep="\t", header=False, index=False)

    print(f"Completed: Results have been saved to {output_file}")
    return output_file


def filter_depth_to_bed(base_depth, filtered_nucleotides, lower_limit, upper_limit):
    """
    Only select the position with depth between lower_limit and upper_limit.
    input_file: base_depth, f"{main_path}/depth_calculation"
    output_file: filtered_depth_nucleotides

    """

    with open(base_depth, "r") as infile, open(filtered_nucleotides, "w") as outfile:
        for line in infile:
            chrom, pos, depth = line.strip().split("\t")
            pos = int(pos)
            depth = float(depth)

            if lower_limit <= depth <= upper_limit:
                outfile.write(f"{chrom}\t{pos - 1}\t{pos}\n")
    print(f"Result have been saved to {base_depth}")

    return base_depth

def merge_nucleotides(filtered_nucleotides, merged_genomic_region, interval = "1"):
    """

    :param filtered_nucleotides:
    :param merged_genomic_region:
    :param interval:
    :return:
    """
    cmd = [
        "bedtools", "merge", "-d", interval,
        "-i", filtered_nucleotides
    ]

    # Open output file for writing
    with open(merged_genomic_region, "w") as outfile:
        # Run the command and redirect stdout to the file
        subprocess.run(cmd, stdout=outfile, check=True)
    print(f"Result have been saved to {merged_genomic_region}")

    return merged_genomic_region

def filter_region_length(merged_genomic_region, filtered_genomic_region, minimal_length = 100):
    """

    :param merged_genomic_region:
    :param filtered_genomic_region:
    :param minimal_length:
    :return:
    """
    with open(merged_genomic_region, "r") as infile, open(filtered_genomic_region, "w") as outfile:
        for line in infile:
            cols = line.strip().split()
            if len(cols) < 3:
                continue  # skip

            seqid, start, end = cols[0], int(cols[1]), int(cols[2])
            region_length = end - start - 1

            line_new = f"{seqid}\t{start}\t{end}\t{region_length}\n"

            if region_length > minimal_length:
                outfile.write(line_new)

    print(f"Result have been saved to {filtered_genomic_region}")

    return filtered_genomic_region

def filter_region_nucleotides_not_used(base_depth, filtered_genomic_region,filtered_region_nucleotides):
    """
    currently not used.
    :param base_depth:
    :param filtered_genomic_region:
    :param filtered_region_nucleotides:
    :return:
    """
    dict_single_nucleotides = {}
    with open(base_depth, "r") as infile:
        for line in infile:
            nucleotide_cols = line.strip().split()
            # transfer 0 based bed format to 1-based
            n_seqid, n_position, n_depth = nucleotide_cols[0], int(nucleotide_cols[1]), int(nucleotide_cols[2])
            nucleotide_info = {"seq_id": n_seqid, "position": n_position, "depth": n_depth, "info": line}
            if n_seqid not in dict_single_nucleotides.keys():
                dict_single_nucleotides[n_seqid] = []
            dict_single_nucleotides[n_seqid].append(nucleotide_info)
    #print(dict_single_nucleotides)

    with open(filtered_genomic_region, "r") as region_file:
        for region in region_file:
            region_cols = region.strip().split()
            # r_start is 0-based, end is 1-based
            r_seqid, r_start, r_end = region_cols[0], int(region_cols[1]), int(region_cols[2])
            seq_nucleotides = dict_single_nucleotides[r_seqid]
            with open(filtered_region_nucleotides, "a") as outfile:
                for i in range(r_start, r_end):
                    ouput_line = seq_nucleotides[i]["info"]
                    outfile.write(ouput_line)

    print(f"Result have been saved to {filtered_region_nucleotides}")
    return filtered_region_nucleotides

def filter_region_nucleotides(depth_file, region_file, output_file):
    regions = {}

    # 1. read BED region（0-based, [start, end)）
    with open(region_file) as f:
        for line in f:
            cols = line.strip().split()
            if len(cols) < 3:
                continue
            seqid, start, end = cols[0], int(cols[1]), int(cols[2])
            regions.setdefault(seqid, []).append((start, end))

    # 2. sort region with start
    for seqid in regions:
        regions[seqid].sort()

    # 3. filter depth_file
    with open(depth_file) as infile, open(output_file, "w") as out:
        for line in infile:
            seqid, pos, depth = line.split()
            pos = int(pos)  # 1-based

            if seqid not in regions:
                continue

            # BED [start, end) vs depth 1-based
            for start, end in regions[seqid]:
                if start < pos <= end:
                    out.write(line)
                    break   # each line only write one time

    print(f"Temporary result of filtered bases are saved to {output_file}")

def create_filter_region_nucleotides(tmp_path, lower_limit, upper_limit, interval, minimal_length):
    """
    Filtering the bases in genomic regions with depth below lower limit and above upper limit.


    :param main_path:
    :param lower_limit:
    :param upper_limit:
    :param interval:
    :param minimal_length:
    :return:
    """
    base_depth = f"{tmp_path}/single_nucleotides_depth.txt"
    filtered_nucleotides = f"{tmp_path}/filtered_nucleotides.bed"
    merged_genomic_region = f"{tmp_path}/merged_genomic_region.bed"
    filtered_genomic_region = f"{tmp_path}/filtered_genomic_region.bed"
    filtered_region_nucleotides = f"{tmp_path}/filtered_region_nucleotides.txt"

    # first select the bases with depth between lower limit and upper limit
    filter_depth_to_bed(base_depth, filtered_nucleotides, lower_limit, upper_limit)

    # merge the bases to generate genomic region
    merge_nucleotides(filtered_nucleotides, merged_genomic_region, interval)

    # filter the regions based on their length
    filter_region_length(merged_genomic_region, filtered_genomic_region, minimal_length)

    #filter the depth for bases file, only keep those bases overlap with the merged balancing selection regions
    filter_region_nucleotides(base_depth, filtered_genomic_region, filtered_region_nucleotides)

    return filtered_region_nucleotides


def calculate_depth_all(gene_depth, gene_region_depth, bam_path, main_path, gff_file, lower_limit, upper_limit, base_interval, minimal_length):
    """

    :param gene_depth:
    :param gene_region_depth:
    :param bam_path:
    :param main_path:
    :param gff_file:
    :param lower_limit:
    :param upper_limit:
    :param base_interval:
    :param minimal_length:
    :return:
    """
    # first calculate the depth of each base in the reference genome
    base_depth = calculate_base_depth(bam_path, main_path)

    # calculate average depth for each gene in the gff file.
    # this is only used for reference gene selection in analyze_all_candidate_position()
    calculate_average_depth(base_depth, gff_file, gene_depth)

    # select bases in balancing selection region
    tmp_path = f"{main_path}/depth_calculation/tmp"
    filtered_region_nucleotides = create_filter_region_nucleotides(tmp_path, lower_limit, upper_limit, base_interval, minimal_length)

    # calculate the average depth of the balancing selection genomic region in each gene
    calculate_average_depth(filtered_region_nucleotides, gff_file, gene_region_depth)

    return True



if __name__ == "__main__":
    base_path = "/lustre/BIF/nobackup/leng010/test"
    reference_genome = "GCF_000002655.1"  # genome annotation should be GCF version
    species = "aspergillus_fumigatus"

    main_path = f"{base_path}/{species}"
    bam_file = f"{main_path}/alignment/alignment_{species}.sorted.bam"
    ref_gff_augustus = f"{main_path}/reference_genome/{reference_genome}_genomic_AUGUSTUS.gff"
    gff_augustus_filtered = f"{main_path}/reference_genome/{reference_genome}_genomic_AUGUSTUS_transcript.gff"

    base_depth = f"{main_path}/depth_calculation/tmp/single_nucleotides_depth.txt"
    filtered_nucleotides = f"{main_path}/depth_calculation/tmp/filtered_nucleotides.bed"
    merged_genomic_region = f"{main_path}/depth_calculation/tmp/merged_genomic_region.bed"
    filtered_genomic_region = f"{main_path}/depth_calculation/tmp/filtered_genomic_region.bed"
    filtered_region_nucleotides = f"{main_path}/depth_calculation/tmp/filtered_region_nucleotides.txt"
    lower_limit = 4
    upper_limit = 6
    minimal_length = 100
    base_interval = "2"

    filter_region_nucleotides(base_depth, filtered_genomic_region, filtered_region_nucleotides)

    gene_region_depth = f"{main_path}/depth_calculation/mean_depth_region.txt"
    # calculate the average depth of the balancing selection genomic region in each gene
    calculate_average_depth(filtered_region_nucleotides, gff_augustus_filtered, gene_region_depth)