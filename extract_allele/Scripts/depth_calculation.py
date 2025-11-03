import pandas as pd
import subprocess
import os

def calculate_genome_coverage(bam_path,main_path):
    """
    Run bedtools genomecov to calculate genome coverage depth
    """
    output_path = f"{main_path}/depth_calculation"
    # Ensure output directory exists
    os.makedirs(output_path, exist_ok=True)

    depth_single_nucleotides = f"{output_path}/single_nucleotides_depth.txt"

    # Define command as a list for subprocess
    cmd = [
        "bedtools", "genomecov",
        "-ibam", bam_path,
        "-d"
    ]

    # Open output file for writing
    with open(depth_single_nucleotides, "w") as outfile:
        # Run the command and redirect stdout to the file
        subprocess.run(cmd, stdout=outfile, check=True)

    print(f"depth for each nucleotide are saved to {depth_single_nucleotides}")
    return depth_single_nucleotides

def calculate_average_depth(depth_file, gff_file, output_file):
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
    for _, row in gff_df.iterrows():
        seqid = row["seqid"]
        start, end = row["start"], row["end"]

        # Extract intervals from depth data corresponding to chromosomes
        if seqid in depth_groups:
            df = depth_groups[seqid]
            region_depths = df[(df["pos"] >= start) & (df["pos"] <= end)]
            if len(region_depths) > 0:
                avg_depth = region_depths["depth"].mean()
            else:
                avg_depth = 0
        else:
            avg_depth = 0

        avg_depths.append(avg_depth)

    gff_df["average_depth"] = avg_depths

    # 保存结果
    gff_df.to_csv(output_file, sep="\t", header=False, index=False)

    print(f"✅ Completed: Results have been saved to {output_file}")
    return output_file

def calculate_depth_all(bam_path, main_path, gff_file):
    depth_single_nucleotides = calculate_genome_coverage(bam_path, main_path)
    depth_path = f"{main_path}/depth_calculation/mean_depth.txt"
    calculate_average_depth(
        depth_single_nucleotides,
        gff_file,
        depth_path
    )
    return depth_path


# 示例调用
if __name__ == "__main__":
    base_path = "/lustre/BIF/nobackup/leng010/test"
    species = "aspergillus_fumigatus"
    reference_genome = "GCF_000002655.1"

    main_path = f"{base_path}/{species}"
    bam_path = f"{main_path}/alignment/alignment_aspergillus_fumigatus.sorted.bam"
    gff_file = f"{main_path}/genome_assemblies/reference_genome/GCF_000002655.1_genomic_mRNA.gff"


    depth_path = calculate_depth_all(bam_path, main_path, gff_file)



