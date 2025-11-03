import pandas as pd
import subprocess


def calculate_genome_coverage():
    """
    Run bedtools genomecov to calculate genome coverage depth
    """
    # Define command as a list for subprocess
    cmd = [
        "bedtools", "genomecov",
        "-ibam", "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/alignment/alignment_aspergillus_fumigatus.sorted.bam",
        "-d"
    ]

    # Open output file for writing
    with open("/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/depth_calculation/genome_depth.txt", "w") as outfile:
        # Run the command and redirect stdout to the file
        subprocess.run(cmd, stdout=outfile, check=True)

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


# 示例调用
if __name__ == "__main__":
    calculate_average_depth(
        depth_file="/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/depth_calculation/genome_depth.txt",
        gff_file="/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/genome_assemblies/reference_genome/GCF_000002655.1_genomic_mRNA.gff",
        output_file="/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/depth_calculation/GCF_000002655.1_genomic_mRNA_with_depth.gff"
    )




