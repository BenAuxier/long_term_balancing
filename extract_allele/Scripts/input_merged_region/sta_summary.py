import os
import subprocess
import pandas as pd
from collections import Counter

def count_base_depth(base_depth, depth_hist_file):
    counter = Counter()

    with open(base_depth) as f:
        for line in f:
            depth = int(line.split()[2])
            counter[depth] += 1

    with open(depth_hist_file, "w") as out:
        out.write("depth\tcount\n")
        for d in sorted(counter):
            out.write(f"{d}\t{counter[d]}\n")

def prepare_statistic_data(base_depth, statistics_path, gene_region_depth, summary_out_file, region_output_file):

    # 1. copy the depth of each single nucleotide
    subprocess.run(["cp", base_depth, statistics_path])

    # 2. copy the depth of each gene
    #subprocess.run(["cp", gene_region_depth, statistics_path])
    depth_hist_file = f"{statistics_path}/depth_hist.tsv"
    count_base_depth(base_depth, depth_hist_file)

    # 3. copy the gene positions
    subprocess.run(["cp", summary_out_file, statistics_path])

    # 4. copy the region info and included gene info
    subprocess.run(["cp", region_output_file, statistics_path])



if __name__ == "__main__":
    reference_genome = "GCF_000002655.1"  # genome annotation should be GCF version
    species = "aspergillus_fumigatus"

    #reference_genome = "GCF_000184455.2"  # genome annotation should be GCF version (RefSeq)
    #species = "aspergillus_oryzae"  # the species

    main_path = f"/lustre/BIF/nobackup/leng010/test/{species}"
    results_path = f"{main_path}/results"

    base_depth = f"{main_path}/depth_calculation/tmp/single_nucleotides_depth.txt"
    gene_region_depth = f"{main_path}/depth_calculation/mean_depth_region.txt"

    # extract the length of chromosomes form the fai file
    ref_assembly = f"{main_path}/reference_genome/{reference_genome}_genomic.fna"
    ref_fai = f"{ref_assembly}.fai"

    summary_out_file = f"{main_path}/results/{species}_final_candidates.xlsx"

    region_output_file = f"{results_path}/genomic_region_genes.csv"

    statistics_path = f"{main_path}/results/statistics_data"
    os.makedirs(statistics_path, exist_ok=True)

    prepare_statistic_data(base_depth, statistics_path, gene_region_depth, summary_out_file, region_output_file)


