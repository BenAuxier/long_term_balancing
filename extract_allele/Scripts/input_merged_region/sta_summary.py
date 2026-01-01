import os
import subprocess
import pandas as pd
from collections import Counter

import csv

import csv




def whitespace_txt_to_csv(txt_file, csv_file, num_lines):
    """
    Convert a whitespace-delimited txt file to CSV and
    add depth_percentage = mean_depth / num_lines.

    Parameters
    ----------
    txt_file : str
        Input whitespace-delimited text file
    csv_file : str
        Output CSV file
    num_lines : int
        Total number of assemblies
    header : list or None
        Column names (WITHOUT depth_percentage); if None, no header is written
    """
    header = ["seqid", "source", "feature", "start", "end", "score",
              "strand", "phase", "id", "mean_depth", "covered_bases"]

    with open(txt_file) as fin, open(csv_file, "w", newline="") as fout:
        writer = csv.writer(fout)
        writer.writerow(header + ["depth_percentage"])

        for line in fin:
            line = line.strip()
            if not line:
                continue

            fields = line.split()

            # mean_depth is the 10th column (index 9)
            mean_depth = float(fields[9])
            depth_percentage = mean_depth * 100 / num_lines

            # gene length
            writer.writerow(fields + [depth_percentage])


def count_base_depth(base_depth, depth_hist_file, num_lines):
    counter = Counter()

    with open(base_depth, "r", encoding="utf-8") as f:
        lines = f.readlines()
        num_count = len(lines)

    with open(base_depth) as f:
        for line in f:
            depth = int(line.split()[2])
            counter[depth] += 1

    with open(depth_hist_file, "w") as out:
        out.write("depth,depth_percentage,count,count_percentage\n")
        for d in sorted(counter):
            normalized_depth = d * 100 / num_lines
            normalized_count = counter[d] * 100 / num_count
            out.write(f"{d},{normalized_depth},{counter[d]},{normalized_count}\n")

    print(num_lines, num_count)


def prepare_statistic_data(base_depth, statistics_path, gene_bs_region_depth, all_gene_depth,
                           summary_out_file, region_output_file, assembly_list, refseq_candidate_file):

    # 0. calculate number of assemblies
    with open(assembly_list, "r", encoding="utf-8") as f:
        lines = f.readlines()
        num_lines = len(lines)

    # 1. copy the depth of each single nucleotide
    #subprocess.run(["cp", base_depth, statistics_path])
    depth_hist_file = f"{statistics_path}/depth_hist.csv"
    count_base_depth(base_depth, depth_hist_file, num_lines)

    # 2. copy the depth of balancing selected genes (only the genes overlap with balancing selection regions)
    #subprocess.run(["cp", all_gene_depth, statistics_path])
    gene_depth_file = f"{statistics_path}/filtered_gene_depth.csv"
    whitespace_txt_to_csv(gene_bs_region_depth, gene_depth_file, num_lines)

    # 3. copy the depth of each gene
    all_depth_file = f"{statistics_path}/all_gene_depth.csv"
    whitespace_txt_to_csv(all_gene_depth, all_depth_file, num_lines)

    # 3. copy the gene positions
    subprocess.run(["cp", summary_out_file, statistics_path])

    # 4. copy the region info and included gene info
    subprocess.run(["cp", region_output_file, statistics_path])

    # 5. copy candidate gene list
    subprocess.run(["cp", refseq_candidate_file, statistics_path])


if __name__ == "__main__":
    reference_genome = "GCF_000002655.1"  # genome annotation should be GCF version
    species = "aspergillus_fumigatus"

    #reference_genome = "GCF_000184455.2"  # genome annotation should be GCF version (RefSeq)
    #species = "aspergillus_oryzae"  # the species

    base_path = "/lustre/BIF/nobackup/leng010/test"
    main_path = f"/lustre/BIF/nobackup/leng010/test/{species}"
    results_path = f"{main_path}/results"

    base_depth = f"{main_path}/depth_calculation/tmp/single_nucleotides_depth.txt"
    gene_bs_region_depth = f"{main_path}/depth_calculation/mean_depth_region.txt"
    all_gene_depth = f"{main_path}/depth_calculation/mean_depth_gene.txt"

    # extract the length of chromosomes form the fai file
    ref_assembly = f"{main_path}/reference_genome/{reference_genome}_genomic.fna"
    ref_fai = f"{ref_assembly}.fai"

    summary_out_file = f"{main_path}/results/{species}_final_candidates.xlsx"

    region_output_file = f"{results_path}/genomic_region_genes.csv"

    assembly_list = f"{base_path}/genome_accessions/{species}.txt"

    statistics_path = f"{main_path}/results/statistics_data"
    os.makedirs(statistics_path, exist_ok=True)

    go_path = f"{main_path}/go_analysis"
    refseq_candidate_file = f"{go_path}/all_candidate_genes.txt"

    prepare_statistic_data(base_depth, statistics_path, gene_bs_region_depth, all_gene_depth, summary_out_file, region_output_file, assembly_list, refseq_candidate_file)


