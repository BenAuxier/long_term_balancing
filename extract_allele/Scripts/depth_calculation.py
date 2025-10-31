#!/usr/bin/env python3
import sys
from collections import defaultdict

def load_depth(depth_file):
    """读取 samtools depth 输出"""
    depth_dict = defaultdict(dict)
    with open(depth_file) as f:
        for line in f:
            chrom, pos, depth = line.strip().split('\t')
            depth_dict[chrom][int(pos)] = int(depth)
    return depth_dict


def calc_gff_depth(gff_file, depth_dict, output_file):
    """根据 GFF 区间计算平均 depth"""
    with open(gff_file) as gff, open(output_file, "w") as out:
        for line in gff:
            if line.startswith("#") or not line.strip():
                out.write(line)
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                out.write(line)
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])

            region_depths = [
                depth_dict[chrom][pos]
                for pos in range(start, end + 1)
                if pos in depth_dict[chrom]
            ]

            avg_depth = sum(region_depths) / len(region_depths) if region_depths else 0.0
            fields.append(f"AvgDepth={avg_depth:.2f}")
            out.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    depth_file = "depth.txt"        # samtools depth 输出
    gff_file = "input.gff"
    output_file = "output_with_depth.gff"

    depth_dict = load_depth(depth_file)
    calc_gff_depth(gff_file, depth_dict, output_file)
