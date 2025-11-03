#!/usr/bin/env python3
import re

from prepare_alignment import run_bedtools_depth
# information
reference_genome = "GCF_000002655.1" # genome annotation should be GCF version
species = "aspergillus_fumigatus"
augustus_species = "aspergillus_fumigatus"
type_annotation = "mRNA" # type of annotation used in depth calculation, the third column
ID_label = "transcript_id" # this is the key that the gene/mRNA id follows
key_words = None

##########################################################################
# file paths, including all input files
base_path = "/lustre/BIF/nobackup/leng010/test"

#assembly_list, this file need to create manually
assembly_list = f"{base_path}/genome_accessions/{species}.txt"
##########################################################################
# path to specific species
main_path = f"{base_path}/{species}"
assembly_dir= f"{main_path}/genome_assemblies"
ref_assembly= f"{assembly_dir}/reference_genome/{reference_genome}_genomic.fna"
ref_gff= f"{assembly_dir}/reference_genome/{reference_genome}_genomic.gff"
gff_filtered= f"{assembly_dir}/reference_genome/{reference_genome}_genomic_{type_annotation}.gff"
bam_path = f"{main_path}/alignment/alignment_{species}.sorted.bam"



def gff_to_bed(gff_file, bed_file, feature_type=None):
    """
    将 GFF/GTF 文件转换为 BED6 格式文件。

    参数:
        gff_file: str, 输入的 GFF/GTF 文件路径
        bed_file: str, 输出的 BED 文件路径
        feature_type: str 或 None, 若指定，则仅转换该类型的 feature (如 'gene', 'exon')
    """
    with open(gff_file) as fin, open(bed_file, "w") as fout:
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            chrom, source, ftype, start, end, score, strand, phase, attributes = fields

            # 若指定类型，则只保留该类型
            if feature_type and ftype != feature_type:
                continue

            # 坐标转换: GFF 为 1-based, BED 为 0-based 且右端不包含
            bed_start = int(start) - 1
            bed_end = int(end)

            # 提取 name（优先 gene_id 或 ID 或 Name）
            name_match = re.search(r'locus_tag=([^;]+)', attributes)
            name_match = name_match or re.search(r'gene_id "([^"]+)"', attributes)
            name_match = name_match or re.search(r'Name=([^;]+)', attributes)
            name = name_match.group(1) if name_match else ftype

            # 若 score 无效，则替换为 0
            if score == ".":
                score = "0"

            fout.write(f"{chrom}\t{bed_start}\t{bed_end}\t{name}\t{score}\t{strand}\n")

    print(f"✅ 已生成 BED 文件: {bed_file}")


# 示例用法
if __name__ == "__main__":
    gff_to_bed(gff_filtered,
               "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/depth_calculation/output.bed",
               feature_type="mRNA")
