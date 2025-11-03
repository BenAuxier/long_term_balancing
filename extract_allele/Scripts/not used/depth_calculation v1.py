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


#depth_path = run_bedtools_depth(gff_filtered, main_path, species, reference_genome)
#!/usr/bin/env python3
import pysam

def calc_avg_depth(bam_path, gff_path, output_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    out = open(output_path, "w")

    with open(gff_path) as gff:
        for line in gff:
            if line.startswith("#") or not line.strip():
                out.write(line)
                continue

            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue

            chrom = fields[0]
            start = int(fields[3]) - 1  # GFF 是 1-based, pysam 是 0-based
            end = int(fields[4])

            # 使用 pileup 统计 depth，不包括 deletion/skip
            total_depth = 0
            total_bases = 0

            for pileupcolumn in bam.pileup(
                chrom, start, end,
                truncate=True,              # 仅限区间内
                stepper="nofilter",          # 不跳过任何reads
                ignore_overlaps=False,
                ignore_orphans=False,
                min_base_quality=0
            ):
                pos = pileupcolumn.reference_pos
                if pos < start or pos >= end:
                    continue

                # 只统计非 deletion/skip 的 reads
                n = sum(1 for pileupread in pileupcolumn.pileups
                        if not pileupread.is_del and not pileupread.is_refskip)

                total_depth += n
                total_bases += 1

            avg_depth = total_depth / total_bases if total_bases > 0 else 0.0

            # 添加到 GFF 末尾
            fields.append(f"AvgDepth={avg_depth:.2f}")
            out.write("\t".join(fields) + "\n")

    out.close()
    bam.close()


calc_avg_depth(
        bam_path,
        gff_filtered,
        output_path="/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/depth_calculation/output_with_depth.gff"
    )














