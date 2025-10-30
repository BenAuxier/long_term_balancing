from analyze_position import extract_align_seq

def check_assembly(bam_path, seq, XM_001481378_start, XM_001481378_end, full_cover=True):
    alignment = extract_align_seq(bam_path, seq, XM_001481378_start, XM_001481378_end, full_cover=True)
    seq = []
    for info in alignment['seq_info']:
        seq.append((info['Genome_accession'],info["seq_ID"]))
    return seq


bam_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/alignment/alignment_aspergillus_fumigatus.sorted.bam"
seq = "NC_007200.1"
XM_001481378_start = 1178886
XM_001481378_end = 1179356

upstream_start = 1178200
upstream_end = 1178400

downstream_start = 1179800
downstream_end = 1180000

"""XM_001481378_seq = check_assembly(bam_path, seq, XM_001481378_start, XM_001481378_end, full_cover=True)

upstream_seq = check_assembly(bam_path, seq, upstream_start, upstream_end, full_cover=True)
downstream_seq = check_assembly(bam_path, seq, downstream_start, downstream_end, full_cover=True)

up_down_intersection = list(set(upstream_seq) & set(downstream_seq))
diff = list(set(up_down_intersection) - set(XM_001481378_seq))
#print(diff)"""


seq = "NC_007194.1"
start = 4350000
end = 4360000
seq = check_assembly(bam_path, seq, start, end, full_cover=True)
print(len(seq))
