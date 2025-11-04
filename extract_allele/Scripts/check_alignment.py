from analyze_position import extract_align_seq

def check_assembly(bam_path, seq, start, end, full_cover=True):
    alignment = extract_align_seq(bam_path, seq, start, end, full_cover=True)
    seq = []
    for info in alignment['seq_info']:
        seq.append((info['Genome_accession'],info["seq_ID"]))
    return seq


bam_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/alignment/alignment_aspergillus_oryzae.sorted.bam"


seq = "NC_036440.1"
start = 1166960
end = 1168770
seq = check_assembly(bam_path, seq, start, end, full_cover=True)
print(len(seq))
