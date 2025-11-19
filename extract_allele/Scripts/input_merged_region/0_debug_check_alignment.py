from analyze_position import extract_align_seq

def check_assembly(bam_path, seq, start, end, full_cover=True):
    alignment = extract_align_seq(bam_path, seq, start, end, full_cover=True)
    seq = []
    for info in alignment['seq_info']:
        seq.append((info['Genome_accession'],info["seq_ID"]))
        #print(info)
    return seq


bam_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/alignment/alignment_aspergillus_oryzae.sorted.bam"


seq = "NC_036442.1"
start = 305024
end = 305657
sequence = check_assembly(bam_path, seq, start, end, full_cover=True)


up_start = 301000
up_end = 302000
up_seq = check_assembly(bam_path, seq, up_start, up_end, full_cover=True)

down_start = 315000
down_end = 316000
down_seq = check_assembly(bam_path, seq, down_start, down_end, full_cover=True)

sequences = set(up_seq) & set(down_seq) & set(sequence)
print(sequences)

