"""
Extract the allele of each gene in multiple genomes

"""
import subprocess
import csv
import pysam
import pandas as pd
from collections import defaultdict
import random
import os
import glob
from Bio import SeqIO

#from pysam.libchtslib import CFalse


####################################
# Generate covered and uncovered genomes at one position.
# Check sequences covered or uncovered on specific genomic region of BAM file
def extract_align_seq(bam_path, seq, start, end, full_cover = True):
    """
    Extract the sequence name that are aligned to one genomic position (chr, start, end) of the BAM file.
    The sequences are from different genome assemblies.
    :arg bam_path: Path to the BAM file.
    :arg seq: Sequence name (in reference genome) of the genomic position.
    :arg start: Start position.
    :arg end: End position.

    :return: align_info
    a dictionary that include list of the sequence names, each in a dictionary.
    Information including NCBI accession + assembly name + sequence name + position in that assembly.
    """

    align_info = {}
    # the position in the reference genome
    align_info["pos_info"] = {"seq": seq, "start": start, "end": end}
    # the list to store the assembly names
    align_info["seq_info"] = []

    bam = pysam.AlignmentFile(bam_path, "rb")

    for read in bam.fetch(seq, start, end):
        # Check if the read overlaps the region of interest
        # read.query_alignment_start / read.query_alignment_end are the aligned
        # positions on the read sequence
        # Coordinates are 0-based, left-inclusive and right-exclusive
        #if read.is_unmapped or read.is_secondary or read.is_supplementary:
            #continue

        # Only reads that cover the interval 100% are retained
        if full_cover:
            if read.reference_start > start or read.reference_end < end:
                continue

        # Create a dictionary to store read alignment information
        dic_align = {}
        name_all = read.query_name
        name_parse = name_all.split("_")
        ncbi_accession = name_parse[0] + "_" + name_parse[1]
        seq_id = name_parse[-1]

        if read.is_reverse == False:
            orientation = "forward"
            # the reads corridinates that align to the start-end of reference.
            read_corridinate_start = find_read_corridinate(read, start)
            read_corridinate_end = find_read_corridinate(read, end)

        else:
            orientation = "reverse"
            # the reads corridinates that align to the start-end of reference.
            read_corridinate_start = find_read_corridinate(read, end)
            read_corridinate_end = find_read_corridinate(read, start)

        # dic_align["name"] = name_all      # Complete read name
        dic_align["Genome_accession"] = ncbi_accession  # Read NCBI accession number
        dic_align["seq_ID"] = seq_id  # ID of the sequence, such as chromosome ID
        # if reverse then
        dic_align["start"] = read_corridinate_start  # Alignment start on read
        dic_align["end"] = read_corridinate_end  # Alignment end on read
        dic_align["orientation"] = orientation  # True if read aligns to reverse strand

        align_info["seq_info"].append(dic_align)
    #print(len(align_info["seq_info"]))
    bam.close()
    return align_info

def find_read_corridinate(read, position):
    """
    finding the corridinate of read on the position of the reference genome.

    :param read:
    :param position:
    :return:
    """
    read_self_start = read.query_alignment_start
    read_reference_start = read.reference_start  # aligned reads on reference genome

    read_cigar_list = read.cigartuples

    insertion = 0
    deletion = 0
    length_reference = 0
    hard_mask = 0
    length_ref_0 = position - read_reference_start + 1

    read_length_hard = 0
    for info in read_cigar_list:
        if info[0] == 5:  # hard mask
            read_length_hard += info[1]
    read_length = read.query_length + read_length_hard

    read_length_2 = 0
    for info in read_cigar_list:
        info = list(info)
        type_seq = str(info[0])
        length_seq = info[1]

        if type_seq == "2":  # deletion
            continue
        else:
            read_length_2 += length_seq

    for info in read_cigar_list:
        info = list(info)
        type_seq = str(info[0])
        length_seq = info[1]
        # print(type_seq,length_seq)
        if type_seq == "5":  # hard mask
            hard_mask = length_seq
        elif type_seq == "4":  # soft mask
            continue
        elif type_seq == "1":  # insertion
            insertion += length_seq
        elif type_seq == "2":  # deletion
            deletion += length_seq
            length_reference += length_seq
        elif type_seq in ["0", "6", "7", "8", "9"]:  # other situations, there are also type > 5, not considered
            length_reference += length_seq

        # what if 3, skipped region?

        if length_reference > length_ref_0:
            break

    read_corridinate = hard_mask + read_self_start + length_ref_0 + insertion - deletion

    if read.is_reverse:
        read_corridinate = read_length - read_corridinate + 1

    return (read_corridinate)


def extract_align_seq_from_dict(bam_path, position_info):
    """
    using function extract_align_seq() with a dictionary that extracted from annotation.
    :param bam_path:
    :param position_info:
    :return: align_info_B
    """
    seq = position_info["seq_name"]
    start = position_info["start"]
    end = position_info["end"]
    align_info_B = extract_align_seq(bam_path, seq, start, end)
    return align_info_B

def not_align_seq(align_info, assembly_path):
    """
    Find the genome assemblies that are not aligned to one genomic position.
    :param align_info: The output of extract_align_seq, with the list (align_info["seq_info"])
    consist of dictionaries of the mapped genome assemblies information at certain position.
    :param assembly_path: the path to a txt file, each row is a NCBI accession number of genome assembly.
    A list consisting of all genome assemblies that were used in the alignment.
    :return: not_align_list, A list including all uncovered genome assemblies at the genomic position.
    """
    align_genome = []
    for seq_info in align_info["seq_info"]:
        ID = seq_info["Genome_accession"]
        align_genome.append(ID)

    #load the list including all genome assemblies
    with open(assembly_path, "r", encoding="utf-8") as f:
        all_genome = [line.strip() for line in f]
    #print(all_genome)
    not_align_list = list(set(all_genome) - set(align_genome))
    position = align_info["pos_info"]
    print(f"Position on reference genome: sequence ID: {position['seq']}, "
          f"start position: {position['start']}, "
          f"end position: {position['end']}.\n"
          f"Overall {len(all_genome)} assemblies, aligned {len(align_genome)} assemblies, "
          f"not aligned {len(not_align_list)} assemblies.")
    return not_align_list
#######################################

# Process data from txt file, and filter it between up and down limits
def process_data(input_file):
    """
    process the data from a txt file.
    :param input_file: path to the analysis result (depth of bam file at mRNA positions).
    :return: rows: rows of the file (a list including dictionaries),
    fieldnames: keys/headers of the file
    """
    # Read and process data
    rows = []
    all_keys = set()

    with open(input_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")

            # Basic columns (first 8 columns + attributes + depth)
            base_cols = [
                "seqid", "source", "type", "start", "end",
                "score", "strand", "phase", "attributes", "depth"
            ]
            row = dict(zip(base_cols, parts))

            # Parse the 'attributes' column
            attributes = row["attributes"].split(";")
            for attr in attributes:
                attr = attr.strip()
                if not attr:
                    continue
                if "=" in attr:
                    key, value = attr.split("=", 1)
                else:
                    continue
                row[key] = value
                all_keys.add(key)

            rows.append(row)

        # Construct output column names: basic columns + all expanded attribute keys
        fieldnames = base_cols + sorted(all_keys)
        #print(fieldnames)

    return rows,fieldnames

def select_depth(rows,lower_limit,upper_limit):
    """
    Only keeps the rows with depth between lower_limit and upper_limit.
    :param rows:
    :param lower_limit:
    :param upper_limit:
    :return: filtered_rows, a list including dictionaries for each mRNA position
    """

    filtered_rows = []
    for row in rows:
        depth = row["depth"]
        depth = float(depth)
        if lower_limit < depth < upper_limit:
            # print(row)
            filtered_rows.append(row)
    return filtered_rows

def save_rows(output_file, filtered_rows):
    """
    Save the filtered rows to a csv file.
    :param output_file:
    :param filtered_rows:
    :return:
    """
    # Write to CSV file
    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in filtered_rows:
            writer.writerow(row)

    print("Completed! Output file:", output_file)

def extract_candidate_position(filtered_rows):
    """
    Extract the necessary position information from the filtered data.
    :param filtered_rows: A list including each rows as a dictionary.
    :return: candidate_data, a dictionary, including each mRNA candidates and its information as key-value pairs.
    """
    candidate_data = {}
    for row in filtered_rows:
        row_data = {
            "seqid": row["seqid"],
            "start": row["start"],
            "end": row["end"],
            "depth": row["depth"],
            "id": row["Name"],
            "locus_tag": row["locus_tag"]
        }
        candidate_data[row["Name"]] = row_data
    return candidate_data

####################################
#Select up and downstream positions
def read_gtf(gtf_path, keep_type=None):
    """
    Read a GTF/GFF file and organize it by seq_name.

    :param gtf_path: path to the GTF/GFF file
    :param keep_type: optional, only keep rows with this feature type (e.g., "exon")
    :return: dict, {seq_name: [row_dict, ...]} sorted by start
    """
    gtf_dict = defaultdict(list)

    with open(gtf_path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t", fieldnames=[
            "seq_name", "source", "type", "start", "end", "score", "strand", "phase", "attributes"
        ])

        for row in reader:
            if row["seq_name"].startswith("#"):
                continue  # skip header/comment lines
            if keep_type and row["type"] != keep_type:
                continue

            # parse attributes
            attr_dict = {}
            for attr in row["attributes"].split(";"):
                if attr.strip() == "":
                    continue
                key_value = attr.strip().split(" ", 1)
                if len(key_value) == 2:
                    key, value = key_value
                    attr_dict[key] = value.strip('"')

            row_data = {
                "seq_name": row["seq_name"],
                "type": row["type"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "gene_id": attr_dict.get("gene_id"),
                "transcript_id": attr_dict.get("transcript_id")
            }

            gtf_dict[row["seq_name"]].append(row_data)

    # sort each row with the order of start position
    for seq in gtf_dict:
        gtf_dict[seq].sort(key=lambda x: x["start"])

    return dict(gtf_dict)

def read_gff(gff_path, keep_type=None):
    """
    Read a GFF3 file and organize it by seq_name.

    :param gff_path: path to the GFF3 file
    :param keep_type: optional, only keep rows with this feature type (e.g., "exon")
    :return: dict, {seq_name: [row_dict, ...]} sorted by start
    """
    gff_dict = defaultdict(list)

    with open(gff_path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t", fieldnames=[
            "seq_name", "source", "type", "start", "end", "score", "strand", "phase", "attributes"
        ])

        for row in reader:
            if row["seq_name"].startswith("#"):
                continue  # skip header/comment lines
            if keep_type and row["type"] != keep_type:
                continue

            # parse attributes (GFF3 format: key=value;key=value;...)
            attr_dict = {}
            for attr in row["attributes"].split(";"):
                if attr.strip() == "":
                    continue
                if "=" in attr:
                    key, value = attr.strip().split("=", 1)
                    attr_dict[key] = value

            row_data = {
                "seq_name": row["seq_name"],
                "type": row["type"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "id": attr_dict.get("ID"),
                "parent": attr_dict.get("Parent"),
                "gene_id": attr_dict.get("gene_id"),         # 有些GFF也可能保留 GTF风格
                "transcript_id": attr_dict.get("transcript_id")
            }

            gff_dict[row["seq_name"]].append(row_data)
            #print(row_data)

    # sort each list by start position
    for seq in gff_dict:
        gff_dict[seq].sort(key=lambda x: x["start"])

    return dict(gff_dict)

def dict_rows_transfer(rows):
    dict_rows = {}
    for row in rows:
        name = row["Name"]
        dict_rows[name] = row
    return dict_rows

def find_position(transcript_id, sequence_id, annotation_sorted, up_num, down_num):
    """
    Retrieve the positional information of the m upstream and n downstream mRNAs for a specified mRNA.
    :param transcript_id: NCBI id, such as "XM_742446.1"
    :param sequence_id: sequence of the transcript, NCBI id, such as "NC_007194.1"
    :param annotation_sorted: the sorted annotation dictionary.
    :param n: number of positions to identify at the up and down-stream positions.
    :return: up_down_locations: a dictionary with the location of input mRNA, its up and down-stream mRNA positions.
    """
    annotation_sequences = annotation_sorted[sequence_id]
    annotation_length = len(annotation_sequences)
    for i in range(0, annotation_length):
        annotation = annotation_sequences[i]
        annotation_ID = annotation["transcript_id"]
        if annotation_ID == transcript_id:
            transcript_position = i
            position_info = annotation_sequences[transcript_position]

            number_upstream = transcript_position - 0
            number_downstream = annotation_length - transcript_position -1

            # Extract
            upstream_position = []
            downstream_position = []
            # extract upstream positions
            if number_upstream >= up_num:
                for j in range(transcript_position-up_num,i): # m upstream positions
                    upstream_position.append(annotation_sequences[j])
            elif up_num > number_upstream >= 1:
                for j in range(0,transcript_position): # m upstream positions
                    upstream_position.append(annotation_sequences[j])
            # extract downstream positions
            if number_downstream >= down_num:
                for j in range(transcript_position+1,transcript_position+down_num+1): # n downstream positions
                    downstream_position.append(annotation_sequences[j])
            elif down_num > number_downstream >= 1:
                for j in range(transcript_position+1,annotation_length): # n downstream positions
                    downstream_position.append(annotation_sequences[j])

            up_down_locations = {
                "position": position_info,
                "upstream_position": upstream_position,
                "downstream_position": downstream_position
            }
            return up_down_locations

        elif i == annotation_length-1:
            print(f"Warning! {transcript_id} is not found in {sequence_id} of reference genome.")

def find_position_depth(transcript_id, rows_dict):
    #print(transcript_id)
    pos_info = rows_dict[transcript_id]
    pos_depth = float(pos_info["depth"])
    return pos_depth

def filter_up_down_depth(up_down_locations, rows_dict, up_num, down_num, cutoff = 10):
    """

    :param up_down_locations:
    :param rows_dict:
    :param up_num:
    :param down_num:
    :param cutoff:
    :return:
    """
    upstream_position = up_down_locations["upstream_position"]
    downstream_position = up_down_locations["downstream_position"]
    sum_depth_up = 0
    sum_depth_down = 0
    for pos1 in upstream_position:
        #print(pos1)
        transcript_id1 = pos1["transcript_id"]
        sum_depth_up += find_position_depth(transcript_id1, rows_dict)
    for pos2 in downstream_position:
        transcript_id2 = pos2["transcript_id"]
        sum_depth_down += find_position_depth(transcript_id2, rows_dict)
    ave_depth_up = sum_depth_up/up_num
    ave_depth_down = sum_depth_down/down_num
    if ave_depth_up >= cutoff or ave_depth_down >= cutoff:
        return True
    else:
        return False



def find_candidate_align(up_down_locations,bam_path,assembly_path):
    """

    :param up_down_locations: a dictionary with the location of input mRNA, its up and down-stream mRNA positions.
    :return: up_down_alignment:
    position_info_align, alignment and un-alignment lists at the position.
    up_gene_align, alignment at upstream positions.
    down_gene_align, alignment at downstream positions
    """
    up_down_alignment = {}
    position_info = up_down_locations["position"] # one dictionary, analysis both align and unaligned genomes
    upstream_position = up_down_locations["upstream_position"] # list including dictionaries
    downstream_position = up_down_locations["downstream_position"] # list including dictionaries
    #print(position_info, upstream_position, downstream_position)

    position_involved_assembly = extract_align_seq_from_dict(bam_path, position_info)
    position_uninvolved_assembly = not_align_seq(position_involved_assembly, assembly_path)

    upstream_position_involved_assembly=[]
    for pos in upstream_position:
        position_involved_assembly = extract_align_seq_from_dict(bam_path, pos)
        upstream_position_involved_assembly.append(position_involved_assembly)

    downstream_position_involved_assembly=[]
    for pos in downstream_position:
        position_involved_assembly = extract_align_seq_from_dict(bam_path, pos)
        downstream_position_involved_assembly.append(position_involved_assembly)

    up_down_alignment = {
        "position_assembly": {
            "involved_assembly": position_involved_assembly,
            "uninvolved_assembly": position_uninvolved_assembly
        },
        "upstream_position_assembly": upstream_position_involved_assembly,
        "downstream_position_assembly": downstream_position_involved_assembly
    }
    #for key, value in up_down_alignment.items():
        #print(key,value)
    return up_down_alignment

# Function: randomly sample n unaligned sequences from the up_down_alignment of a given gene,
# and then analyze the state at each position for each of them
def random_select_assembly(up_down_alignment, assembly_num):
    """
    Randomly sample n unaligned sequences at the interested mRNA position.
    This function is currently not used!
    :param up_down_alignment:
    :param n:
    :return: selected_assemblies, list of selected assemblies
    """
    uninvolved_assembly_list = up_down_alignment["position_assembly"]["uninvolved_assembly"]
    number_uninvolved = len(uninvolved_assembly_list)
    if 0 < number_uninvolved < assembly_num:
        assembly_num = number_uninvolved
    selected_assemblies = random.sample(uninvolved_assembly_list, assembly_num)
    #print(selected_assemblies)
    return selected_assemblies

# Finding whether the selected assemblies are involved in the positions.
def find_candidate_involvement(up_down_alignment):
    """

    :param up_down_alignment:
    :param
    :return:
    """
    # get the random assembly list to test status
    selected_assemblies = up_down_alignment["position_assembly"]["uninvolved_assembly"]
    all_status = {}
    upstream_alignment = up_down_alignment["upstream_position_assembly"]
    downstream_alignment = up_down_alignment["downstream_position_assembly"]
    for assembly in selected_assemblies: #test each selected assembly in up and downstream positions
        assembly_status = {}
        upstream_status = []
        downstream_status = []
        for positions in upstream_alignment: # analyze each position
            seq_info_list = positions["seq_info"]
            genome_accessions = [d["Genome_accession"] for d in seq_info_list]
            if assembly in genome_accessions:
                upstream_status.append(True)
            else:
                upstream_status.append(False)
        assembly_status["upstream_status"] = upstream_status

        for positions in downstream_alignment: # analyze each position
            seq_info_list = positions["seq_info"]
            genome_accessions = [d["Genome_accession"] for d in seq_info_list]

            if assembly in genome_accessions:
                downstream_status.append(True)
            else:
                downstream_status.append(False)
        assembly_status["downstream_status"] = downstream_status
        all_status[assembly] = assembly_status
    #print(all_status)
    return all_status

def find_up_down_loci(all_status, up_down_locations,up_down_alignment): #
    """

    :param all_status:
    :param up_down_locations:
    :param up_down_alignment:
    :return: up_down_loci
    """
    up_down_loci = {}
    # the locus information related to this locus in th reference genome
    upstream_positions = up_down_locations["upstream_position"]
    downstream_positions = up_down_locations["downstream_position"]
    # the aligned genomes onto the locus
    upstream_align = up_down_alignment["upstream_position_assembly"]
    #print(upstream_align)
    downstream_align = up_down_alignment["downstream_position_assembly"]


    for genome_id, status in all_status.items():

        #print(genome_id,status)
        loci_selected = {}
        upstream_status = status["upstream_status"]
        downstream_status = status["downstream_status"]
        up_number = len(upstream_status)
        down_number = len(downstream_status)

        # incase situations without any findings
        loci_selected["upstream_rank"] = False
        loci_selected["upstream_ref_position"] = False
        loci_selected["upstream_read_sequence_ID"] = False
        loci_selected["upstream_read_position"] = False
        loci_selected["downstream_rank"] = False
        loci_selected["downstream_ref_position"] = False
        loci_selected["downstream_read_sequence_ID"] = False
        loci_selected["downstream_read_position"] = False


        for i in range(1,up_number+1): # the i th nearest loci in the upstream positions
            # find upstream true locus
            status_loci_up = upstream_status[-i]
            if status_loci_up: # if status_loci == True
                loci_selected["upstream_rank"] = i
                # The position information of the loci on the reference genome
                loci_selected["upstream_ref_position"] = upstream_positions[-i]
                ##
                # extract the information of the position in the selected assembly
                up_align_list = upstream_align[-i]["seq_info"]
                aligned_info_seq = []
                aligned_info_pos = {}
                for genome in up_align_list:
                    if genome["Genome_accession"] == genome_id:
                        # print(genome)
                        sequence_ID = genome["seq_ID"]
                        aligned_info_seq.append(sequence_ID)
                        # Assume that each sequence appears only once
                        aligned_info_pos[sequence_ID] = genome

                #print(aligned_info)
                loci_selected["upstream_read_sequence_ID"] = aligned_info_seq
                loci_selected["upstream_read_position"] = aligned_info_pos

                break  # exit the loop once condition is met

            elif i == up_number:
                loci_selected["upstream_rank"] = False
                loci_selected["upstream_ref_position"] = False
                loci_selected["upstream_read_sequence_ID"] = False
                loci_selected["upstream_read_position"] = False

        for j in range(1,down_number+1): # the i th nearest loci in the downstream positions
            # find downstream true locus
            status_loci_down = downstream_status[j-1]
            if status_loci_down: # if status_loci == True
                loci_selected["downstream_rank"] = j
                loci_selected["downstream_ref_position"] = downstream_positions[j-1]

                # extract the information of the position in the selected assembly
                down_align_list = downstream_align[j-1]["seq_info"]
                aligned_info_seq = []
                aligned_info_pos = {}
                for genome in down_align_list:
                    if genome["Genome_accession"] == genome_id:
                        #print(genome)
                        sequence_ID = genome["seq_ID"]
                        aligned_info_seq.append(sequence_ID)
                        # Assume that each sequence appears only once
                        aligned_info_pos[sequence_ID] = genome

                # print(aligned_info)
                # :25 error
                loci_selected["downstream_read_sequence_ID"] = aligned_info_seq
                #print("aligned_info_seq",aligned_info_seq)
                loci_selected["downstream_read_position"] = aligned_info_pos

                break  # exit the loop once condition is met

            elif j == down_number:
                loci_selected["downstream_rank"] = False
                loci_selected["downstream_ref_position"] = False
                loci_selected["downstream_read_sequence_ID"] = False
                loci_selected["downstream_read_position"] = False

        up_seq = loci_selected["upstream_read_sequence_ID"]
        down_seq = loci_selected["downstream_read_sequence_ID"]

        loci_selected["seq_chromosome"] = False
        if up_seq and down_seq:
            common_seq = set(up_seq).intersection(down_seq)

            if len(common_seq) > 0:
                loci_selected["seq_chromosome"] = common_seq

                # Assume there is only one shared seq in common_seq
                seq = list(common_seq)[0]

                # calculate interval between two reference genes in read genome
                read_up_position = loci_selected["upstream_read_position"][seq]
                read_down_position = loci_selected["downstream_read_position"][seq]
                orientation = read_up_position["orientation"]

                read_interval = abs(read_down_position["end"] - read_up_position["start"])
                """if orientation == "forward":
                    read_interval = read_down_position["end"] - read_up_position["start"]
                else:
                    read_interval = read_up_position["end"] - read_down_position["start"]"""

                loci_selected["read_interval"] = read_interval
                #print(read_up_position,'\t',read_down_position, interval)

                # calculate interval between two reference genes in reference genome
                ref_up_position = loci_selected["upstream_ref_position"]
                ref_down_position = loci_selected["downstream_ref_position"]
                ref_interval = ref_down_position["end"] - ref_up_position["start"]
                loci_selected["ref_interval"] = ref_interval

                interval_diff_pct = abs(ref_interval-read_interval)*100/ref_interval
                loci_selected["interval difference(%)"] = interval_diff_pct

        up_down_loci[genome_id] = loci_selected
        #print(loci_selected)

    return up_down_loci

def analyze_all_candidate_position(candidate_data,annotation_sorted,bam_path,assembly_path, up_num, down_num, lower_limit):
    """
    Find the position data for each of the candidate positions.
    :param candidate_data:
    :param annotation_sorted:
    :return: candidate_data_positions:
    """
    candidate_data_summary = {}
    for transcript_id,info in candidate_data.items():
        sequence_id = info["seqid"]
        # finding the up and downstream positions
        up_down_locations = find_position(transcript_id, sequence_id, annotation_sorted, up_num, down_num)

        # check the mean depth of the positions
        depth_status = filter_up_down_depth(up_down_locations, rows_dict, up_num, down_num, lower_limit)
        if depth_status == False:
            continue

        # finding the genes aligned to the positions
        up_down_alignment = find_candidate_align(up_down_locations,bam_path,assembly_path)
        # identify status of random genome assemblies at the positions
        all_status = find_candidate_involvement(up_down_alignment)
        # identify the nearest upstream and downstream loci
        up_down_loci = find_up_down_loci(all_status, up_down_locations,up_down_alignment)

        candidate_data_summary[transcript_id] = {
            "position_info": up_down_locations,
            "align_info": up_down_alignment,
            "status_info":all_status,
            "up_down_loci": up_down_loci
        }

    for key,value in candidate_data_summary.items():
        for key2, value2 in value["up_down_loci"].items():
            #print(key, key2, value2)
            continue

    return candidate_data_summary

def extract_allele_sequence(assembly_dir, candidate_gene, genome_accession, contig, start, end, orientation, type_seq, output_path):
    """
    Extract allele sequence for a candidate gene from genome assembly.

    Parameters:
        candidate_gene (str): gene ID (e.g., "XM_743603.1")
        genome_accession (str): genome accession (e.g., "GCA_051225625.1")
        contig (str): contig name (e.g., "CP097565.1")
        start (int): start position
        end (int): end position
        output_path (str): output folder path

    Returns:
        str: path of saved fasta file, or warning message if failed
    """
    # 0. finding genome assembly
    assembly_pattern = os.path.join(assembly_dir, f"{genome_accession}*.fna")
    matches = glob.glob(assembly_pattern)

    if not matches:
        warning = f"Warning: No genome assembly found for {genome_accession}"
        print(warning)
        #return warning
    genome_assembly_path = matches[0]

    # 1. scale
    seq_start = max(1, start - 5000)
    seq_end = end + 5000

    # 2. extract sequence
    target_seq = None
    for record in SeqIO.parse(genome_assembly_path, "fasta"):
        #print(record.id)
        if contig in record.id:
            #print("find")
            target_seq = record.seq[seq_start - 1:seq_end]  # Biopython 序列是0-based
            break

    if target_seq is None:
        warning = f"Warning: Contig {contig} not found in {genome_accession}"
        print(warning)
        return warning

    if orientation.lower() == "reverse":
        target_seq = target_seq.reverse_complement()

    # 3. output path
    output_dir = os.path.join(output_path, candidate_gene)
    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, f"{candidate_gene}_{type_seq}_{genome_accession}_{contig}-{seq_start}-{seq_end}.fa")

    # 4. save results
    with open(output_file, "w") as f:
        f.write(f">{genome_accession}_{contig}:{seq_start}-{seq_end}\n")
        f.write(str(target_seq) + "\n")

    print(f"Sequence saved to {output_file}")
    return output_file

def find_allele_sequence_inbetween(bam_path, reference_genome, assembly_dir,candidate_data_summary,output_path, assembly_num):
    """

    :param reference_genome:
    :param assembly_dir:
    :param candidate_data_summary:
    :param output_path:
    :return:
    """
    for gene,summary in candidate_data_summary.items():

        seq_list = []
        m = 0
        for genome, info in summary["up_down_loci"].items():
            seq = info["seq_chromosome"]

            if seq == False:
                continue

            if info["interval difference(%)"] > 200:
                continue

            seq_list.append(seq) # check the number of chromosomes

            seq_info = list(seq)[0] # assume there is only one seq

            upstream_rank = info["upstream_rank"]
            downstream_rank = info["downstream_rank"]
            upstream_read_position = info["upstream_read_position"][seq_info]
            downstream_read_position = info["downstream_read_position"][seq_info]

            start_read = min(upstream_read_position["start"], downstream_read_position["start"])
            end_read = max(upstream_read_position["end"], downstream_read_position["end"])
            orientation = upstream_read_position["orientation"]

            if end_read - start_read > 100000:
                continue

            extract_allele_sequence(
                assembly_dir,
                gene,
                genome,
                seq_info,
                start_read,
                end_read,
                orientation,
                "allele",
                output_path
            )
            m += 1
            if m == assembly_num:
                break

        if seq_list != []:
            # also extract the reference sequence
            up_down_alignment = summary["align_info"]
            seq_info_ref = up_down_alignment["position_assembly"]["involved_assembly"]["pos_info"]["seq"]
            # for a, b in up_down_alignment.items():
            # print(a,b)

            if up_down_alignment["upstream_position_assembly"] == [] or up_down_alignment[
                "downstream_position_assembly"] == []:
                continue

            start_ref = up_down_alignment["upstream_position_assembly"][0]["pos_info"]["start"]
            end_ref = up_down_alignment["downstream_position_assembly"][-1]["pos_info"]["end"]
            extract_allele_sequence(
                assembly_dir,
                gene,
                reference_genome,
                seq_info_ref,
                start_ref,
                end_ref,
                "forward",
                "reference",
                output_path
            )

            #extract other allele sequence similar to the reference
            align_reference_allele = extract_align_seq(bam_path, seq_info_ref, start_ref, end_ref, full_cover=True)
            seq_reference_alleles = align_reference_allele["seq_info"]

            n = 0
            for seq_reference_allele in seq_reference_alleles:
                genome_accession_reference_allele = seq_reference_allele["Genome_accession"]
                seq_ID_reference_allele = seq_reference_allele["seq_ID"]
                start_reference_allele = seq_reference_allele["start"]
                end_reference_allele = seq_reference_allele["end"]
                orientation_reference_allele =seq_reference_allele["orientation"]
                extract_allele_sequence(
                    assembly_dir,
                    gene,
                    genome_accession_reference_allele,
                    seq_ID_reference_allele,
                    start_reference_allele,
                    end_reference_allele,
                    orientation_reference_allele,
                    "reference_alleles",
                    output_path
                )
                n += 1
                if n == assembly_num:
                    break


#file paths
# The modified csv file that contains the mRNA with depth of interests.
depth_path = r"/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/check_coverage/test3/all51_to_GCF_000002655.1_meandepth.txt"
# Reference genome annotation of the BAM file
gtf_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/reference/GCF_000002655.1_ASM265v1_genomic.gtf"
gff_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/reference/GCF_000002655.1_ASM265v1_genomic.gff"

# Bam file, used multiple genome assemblies align to the reference genome
bam_path = ("/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/multi_align_2/40_assembly/normal_align/"
            "all51_to_GCF_000002655.1_asm5.sorted.bam")
# The txt file that includes the genome assemblies used in alignment.
assembly_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/multi_align_2/40_assembly/accessions.txt"
# path to the 53 assemblies
assembly_dir = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/multi_align_2/40_assembly/genome_assemblies"
# output base path of the fasta files
output_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/check_coverage/test4_extract_seq/extract_allele"


# settings
up_num = 5
down_num = 5
assembly_num = 7
lower_limit = 10
upper_limit = 40
reference_genome = "GCF_000002655.1"

# processes the input candidate mRNAs
# Input and output file paths
rows,fieldnames = process_data(depth_path)
print("rows",len(rows))
filtered_rows = select_depth(rows,lower_limit,upper_limit) # the candidate mRNA data
print("filtered_rows",len(filtered_rows))
rows_dict = dict_rows_transfer(rows)


candidate_data = extract_candidate_position(filtered_rows)

# Select up and downstream n th mRNA position of the candidate position
#annotation_sorted = read_gtf(annotation_path, keep_type="transcript")
annotation_sorted = read_gff(gff_path, keep_type="mRNA")

"""#extract the up dan downstream positions
up_down_locations = find_position("XM_749896.2", "NC_007196.1", annotation_sorted, up_num, down_num)
#print(up_down_locations)
up_down_alignment = find_candidate_align(up_down_locations,bam_path,assembly_path)
selected_assemblies = random_select_assembly(up_down_alignment)
all_status = find_candidate_involvement(up_down_alignment)
up_down_loci = find_up_down_loci(all_status, up_down_locations,up_down_alignment)
filter_up_down_depth(up_down_locations, rows_dict, up_num, down_num, lower_limit)"""

#test the main code
candidate_data_test = dict(list(candidate_data.items())[0:])
#print(candidate_data_test)
candidate_data_summary = analyze_all_candidate_position(candidate_data_test,
            annotation_sorted,bam_path,assembly_path, up_num, down_num, lower_limit)

gene_between = find_allele_sequence_inbetween(bam_path, reference_genome, assembly_dir,candidate_data_summary,output_path, assembly_num)


"""#test
output_path = "/lustre/BIF/nobackup/leng010/test/Asp_fumigatus/check_coverage/test_seq_5/MAT1-2-4_test2"
candidate_data_test = {'XM_749897.2':
                      {'seqid': 'NC_007196.1',
                       'start': '1523105',
                       'end': '1524006',
                       'depth': '35.0000000',
                       'id': 'XM_749897.2',
                       'locus_tag': 'AFUA_3G06160'}}

annotation_sorted = read_gff(gff_path, keep_type="mRNA")
candidate_data_summary = analyze_all_candidate_position(candidate_data_test,
            annotation_sorted,bam_path,assembly_path, up_num, down_num, lower_limit)

gene_between = find_allele_sequence_inbetween(bam_path, reference_genome, assembly_dir,candidate_data_summary,output_path, assembly_num)
"""

