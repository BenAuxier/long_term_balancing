import csv
import pysam
from collections import defaultdict
import random
from merge_region import dict_gene_data


# Generate covered and uncovered genomes at one position.
# Check sequences covered or uncovered on specific genomic region of BAM file
def extract_align_seq(bam_path, seq, start, end, full_cover=True):
    """
    Extract the sequence name that are aligned to one genomic position (chr, start, end) of the BAM file.
    The sequences are from divergent genome assemblies.
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
        # if read.is_unmapped or read.is_secondary or read.is_supplementary:
        # continue

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
            # the reads coordinate that align to the start-end of reference.
            read_coordinate_start = find_read_coordinate(read, start)
            read_coordinate_end = find_read_coordinate(read, end)

        else:
            orientation = "reverse"
            # the reads coordinates that align to the start-end of reference.
            read_coordinate_start = find_read_coordinate(read, end)
            read_coordinate_end = find_read_coordinate(read, start)

        # if the read exhibits deletion between start and end
        read_length = abs(read_coordinate_end-read_coordinate_start+1)
        ref_lenth = abs(end-start+1)
        #print(read_length,ref_lenth)

        if read_length < 0.9 * ref_lenth: # if the read have larger than 1/10 deletion
            continue

        # dic_align["name"] = name_all      # Complete read name
        dic_align["Genome_accession"] = ncbi_accession  # Read NCBI accession number
        dic_align["seq_ID"] = seq_id  # ID of the sequence, such as chromosome ID
        # if reverse then
        dic_align["start"] = read_coordinate_start  # Alignment start on read
        dic_align["end"] = read_coordinate_end  # Alignment end on read
        dic_align["orientation"] = orientation  # True if read aligns to reverse strand

        align_info["seq_info"].append(dic_align)
    # print(len(align_info["seq_info"]))
    bam.close()
    return align_info


def find_read_coordinate(read, position):
    """
    finding the coordinate of read on the position of the reference genome.

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
            if type_seq == "2":  # deletion
                deletion = deletion - (length_reference-length_ref_0)
            break

    read_coordinate = hard_mask + read_self_start + length_ref_0 + insertion - deletion

    if read.is_reverse:
        read_coordinate = read_length - read_coordinate + 1

    return (read_coordinate)


def extract_align_seq_from_dict(bam_path, position_info):
    """
    using function extract_align_seq() with a dictionary that extracted from annotation.
    :param bam_path:
    :param position_info:
    :return: align_info_B
    """
    seq = position_info["seq_ID"]
    start = position_info["start"]
    end = position_info["end"]
    align_info_B = extract_align_seq(bam_path, seq, start, end)
    return align_info_B


def find_aligned_assembly(align_info):
    aligned_genome = []
    for seq_info in align_info["seq_info"]:
        seq_ID = seq_info["Genome_accession"]
        aligned_genome.append(seq_ID)

    return aligned_genome


def find_not_aligned_assembly(aligned_genome, align_info, assembly_path):
    """
    Find the genome assemblies that are not aligned to one genomic position.
    :param align_info: The output of extract_align_seq, with the list (align_info["seq_info"])
    consist of dictionaries of the mapped genome assemblies information at certain position.
    :param assembly_path: the path to a txt file, each row is a NCBI accession number of genome assembly.
    A list consisting of all genome assemblies that were used in the alignment.
    :return: not_align_list, A list including all uncovered genome assemblies at the genomic position.
    """

    # load the list including all genome assemblies
    with open(assembly_path, "r", encoding="utf-8") as f:
        all_genome = [line.strip() for line in f]
    # print(all_genome)
    not_align_list = list(set(all_genome) - set(aligned_genome))
    position = align_info["pos_info"]
    # print(f"Position on reference genome: sequence ID: {position['seq']}, "
    # f"start position: {position['start']}, "
    # f"end position: {position['end']}.\n"
    # f"Overall {len(all_genome)} assemblies, aligned {len(aligned_genome)} assemblies, "
    # f"not aligned {len(not_align_list)} assemblies.")
    return not_align_list

####################################
# data analysis
def find_position_seq(sequence_id, start, end, annotation_sorted, up_num, down_num, type_annotation):
    """
    Retrieve the positional information of the m upstream and n downstream mRNAs for a specified mRNA.
    :param id: NCBI id, such as "XM_742446.1"
    :param sequence_id: sequence of the NCBI id, such as "NC_007194.1"
    :param annotation_sorted: the sorted annotation dictionary.
    :param n: number of positions to identify at the up and down-stream positions.
    :return: up_down_locations: a dictionary with the location of input mRNA, its up and down-stream mRNA positions.
    {
        'position': {'seq_ID': 'NC_007195.1', 'start': 4650883, 'end': 4653737},
        'type': 'mRNA',
        'upstream_position': [
            {'seq_ID': 'NC_007195.1', 'type': 'mRNA', 'start': 4643751, 'end': 4644257, 'id': 'rna-XM_750980.1',
             'parent': 'gene-AFUA_2G17380', 'gene_id': None, 'transcript_id': 'XM_750980.1', 'rank': 1555}, ...],
        'downstream_position': [
            {'seq_ID': 'NC_007195.1', 'type': 'mRNA', 'start': 4653899, 'end': 4654978, 'id': 'rna-XM_750985.1',
             'parent': 'gene-AFUA_2G17430', 'gene_id': None, 'transcript_id': 'XM_750985.1', 'rank': 1560}, ...]
    }
    """
    annotation_sequences = annotation_sorted[sequence_id]
    annotation_length = len(annotation_sequences)
    up_position = 0
    down_position = 0

    for i in range(1, annotation_length):
        if annotation_sequences[i - 1]["start"] <= start <= annotation_sequences[i]["start"]:
            up_position = i
        if annotation_sequences[i - 1]["end"] <= end <= annotation_sequences[i]["end"]:
            down_position = i - 1
        if i == annotation_length - 1:
            break

    if up_position > 0 and down_position > 0:
        # overall number of positions in up and down stream positions
        number_upstream = up_position
        number_downstream = annotation_length - down_position - 1

        # Extract
        upstream_position = []
        downstream_position = []
        # extract upstream positions
        if number_upstream >= up_num:
            for j in range(up_position - up_num, up_position):  # m upstream positions
                upstream_position.append(annotation_sequences[j])
        elif up_num > number_upstream >= 1:
            for j in range(0, up_position):  # m upstream positions
                upstream_position.append(annotation_sequences[j])

        # extract downstream positions
        if number_downstream >= down_num:
            for j in range(down_position + 1, down_position + down_num + 1):  # n downstream positions
                downstream_position.append(annotation_sequences[j])
        elif down_num > number_downstream >= 1:
            for j in range(down_position + 1, annotation_length):  # n downstream positions
                downstream_position.append(annotation_sequences[j])
        position_info = {
            'seq_ID': sequence_id,
            'start': start,
            'end': end
        }

        up_down_locations_seq = {
            "position": position_info,
            'type': type_annotation,
            "upstream_position": upstream_position,
            "downstream_position": downstream_position
        }
        return up_down_locations_seq

    else:
        print(
            f"Warning! sequence {sequence_id},start {start}, end {end} is not found in {sequence_id} of reference genome.")


def find_position_depth(id, candidate_data_dict):
    pos_info = candidate_data_dict[id]
    pos_depth = float(pos_info["depth"])
    return pos_depth


def filter_up_down_depth(up_down_locations, candidate_data_dict, up_num, down_num, cutoff):
    """

    :param up_down_locations:
    :param candidate_data_dict:
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
        # print(pos1)
        id1 = pos1["id"]
        sum_depth_up += find_position_depth(id1, candidate_data_dict)
    for pos2 in downstream_position:
        id2 = pos2["id"]
        sum_depth_down += find_position_depth(id2, candidate_data_dict)
    ave_depth_up = sum_depth_up / up_num
    ave_depth_down = sum_depth_down / down_num
    if ave_depth_up >= cutoff or ave_depth_down >= cutoff:
        return True
    else:
        return False


def find_candidate_align(up_down_locations, bam_path, assembly_path):
    """

    :param up_down_locations: a dictionary with the location of input mRNA, its up and down-stream mRNA positions.
    :return: up_down_alignment:
    position_info_align, alignment and un-alignment lists at the position.
    up_gene_align, alignment at upstream positions.
    down_gene_align, alignment at downstream positions
    """
    up_down_alignment = {}
    position_info = up_down_locations["position"]  # one dictionary, analysis both align and unaligned genomes
    upstream_position = up_down_locations["upstream_position"]  # list including dictionaries
    downstream_position = up_down_locations["downstream_position"]  # list including dictionaries

    # at the candidate position
    position_align_info = extract_align_seq_from_dict(bam_path, position_info)
    position_involved_assembly = find_aligned_assembly(position_align_info)
    position_uninvolved_assembly = find_not_aligned_assembly(position_involved_assembly, position_align_info,
                                                             assembly_path)

    upstream_position_involved_assembly = []
    for pos_up in upstream_position:
        position_align_info_up = extract_align_seq_from_dict(bam_path, pos_up)
        upstream_position_involved_assembly.append(position_align_info_up)

    downstream_position_involved_assembly = []
    for pos_down in downstream_position:
        position_align_info_down = extract_align_seq_from_dict(bam_path, pos_down)
        downstream_position_involved_assembly.append(position_align_info_down)

    up_down_alignment = {
        "position_assembly":  # at the candidate position
            {"involved_assembly": position_involved_assembly,
             "uninvolved_assembly": position_uninvolved_assembly
             },
        # aligned assemblies at the up and down positions
        "upstream_position_assembly": upstream_position_involved_assembly,
        "downstream_position_assembly": downstream_position_involved_assembly
    }
    return up_down_alignment


# Function: randomly sample n unaligned sequences from the up_down_alignment of a given gene,
# and then analyze the state at each position for each of them
def random_select_assembly(info_list, assembly_num, genome_exclusion = []):
    """
    Randomly sample n unaligned sequences at the interested mRNA position.
    This function is currently not used!
    :param up_down_alignment:
    :param n:
    :return: selected_assemblies, list of selected assemblies
    """
    assemblies = list(info_list.keys())

    assemblies_filtered = [x for x in assemblies if x not in genome_exclusion]

    number_assemblies = len(assemblies_filtered)
    if number_assemblies == 0:
        info_selected = {}

    elif 0 < number_assemblies <= assembly_num:
        info_selected = info_list

    elif number_assemblies > assembly_num:
        assemblies_selected = random.sample(assemblies_filtered, assembly_num)

        info_selected = {}

        for key, value in info_list.items():
            if key in assemblies_selected:
                info_selected[key] = value

    return info_selected


def find_candidate_status(assembly, reads):
    status = []
    for positions in reads:  # analyze each position
        seq_info_list = positions["seq_info"]
        genome_accessions = [d["Genome_accession"] for d in seq_info_list]
        if assembly in genome_accessions:
            status.append(True)
        else:
            status.append(False)
    return status


def find_assembly_status(selected_assemblies, upstream_reads, downstream_reads):
    assemblies_status = {}
    for assembly in selected_assemblies:  # test each selected assembly in up and downstream positions
        assembly_status = {
            "upstream_status": find_candidate_status(assembly, upstream_reads),
            "downstream_status": find_candidate_status(assembly, downstream_reads)
        }
        assemblies_status[assembly] = assembly_status

    return assemblies_status


# Finding whether the selected assemblies are involved in the positions.
def find_candidate_involvement(up_down_alignment):
    """

    :param up_down_alignment:
    :param
    :return:
    """
    # get the random assembly list to test status
    assemblies_ref_allele = up_down_alignment["position_assembly"]["involved_assembly"]
    assemblies_diver_allele = up_down_alignment["position_assembly"]["uninvolved_assembly"]

    # get the position information
    upstream_reads = up_down_alignment["upstream_position_assembly"]
    downstream_reads = up_down_alignment["downstream_position_assembly"]

    all_status = {
        "ref_allele": find_assembly_status(assemblies_ref_allele, upstream_reads, downstream_reads),
        "diver_allele": find_assembly_status(assemblies_diver_allele, upstream_reads, downstream_reads)
    }
    return all_status


def find_up_down_loci_one_status(one_status, up_down_locations, up_down_alignment):  #
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
    # print(upstream_align)
    downstream_align = up_down_alignment["downstream_position_assembly"]

    for genome_id, status in one_status.items():

        # print(genome_id,status)
        loci_selected = {}
        upstream_status = status["upstream_status"]
        downstream_status = status["downstream_status"]
        up_number = len(upstream_status)
        down_number = len(downstream_status)

        # incase situations without any findings
        all_keys = [
            'upstream_rank',
            'upstream_ref_position',
            'upstream_read_sequence_ID',
            'upstream_read_position',
            'downstream_rank',
            'downstream_ref_position',
            'downstream_read_sequence_ID',
            'downstream_read_position',
            'interval difference(%)',
            'seq_chromosome',
            'read_interval',
            'ref_interval']

        for key in all_keys:
            loci_selected[key] = False

        for j in range(1, up_number + 1):  # the i th nearest loci in the upstream positions
            # find upstream true locus
            status_loci_up = upstream_status[j - 1]
            if status_loci_up:  # if status_loci == True
                loci_selected["upstream_rank"] = j
                loci_selected["upstream_ref_position"] = upstream_positions[j - 1]

                # extract the information of the position in the selected assembly
                up_align_list = upstream_align[j - 1]["seq_info"]
                aligned_info_seq = []
                aligned_info_pos = {}
                for genome in up_align_list:
                    if genome["Genome_accession"] == genome_id:
                        # print(genome)
                        sequence_ID = genome["seq_ID"]
                        aligned_info_seq.append(sequence_ID)
                        # Assume that each sequence appears only once
                        aligned_info_pos[sequence_ID] = genome

                # print(aligned_info)
                loci_selected["upstream_read_sequence_ID"] = aligned_info_seq
                # print("aligned_info_seq",aligned_info_seq)
                loci_selected["upstream_read_position"] = aligned_info_pos
                break  # exit the loop once condition is met

        for i in range(1, down_number + 1):  # the i th nearest loci in the downstream positions
            # find downstream true locus
            status_loci_down = downstream_status[-i]
            if status_loci_down:  # if status_loci == True
                loci_selected["downstream_rank"] = i
                # The position information of the loci on the reference genome
                loci_selected["downstream_ref_position"] = downstream_positions[-i]
                ##
                # extract the information of the position in the selected assembly
                down_align_list = downstream_align[-i]["seq_info"]
                aligned_info_seq = []
                aligned_info_pos = {}
                for genome in down_align_list:
                    if genome["Genome_accession"] == genome_id:
                        # print(genome)
                        sequence_ID = genome["seq_ID"]
                        aligned_info_seq.append(sequence_ID)
                        # Assume that each sequence appears only once
                        aligned_info_pos[sequence_ID] = genome

                # print(aligned_info)
                loci_selected["downstream_read_sequence_ID"] = aligned_info_seq
                loci_selected["downstream_read_position"] = aligned_info_pos

                break  # exit the loop once condition is met

        up_seq = loci_selected["upstream_read_sequence_ID"]
        down_seq = loci_selected["downstream_read_sequence_ID"]

        if not up_seq or not down_seq:
            continue

        common_seq = set(up_seq).intersection(down_seq)

        if len(common_seq) == 0:
            continue

        loci_selected["seq_chromosome"] = common_seq

        # Assume there is only one shared seq in common_seq
        seq = list(common_seq)[0]

        # calculate interval between two reference genes in read genome
        read_up_position = loci_selected["upstream_read_position"][seq]
        read_down_position = loci_selected["downstream_read_position"][seq]
        orientation = read_up_position["orientation"]

        read_interval = abs(read_down_position["end"] - read_up_position["start"])
        loci_selected["read_interval"] = read_interval
        # print(read_up_position,'\t',read_down_position, interval)

        # calculate interval between two reference genes in reference genome
        ref_up_position = loci_selected["upstream_ref_position"]
        ref_down_position = loci_selected["downstream_ref_position"]
        ref_interval = ref_down_position["end"] - ref_up_position["start"]
        loci_selected["ref_interval"] = ref_interval

        interval_diff_pct = abs(ref_interval - read_interval) * 100 / ref_interval
        loci_selected["interval difference(%)"] = interval_diff_pct

        up_down_loci[genome_id] = loci_selected

        # print(loci_selected)

    return up_down_loci


def filter_up_down_loci(up_down_loci, interval_difference=250):
    filtered_up_down_loci = {}
    for genome_id, loci_selected in up_down_loci.items():
        if loci_selected["interval difference(%)"] is False:
            continue

        if loci_selected["interval difference(%)"] >= interval_difference:
            continue

        filtered_up_down_loci[genome_id] = loci_selected

    return filtered_up_down_loci


def find_up_down_loci(all_status, up_down_locations, up_down_alignment):
    ref_allele_status = all_status["ref_allele"]
    diver_allele_status = all_status["diver_allele"]

    ref_up_down_loci = find_up_down_loci_one_status(ref_allele_status, up_down_locations, up_down_alignment)
    ref_up_down_loci_filtered = filter_up_down_loci(ref_up_down_loci, 250)

    diver_up_down_loci = find_up_down_loci_one_status(diver_allele_status, up_down_locations, up_down_alignment)
    diver_up_down_loci_filtered = filter_up_down_loci(diver_up_down_loci, 250)
    all_up_down_loci = {
        "ref_up_down_loci": ref_up_down_loci_filtered,
        "diver_up_down_loci": diver_up_down_loci_filtered
    }

    return all_up_down_loci


def count_aligned_reads(all_up_down_loci):
    # count the number of reads completely aligned to the interval.
    ref_up_down_loci = all_up_down_loci["ref_up_down_loci"]
    diver_up_down_loci = all_up_down_loci["diver_up_down_loci"]
    ref_assembly_number = len(ref_up_down_loci.keys())
    diver_assembly_number = len(diver_up_down_loci.keys())
    aligned_reads_number = {
        "ref_assembly_number": ref_assembly_number,
        "diver_assembly_number": diver_assembly_number,
        "all_assembly_number": ref_assembly_number + diver_assembly_number
    }

    return aligned_reads_number


def analyze_all_candidate_position(selected_data, annotation_sorted, gene_data, bam_path, assembly_path, up_num, down_num,
                                   lower_limit,minimal_alignment,type_annotation):
    """
    Analyze the position data for each of the candidate positions.
    :param selected_data: the merged candidate data.
    :param annotation_sorted: genome annotation.
    :return: candidate_data_summary: a list, including analysis results of the position data
    [{
        "region_name": XM_750984.1,
        "candidate_data": candidate region, # information of this candidate, {seq: [merged_info1, ...], ...}
        "position_info": up_down_locations, # find upstream and downstream positions, dict
        "align_info": up_down_alignment, # find reads aligned to these positions,dict
        "status_info": all_status, # select assemblies and find whether they are present at these positions
        "up_down_loci": all_up_down_loci, # the selected locations
        "aligned_reads_number": aligned_reads_number
    }, ......]
    """
    candidate_data_summary = []
    for seq_id, candidates in selected_data.items():

        for candidate in candidates:
            sequence_id = seq_id
            start = candidate["start"]
            end = candidate["end"]
            # finding the up and downstream positions
            up_down_locations = find_position_seq(sequence_id, start, end, annotation_sorted, up_num, down_num, type_annotation)

            if up_down_locations == None:
                continue

            # check the mean depth of the positions
            gene_data_dict = dict_gene_data(gene_data)
            depth_status = filter_up_down_depth(up_down_locations, gene_data_dict, up_num, down_num, lower_limit)
            if depth_status == False:
                continue

            # finding the genes aligned to the positions
            up_down_alignment = find_candidate_align(up_down_locations, bam_path, assembly_path)

            # identify status of random genome assemblies at the positions
            all_status = find_candidate_involvement(up_down_alignment)

            # identify the most distant upstream and downstream loci
            all_up_down_loci = find_up_down_loci(all_status, up_down_locations, up_down_alignment)

            # Check whether both alleles exist
            ref_allele_info = all_up_down_loci["ref_up_down_loci"]
            diver_allele_info = all_up_down_loci["diver_up_down_loci"]

            if len(ref_allele_info) == 0 or len(diver_allele_info) == 0:
                continue

            aligned_reads_number = count_aligned_reads(all_up_down_loci)

            # Filter the complicated genomic region where not many reads completely aligned to
            all_aligned_reads_number = aligned_reads_number["all_assembly_number"]
            if all_aligned_reads_number < minimal_alignment:
                continue

            ref_aligned_number = aligned_reads_number["ref_assembly_number"]
            diver_aligned_number = aligned_reads_number["diver_assembly_number"]

            if ref_aligned_number < minimal_alignment/3 or diver_aligned_number < minimal_alignment/3:
                continue

            summary = {
                "region_name": candidate["region_name"],
                "candidate_data": candidate,
                "position_info": up_down_locations,
                "align_info": up_down_alignment,
                "status_info": all_status,
                "up_down_loci": all_up_down_loci,
                "aligned_reads_number": aligned_reads_number
            }

            candidate_data_summary.append(summary)

    return candidate_data_summary