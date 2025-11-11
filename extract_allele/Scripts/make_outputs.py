import csv
import pandas as pd
import os
import glob
from Bio import SeqIO
from load_reference import read_gff
from analyze_position import random_select_assembly
from augustus_annotation import annotate_file_path


#############################################################
# extract the sequence of the allele genomic region for each genes.

def extract_allele_sequence(genome_assembly_path, candidate_gene, genome_accession, seq, start, end, orientation,
                            output_path):
    """
    Extract the sequence of specific position of a genome assembly
    :param assembly_dir:
    :param candidate_gene:
    :param genome_accession:
    :param seq:
    :param start:
    :param end:
    :param orientation:
    :param output_path:
    :return:
    """
    # extract sequence
    target_seq = None
    for record in SeqIO.parse(genome_assembly_path, "fasta"):
        # print(record.id)
        if seq in record.id:
            # print("find")
            target_seq = record.seq[start - 1:end]  # Biopython is 0-based
            break

    if target_seq is None:
        warning = f"Warning: Sequence {seq} not found in {genome_accession}"
        print(warning)
        return warning

    if orientation.lower() == "reverse":
        target_seq = target_seq.reverse_complement()

    # output path
    output_dir = os.path.join(output_path, candidate_gene)
    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, f"{candidate_gene}_{genome_accession}_{seq}-{start}-{end}.fa")

    # 4. save results
    with open(output_file, "w") as f:
        f.write(f">{genome_accession}_{seq}:{start}-{end}\n")
        f.write(str(target_seq) + "\n")

    print(f"Sequence saved to {output_file}")
    return output_file


def find_genome_assembly_path(assembly_dir, genome):
    assembly_pattern = os.path.join(assembly_dir, f"{genome}*.fna")
    matches = glob.glob(assembly_pattern)
    if not matches:
        warning = f"Warning: No genome assembly found for {genome}"
        print(warning)
        # return warning
    genome_assembly_path = matches[0]
    return genome_assembly_path


def extract_region_seq(allele_info, region_name, label, assembly_dir, output_path, extend, assembly_num):
    """
    Extract the sequence of specific position of a genome assembly
    :param allele_info:
    :param region_name:
    :param assembly_dir:
    :param output_path:
    :param extend:
    :param assembly_num:
    :return:
    """

    allele_info_selected = random_select_assembly(allele_info, assembly_num)

    for genome, info in allele_info_selected.items():
        seq = info["seq_chromosome"]
        seq_info = list(seq)[0]  # assume there is only one seq

        upstream_read_position = info["upstream_read_position"][seq_info]
        downstream_read_position = info["downstream_read_position"][seq_info]

        start_read = min(upstream_read_position["start"], downstream_read_position["start"])
        end_read = max(upstream_read_position["end"], downstream_read_position["end"])
        orientation = upstream_read_position["orientation"]

        # if end_read - start_read > 100000:
        # continue

        # extend the start-end interval
        start_read = max(1, start_read - extend)
        end_read = end_read + extend  # require modified

        # finding genome assembly path
        genome_assembly_path = find_genome_assembly_path(assembly_dir, genome)

        extract_allele_sequence(
            genome_assembly_path,
            region_name,
            f"{label}_{genome}",
            seq_info,
            start_read,
            end_read,
            orientation,
            output_path
        )

    return True


def find_allele_sequence_inbetween(assembly_dir, candidate_data_summary, output_path, extend, assembly_num):
    """

    :param reference_genome:
    :param assembly_dir:
    :param candidate_data_summary:
    summary = {
                "region_name": candidate["region_name"],
                "position_info": up_down_locations,
                "align_info": up_down_alignment,
                "status_info": all_status,
                "up_down_loci": all_up_down_loci,
                "aligned_reads_number": aligned_reads_number
            }
    :param output_path:
    :return:
    """
    for summary in candidate_data_summary:  # the summary information of each candidate gene
        region_name = summary["region_name"]

        ref_allele_info = summary["up_down_loci"]["ref_up_down_loci"]
        diff_allele_info = summary["up_down_loci"]["diff_up_down_loci"]

        ref_extract = extract_region_seq(ref_allele_info, region_name, "ref_allele", assembly_dir, output_path, extend,
                                         assembly_num)
        diff_extract = extract_region_seq(diff_allele_info, region_name, "diff_allele", assembly_dir, output_path,
                                          extend, assembly_num)

    return output_path

def extract_reference_allele_Augustus(candidate_data_summary, reference_genome, gff_path, output_path, extend,
                             ref_assembly):
    """

    :param candidate_data_summary:
    :param reference_genome:
    :param annotation_sorted: {seq_ID: [row_dict, ...]}
    :param output_path:
    :return:
    """
    # load annotation
    annotation_sorted = read_gff(gff_path)

    for summary in candidate_data_summary:
        region_name = summary["region_name"]
        up_down_locations = summary["position_info"]

        seq_info_ref = up_down_locations["position"]["seq_ID"]
        start_ref = up_down_locations["upstream_position"][0]["start"]
        end_ref = up_down_locations["downstream_position"][-1]["end"]

        # calculate the start and end position of the extraction region
        # start position
        start = max(1, start_ref - extend)

        # end position
        annotation_info = annotation_sorted[seq_info_ref]
        annotation_end = annotation_info[-1]["end"]
        end = min(end_ref + extend, annotation_end)

        # extract the sequence from the reference genome
        extract_allele_sequence(
            ref_assembly,
            region_name,
            f"reference_genome_Augustus_{reference_genome}",
            seq_info_ref,
            start,
            end,
            "forward",
            output_path
        )

def extract_reference_allele(candidate_data_summary, reference_genome, gff_path, output_path, extend,
                             ref_assembly):
    """

    :param candidate_data_summary:
    :param reference_genome:
    :param annotation_sorted: {seq_ID: [row_dict, ...]}
    :param output_path:
    :return:
    """
    # load annotation
    annotation_sorted = read_gff(gff_path)

    for summary in candidate_data_summary:
        region_name = summary["region_name"]
        up_down_locations = summary["position_info"]

        seq_info_ref = up_down_locations["position"]["seq_ID"]
        start_ref = up_down_locations["upstream_position"][0]["start"]
        end_ref = up_down_locations["downstream_position"][-1]["end"]

        # calculate the start and end position of the extraction region
        # start position
        start = max(1, start_ref - extend)

        # end position
        annotation_info = annotation_sorted[seq_info_ref]
        annotation_end = annotation_info[-1]["end"]
        end = min(end_ref + extend, annotation_end)
        genome_accession = f"reference_genome_{reference_genome}"

        # extract the sequence from the reference genome
        extract_allele_sequence(
            ref_assembly,
            region_name,
            genome_accession,
            seq_info_ref,
            start,
            end,
            "forward",
            output_path
        )

        seq_name = f"{genome_accession}_{seq_info_ref}:{start}-{end}"

        # --- Read and filter GFF annotations ---
        # annotation_sorted {seq_ID: [row_dict, ...]}
        extract_annotation = []
        for annotation in annotation_info:
            start_anno = annotation["start"]
            end_anno = annotation["end"]

            if start_anno >= start and end_anno <= end:
                annotation_new = annotation.copy()
                annotation_new["start"] = annotation["start"] - start + 1
                annotation_new["end"] = annotation["end"] - start + 1
                extract_annotation.append(annotation_new)

        # create output path for gff3 file
        #output_dir = os.path.join(output_path, region_name, "reference_annotation")
        output_dir = os.path.join(output_path, region_name)
        os.makedirs(output_dir, exist_ok=True)
        output_file = f"{output_dir}/{region_name}_{genome_accession}_{seq_info_ref}-{start}-{end}.gff3"

        # write in gff3 formate file
        with open(output_file, "w", encoding="utf-8", newline="") as out:
            out.write(f"##gff-version 3 {seq_name}\n")

            writer = csv.writer(out, delimiter="\t", lineterminator="\n")

            for ann in extract_annotation:
                type_ = ann.get("type", ".")

                if type_ == "mRNA":
                    type_ = "transcript"

                seq_ID = seq_name
                source = ann.get("source", ".")
                start = str(ann.get("start", "."))
                end = str(ann.get("end", "."))
                score = ann.get("score", ".")
                strand = ann.get("strand", ".")
                phase = ann.get("phase", ".")
                attributes = ann.get("attributes", ".")

                # write in file
                writer.writerow([seq_ID, source, type_, start, end, score, strand, phase, attributes])

        print(f"GFF3 file successfully saved: {output_file}")


def find_final_candidates(CDS_dict, candidate_data_summary, candidate_data, genome_num):
    """

    :param candidate_data_summary:
    :param candidate_data:
    :return:
    """
    final_candidates = []

    for summary in candidate_data_summary:
        region_info = summary["position_info"]["position"]
        region_name = summary["region_name"]
        # print("region_info", region_info)
        region_seq = region_info["seq_ID"]
        region_start = region_info["start"]
        region_end = region_info["end"]

        for row in candidate_data:
            candidate_seq = row["seq_ID"]
            candidate_start = row["start"]
            candidate_end = row["end"]

            if candidate_seq == region_seq:
                if region_start <= candidate_start <= region_end and region_start <= candidate_end <= region_end:
                    candidate_info = {
                        "genomic_region": region_name,
                        "id": row["id"],
                        "protein_id": CDS_dict[row["id"]],
                        "locus_tag": row["locus_tag"],
                        "type": row["type"],
                        "seq_ID": row["seq_ID"],
                        "start": row["start"],
                        "end": row["end"],
                        "gene_length": row["end"]-row["start"],
                        "mean_depth": row["depth"],
                        "depth_ratio (%)": row["depth"] * 100 / genome_num,
                        "gene_info": row
                    }
                    final_candidates.append(candidate_info)

    return final_candidates


def save_final_candidates(final_candidates, output_path):
    """
    Save a list of candidate dictionaries to an Excel file.
    Each dictionary becomes one row, and 'gene_info' (a nested dict)
    is converted into a "key:value" comma-separated string.
    """

    # Ensure the output directory exists
    os.makedirs(output_path, exist_ok=True)

    processed_data = []

    for item in final_candidates:
        # Copy the item to avoid modifying the original list
        row = item.copy()

        # Convert 'gene_info' dict to a readable "key:value" string
        if isinstance(row.get("gene_info"), dict):
            row["gene_info"] = ", ".join([f"{k}: {v}" for k, v in row["gene_info"].items()])

        processed_data.append(row)

    # Create a DataFrame (keys automatically become column headers)
    df = pd.DataFrame(processed_data)

    # Build the output file path
    output_file = os.path.join(output_path, "final_candidates.xlsx")

    # Save to Excel (requires 'openpyxl' installed)
    df.to_excel(output_file, index=False)

    print(f"Final candidate data (genes) saved to: {output_file}")
    return True

def extract_candidates(CDS_dict, candidate_data_summary, main_path, candidate_data, genome_num):

    results_path = f"{main_path}/results"

    # find and save final candidate genes and related information
    final_candidates = find_final_candidates(CDS_dict,candidate_data_summary, candidate_data, genome_num)
    save_final_candidates(final_candidates, results_path)

    return results_path

def extract_sequences(candidate_data_summary, reference_genome, gff_path, main_path, extend, ref_assembly,assembly_dir,assembly_num, augustus_species):

    sequence_path = f"{main_path}/extract_sequences"
    # extract sequence and annotation from other genomes
    sequence_path = find_allele_sequence_inbetween(assembly_dir, candidate_data_summary, sequence_path, extend, assembly_num)

    #extract sequence from reference genome for Augustus annotation
    extract_reference_allele_Augustus(candidate_data_summary, reference_genome, gff_path, sequence_path, extend, ref_assembly)

    # make AUGUSTUS annotation
    annotate_file_path(sequence_path, augustus_species)

    # extract both sequence and annotation from the reference genome
    extract_reference_allele(candidate_data_summary, reference_genome, gff_path, sequence_path, extend, ref_assembly)

    return sequence_path
