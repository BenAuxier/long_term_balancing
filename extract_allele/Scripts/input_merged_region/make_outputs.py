import csv
import pandas as pd
import os
import glob
import subprocess
from Bio import SeqIO
from load_reference import read_gff
from load_reference import read_gff_augustus
from analyze_position import random_select_assembly

#################################################################

def find_region_gene(candidate_data_summary, annotation_refseq, annotation_augustus, region_output_file, refseq_candidate_file):
    """

    :param candidate_data_summary:

    [{
        "region_name": name,
        "candidate_data": candidate region, # information of this candidate, {seq: [merged_info1, ...], ...}
        "position_info": up_down_locations, # find upstream and downstream positions, dict
        "up_down_loci": all_up_down_loci, # the selected locations
    }, ......]

    [{"region_name": candidate["region_name"],
     "position_info": up_down_locations,
     "up_down_loci": all_up_down_loci,
    },...]
    :param annotation_refseq: {seq_ID: [row_dict, ...]}
    :param annotation_augustus: {seq_ID: [row_dict, ...]}
    :param output_path:
    :return:
    """
    # check output path
    os.makedirs(os.path.dirname(region_output_file), exist_ok=True)
    os.makedirs(os.path.dirname(refseq_candidate_file), exist_ok=True)

    summary_genes = []
    all_refseq_genes = []
    for summary in candidate_data_summary:
        region_name = summary["region_name"]
        region_info = summary["position_info"]["position"]
        region_seq = region_info["seq_ID"]
        region_start = region_info["start"]
        region_end = region_info["end"]

        annotation_refseq_seq = annotation_refseq[region_seq]
        annotation_augustus_seq = annotation_augustus[region_seq]

        refseq_gene = []
        augustus_gene = []

        for row in annotation_refseq_seq:
            candidate_start = row["start"]
            candidate_end = row["end"]
            candidate_id = row["id"]

            if (region_start <= candidate_start <= region_end) or (region_start <= candidate_end <= region_end):
                refseq_gene.append(candidate_id)
                all_refseq_genes.extend(refseq_gene)

        for row in annotation_augustus_seq:
            candidate_start = row["start"]
            candidate_end = row["end"]
            candidate_id = row["id"]

            if (region_start <= candidate_start <= region_end) or (region_start <= candidate_end <= region_end):
                augustus_gene.append(candidate_id)

        summary_genes.append([region_name, refseq_gene, augustus_gene])

    with open(region_output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["genomic_region", "RefSeq_genes", "Augustus_genes"])
        writer.writerows(summary_genes)
    print(f"Genomic region file written to {region_output_file}")

    all_refseq_genes = list(set(all_refseq_genes))
    with open(refseq_candidate_file, "w", newline="", encoding="utf-8") as f:
        for item in all_refseq_genes:
            f.write(str(item) + "\n")

    print(f"Refseq candidate genes written to {refseq_candidate_file}")

    return summary_genes



###############################################################
# augustus_annotation
def run_augustus_on_fasta(fa_path, augustus_species, gff3_status = "off", suffix = ""):
    """
    Run AUGUSTUS on a single FASTA file and save the output in GFF3 format.
    """
    base = os.path.splitext(os.path.basename(fa_path))[0]
    dir_path = os.path.dirname(fa_path)

    if gff3_status == "on":
        output_path = os.path.join(dir_path, f"{base}{suffix}.gff3")
    else:
        output_path = os.path.join(dir_path, f"{base}{suffix}.gff")

    # check whether the annotation already exist
    if os.path.exists(output_path):
        print(f"[Skip] {output_path} already exists. Skipping AUGUSTUS run.")
        return output_path

    print(f"Running AUGUSTUS on {fa_path} ...")

    # Build command
    cmd = [
        "augustus",
        f"--species={augustus_species}",
        f"--gff3={gff3_status}",
        fa_path
    ]

    # Run command and write output
    with open(output_path, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, check=True)

    return output_path

def annotate_file_path(input_dir,augustus_species, gff3_status = "off"):
    """
    Recursively find all .fa files under INPUT_DIR and run AUGUSTUS for each.
    """

    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".fa"):
                fa_path = os.path.join(root, file)
                try:
                    run_augustus_on_fasta(fa_path,augustus_species, gff3_status)
                except subprocess.CalledProcessError as e:
                    print(f"Error running AUGUSTUS on {fa_path}: {e}")
    return True


#############################################################
# extract the sequence of the allele genomic region for each genes.

def extract_allele_sequence(genome_assembly_path, genome_accession, region_name,
                            seq, start, end, orientation, output_file):
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

    sequence_name = f"{region_name}__{genome_accession}__{seq}__{start}-{end}"

    # save results
    with open(output_file, "w") as f:
        f.write(f">{sequence_name}\n")
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


def extract_region_seq(allele_info, region_name, label, assembly_dir, output_path, extend, assembly_num, genome_exclusion = []):
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

    allele_info_selected = random_select_assembly(allele_info, assembly_num, genome_exclusion)

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

        output_dir = os.path.join(output_path, region_name)
        os.makedirs(output_dir, exist_ok=True)

        file_name = f"{region_name}__{genome}__{seq_info}__{start_read}-{end_read}__{label}"
        output_file = os.path.join(output_dir, f"{file_name}.fa")

        extract_allele_sequence(
            genome_assembly_path,
            genome,
            region_name,
            seq_info,
            start_read,
            end_read,
            orientation,
            output_file
        )

    return True


def find_allele_sequence_inbetween(reference_genome, assembly_dir, candidate_data_summary, output_path, extend, assembly_num):
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
        diver_allele_info = summary["up_down_loci"]["diver_up_down_loci"]

        # not extract the reference genome
        genome_exclusion = [reference_genome]

        ref_extract = extract_region_seq(ref_allele_info, region_name, "ref_allele", assembly_dir, output_path, extend,
                                         assembly_num, genome_exclusion)
        diver_extract = extract_region_seq(diver_allele_info, region_name, "diver_allele", assembly_dir, output_path,
                                          extend, assembly_num, genome_exclusion)

    return output_path

def extract_sequence_interpro(extracted_sequences, assembly_dir, candidate_data_summary, output_path, extend, assembly_num):
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

        region_extracted = extracted_sequences[region_name]

        ref_extracted = region_extracted["ref_allele"]
        diver_extracted = region_extracted["diver_allele"]

        ref_allele_info = summary["up_down_loci"]["ref_up_down_loci"]
        diver_allele_info = summary["up_down_loci"]["diver_up_down_loci"]

        ref_extract = extract_region_seq(ref_allele_info, region_name, "ref_allele", assembly_dir, output_path, extend,
                                         assembly_num, ref_extracted)
        diver_extract = extract_region_seq(diver_allele_info, region_name, "diver_allele", assembly_dir, output_path,
                                          extend, assembly_num,diver_extracted)

    return output_path


def extract_reference_seq_augustus(candidate_data_summary, reference_genome, gff_path, output_path, extend,
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

        output_dir = os.path.join(output_path, region_name)
        os.makedirs(output_dir, exist_ok=True)
        label = "re-annotate_AUGUSTUS"
        file_name = f"{region_name}__{reference_genome}__{seq_info_ref}__{start}-{end}__{label}"

        output_file = os.path.join(output_dir, f"{file_name}.fa")
        extract_allele_sequence(
            ref_assembly,
            reference_genome,
            region_name,
            seq_info_ref,
            start,
            end,
            "forward",
            output_file
        )

def extract_reference_allele(candidate_data_summary, reference_genome, annotation_sorted, output_path, extend,
                             ref_assembly, label):
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

        # rename
        output_dir = os.path.join(output_path, region_name)
        os.makedirs(output_dir, exist_ok=True)

        file_name = f"{region_name}__{reference_genome}__{seq_info_ref}__{start}-{end}__{label}"

        output_file_fa = os.path.join(output_dir,f"{file_name}.fa")
        extract_allele_sequence(
            ref_assembly,
            reference_genome,
            region_name,
            seq_info_ref,
            start,
            end,
            "forward",
            output_file_fa
        )

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
        output_file_gff3 = f"{output_dir}/{file_name}.gff3"

        # write in gff3 formate file
        with open(output_file_gff3, "w", encoding="utf-8", newline="") as out:
            out.write(f"##gff-version 3 {file_name}\n")

            writer = csv.writer(out, delimiter="\t", lineterminator="\n")

            for ann in extract_annotation:
                type_ann = ann.get("type", ".")

                if type_ann == "mRNA":
                    type_ann = "transcript"

                # seq name should be same as the name in fasta file title
                seq_ID = f"{region_name}__{reference_genome}__{seq_info_ref}__{start}-{end}"
                source = ann.get("source", ".")
                start_anno = str(ann.get("start", "."))
                end_anno = str(ann.get("end", "."))
                score = ann.get("score", ".")
                strand = ann.get("strand", ".")
                phase = ann.get("phase", ".")
                attributes = ann.get("attributes", ".")

                # write in file
                writer.writerow([seq_ID, source, type_ann, start_anno, end_anno, score, strand, phase, attributes])

        print(f"GFF3 file successfully saved: {output_file_gff3}")



def extract_reference_allele_ref(candidate_data_summary, reference_genome, gff_path, output_path, extend,
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

    label = "ref_RefSeq"

    extract_reference_allele(candidate_data_summary, reference_genome,
                             annotation_sorted, output_path, extend,
                             ref_assembly, label)


def extract_reference_allele_augustus(candidate_data_summary, reference_genome, gff_path_augustus, output_path, extend,
                             ref_assembly):
    """

    :param candidate_data_summary:
    :param reference_genome:
    :param annotation_sorted: {seq_ID: [row_dict, ...]}
    :param output_path:
    :return:
    """
    # load annotation
    annotation_sorted_augustus = read_gff_augustus(gff_path_augustus)

    label = "ref_AUGUSTUS"

    extract_reference_allele(candidate_data_summary, reference_genome,
                             annotation_sorted_augustus, output_path, extend,
                             ref_assembly, label)


def find_reference_gene(annotation_sorted, seq, start, end, ID_dict):
    """

    :param annotation_sorted:
    :param seq:
    :param start:
    :param end:
    :param ID_dict:
    :return:
    """
    #find the gene id included in a selected genomic region
    seq_annotation = annotation_sorted[seq]

    genes_included = []
    CDS_included = []
    for annotation in seq_annotation:
        annotation_id = annotation["id"]
        annotation_start = annotation["start"]
        annotation_end = annotation["end"]
        if (end < annotation_start) or (annotation_end < start):
            continue

        genes_included.append(annotation_id)

        if not ID_dict:
            continue

        CDS_id = ID_dict[annotation_id]
        CDS_included.append(CDS_id)

    return genes_included, CDS_included

def find_final_candidates(candidate_data_summary, candidate_data, genome_num, annotation_sorted, CDS_dict):
    """

    :param candidate_data_summary:
    :param candidate_data:
    :return:
    """
    final_candidates = []

    for summary in candidate_data_summary:
        region_info = summary["position_info"]["position"]
        region_name = summary["region_name"]
        region_seq = region_info["seq_ID"]
        region_start = region_info["start"]
        region_end = region_info["end"]

        for row in candidate_data:
            candidate_seq = row["seq_ID"]
            candidate_start = row["start"]
            candidate_end = row["end"]

            if candidate_seq == region_seq:
                if (region_start <= candidate_start <= region_end) or (region_start <= candidate_end <= region_end):

                    genes_included, CDS_included = find_reference_gene(annotation_sorted, candidate_seq,
                                                         candidate_start,candidate_end, CDS_dict)

                    candidate_info = {
                        "genomic_region": region_name,
                        "id": row.get("id", "."),
                        "RefSeq_genes": genes_included,
                        "RefSeq_CDS": CDS_included,
                        "RefSeq_genes_number": len(CDS_included),
                        "seq_ID": row.get("seq_ID", "."),
                        "start": row.get("start", "."),
                        "end": row.get("end", "."),
                        "gene_length": row["end"]-row["start"],
                        "mean_depth": row.get("depth", "."),
                        "depth_ratio (%)": row["depth"] * 100 / genome_num,
                        "gene_info": row
                    }

                    final_candidates.append(candidate_info)

    return final_candidates


def save_final_candidates(final_candidates, output_file):
    """
    Save a list of candidate dictionaries to an Excel file.
    Each dictionary becomes one row, and 'gene_info' (a nested dict)
    is converted into a "key:value" comma-separated string.
    """

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

    # Save to Excel (requires 'openpyxl' installed)
    df.to_excel(output_file, index=False)

    print(f"Final candidate data (genes) saved to: {output_file}")
    return output_file

def extract_candidates(candidate_data_summary, output_file, candidate_data, genome_num, annotation_sorted, CDS_dict):
    """

    :param candidate_data_summary:
    :param output_file:
    :param candidate_data:
    :param genome_num:
    :param annotation_sorted:
    :param CDS_dict:
    :return:
    """
    results_path = os.path.dirname(output_file)
    os.makedirs(results_path, exist_ok=True)

    # find and save final candidate genes and related information
    final_candidates = find_final_candidates(candidate_data_summary, candidate_data, genome_num,
                                             annotation_sorted, CDS_dict)

    results_path = save_final_candidates(final_candidates, output_file)

    return results_path

def extract_sequences(candidate_data_summary, reference_genome, ref_gff, ref_gff3_augustus, sequence_path, extend, ref_assembly,assembly_dir,assembly_num, augustus_species):
    """

    :param candidate_data_summary:
    :param reference_genome:
    :param ref_gff:
    :param ref_gff3_augustus:
    :param sequence_path:
    :param extend:
    :param ref_assembly:
    :param assembly_dir:
    :param assembly_num:
    :param augustus_species:
    :return:
    """

    os.makedirs(sequence_path, exist_ok=True)
    # extract sequence and annotation from other genomes
    sequence_path = find_allele_sequence_inbetween(reference_genome, assembly_dir, candidate_data_summary, sequence_path, extend, assembly_num)

    #extract sequence from reference genome for Augustus annotation
    # seems not necessary if AUGUSTUS annotated reference is used
    #extract_reference_seq_augustus(candidate_data_summary, reference_genome, ref_gff, sequence_path, extend, ref_assembly)

    # make AUGUSTUS annotation
    annotate_file_path(sequence_path, augustus_species,"on")

    # extract both sequence and annotation from the AUGUSTUS annotated reference
    extract_reference_allele_augustus(candidate_data_summary, reference_genome, ref_gff3_augustus, sequence_path, extend,ref_assembly)

    # extract both sequence and annotation from the (RefSeq) reference genome
    extract_reference_allele_ref(candidate_data_summary, reference_genome, ref_gff, sequence_path, extend, ref_assembly)

    return sequence_path

def check_sequence(sequence_path):
    """

    :param sequence_path:
    :return: extracted_sequences
    {"g3347.t1-g3348.t1":{'ref_allele': ['GCA_020502525.1', ...],
    'diver_allele': ['GCA_001643665.1', ...], ...}

    """

    extracted_sequences = {}
    sequence_files = os.listdir(sequence_path)
    for genomic_region in sequence_files:
        region_dict = {
            "ref_allele": [],
            "diver_allele": []
        }

        genomic_path = f"{sequence_path}/{genomic_region}"

        # search for the file including genome accession in its name
        assembly_pattern = os.path.join(genomic_path, "*.fa")
        matches = glob.glob(assembly_pattern)

        if not matches:
            warning = f"Warning: No genome assembly found for {genomic_region}"
            print(warning)
            # return warning
        matches = [os.path.basename(f) for f in matches]

        for sequence_file in matches:
            file_names = sequence_file.strip(".fa").split("__")

            assembly = file_names[1]
            status = file_names[-1]
            if status == "ref_allele":
                region_dict["ref_allele"].append(assembly)
                continue
            elif status == "diver_allele":
                region_dict["diver_allele"].append(assembly)
                continue

        extracted_sequences[genomic_region] = region_dict

    return extracted_sequences

def copy_sequences(sequence_path, sequence_interpro):
    sequence_files = os.listdir(sequence_path)
    for genomic_region in sequence_files:
        genomic_path_old = f"{sequence_path}/{genomic_region}"
        genomic_path_new = f"{sequence_interpro}/{genomic_region}"

        if not os.path.exists(genomic_path_old):
            continue

        # if the new path do not exist (when all possible assemblies already extracted)
        if not os.path.exists(genomic_path_new):
            os.makedirs(genomic_path_new, exist_ok=True)

        # search for the file including genome accession in its name
        assembly_pattern = os.path.join(genomic_path_old, "*")
        matches = glob.glob(assembly_pattern)

        for file_path in matches:
            cmd = ["cp", file_path, genomic_path_new]
            subprocess.run(cmd)

def extract_sequences_interpro(sequence_interpro, candidate_data_summary, sequence_path, extend, assembly_dir, assembly_num, augustus_species):
    """

    :param sequence_interpro:
    :param candidate_data_summary:
    :param sequence_path:
    :param extend:
    :param assembly_dir:
    :param assembly_num:
    :param augustus_species:
    :return:
    """
    os.makedirs(sequence_interpro, exist_ok=True)

    #check extracted sequences
    extracted_sequences = check_sequence(sequence_path)

    # extract sequence and annotation from other genomes
    extract_sequence_interpro(extracted_sequences, assembly_dir, candidate_data_summary, sequence_interpro, extend, assembly_num)

    # make AUGUSTUS annotation
    annotate_file_path(sequence_interpro, augustus_species, "on")

    #copy other previously extracted sequences and annotations to filepath
    copy_sequences(sequence_path, sequence_interpro)

    return True




if __name__ == "__main__":
    # information
    reference_genome = "GCF_000002655.1"  # genome annotation should be GCF version
    species = "aspergillus_fumigatus"
    augustus_species = "aspergillus_fumigatus"

    type_annotation_ref = "mRNA"  # type of annotation used in depth calculation, the third column
    type_annotation_augustus = "transcript"  # type of annotation used in depth calculation, the third column

    key_words = None  # the keywords that have to be included in the annotation
    ID_ref_label = "locus_tag"  #
    ID_augustus_label = "gene_id"  # this is the key that the gene/mRNA id follows in gff file

    # file paths, including all input files
    base_path = "/lustre/BIF/nobackup/leng010/test"
    # path to specific species
    main_path = f"{base_path}/{species}"
    assembly_dir = f"{main_path}/genome_assemblies"
    ref_path = f"{main_path}/reference_genome"
    ref_assembly = f"{ref_path}/{reference_genome}_genomic.fna"
    gff_refseq = f"{ref_path}/{reference_genome}_genomic.gff"
    gff_refseq_filtered = f"{gff_refseq[:-4]}_{type_annotation_ref}.gff"
    gff_augustus = f"{ref_path}/{reference_genome}_genomic_AUGUSTUS.gff"
    gff3_augustus = f"{ref_path}/{reference_genome}_genomic_AUGUSTUS.gff3"
    gff_augustus_filtered = f"{gff_augustus[:-4]}_{type_annotation_augustus}.gff"
    bam_path = f"{main_path}/alignment"
    bam_file = f"{main_path}/alignment/alignment_{species}.sorted.bam"
    filtered_region_nucleotides = f"{main_path}/depth_calculation/filtered_region_nucleotides.txt"

    sequence_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_fumigatus/extract_sequences"
    extracted_sequences = check_sequence(sequence_path)

    ###########################################################################

