from visualization_clinker import run_clinker_batch
from doublecheck_alignment import annotate_file_path

results_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae/results"
main_path = "/lustre/BIF/nobackup/leng010/test/aspergillus_oryzae"
sequence_path = f"{main_path}/extract_sequences"

clinker_output_path = f"{results_path}/clinker_results"
#run_clinker_batch(sequence_path, clinker_output_path)

# align sequence onto reference sequence to doublecheck and debug
output_main_path = f"{results_path}/sequence_alignments"
annotate_file_path(sequence_path, output_main_path)