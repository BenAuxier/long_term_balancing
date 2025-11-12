from visualization_clinker import run_clinker_batch
sequence_path = "/lustre/BIF/nobackup/leng010/test/magnaporthe_grisea/extract_sequences"
results_path = "/lustre/BIF/nobackup/leng010/test/magnaporthe_grisea/results"
clinker_output_dir = run_clinker_batch(sequence_path, results_path)

