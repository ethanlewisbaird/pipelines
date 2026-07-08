import os
import subprocess
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description="Run pySCENIC pipeline 50 times.")
parser.add_argument('--output_dir', required=True, help="Base output directory")
parser.add_argument('--loom_file', required=True, help="Input loom file with expression data")
parser.add_argument('--tf_list', required=True, help="List of transcription factors")
parser.add_argument('--db_glob', required=True, help="Path to cisTarget databases")
parser.add_argument('--n_tasks', type=int, default=4, help="Number of parallel tasks")
parser.add_argument('--annotations', required=True, help="Path to annotations file")
args = parser.parse_args()

# Define file paths
loom_file = args.loom_file
tf_list = args.tf_list
db_glob = args.db_glob
base_output_dir = args.output_dir
n_tasks = args.n_tasks
annotations = args.annotations

# Run the pipeline 50 times
for i in range(1, 51):
    print(f"Running pySCENIC iteration {i}...")

    # Create a unique subdirectory for each run
    output_dir = os.path.join(base_output_dir, f'run_{i}')
    os.makedirs(output_dir, exist_ok=True)

    # Define output file paths with iteration number
    grn_output = os.path.join(output_dir, f'adjacencies_{i}.tsv')
    ctx_output = os.path.join(output_dir, f'ctx_{i}.csv')
    auc_output = os.path.join(output_dir, f'auc_mtx_{i}.csv')
    loom_output = os.path.join(output_dir, f'output_{i}.loom')

    # Step 1: GRN Inference
    subprocess.run([
        'pyscenic', 'grn', '--method', 'grnboost2',
        '--num_workers', str(n_tasks),
        '--output', grn_output,
        loom_file,
        tf_list
    ], check=True)

    # Step 2: Regulon Prediction
    subprocess.run([
        'pyscenic', 'ctx',
        '--annotations_fname', annotations,
        '--expression_mtx_fname', loom_file,
        '--output', ctx_output,
        '--mask_dropouts',
        '--num_workers', str(n_tasks),
        '--all_modules',
        grn_output,
        db_glob
    ], check=True)

    # Step 3: AUCell Analysis
    subprocess.run([
        'pyscenic', 'aucell',
        '--output', auc_output,
        '--num_workers', str(n_tasks),
        loom_file,
        ctx_output
    ], check=True)