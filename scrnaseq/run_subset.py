#!/usr/bin/env python3
"""
Wrapper to run subset reanalysis with environment variables.

Usage:
    python run_subset.py --data-dir /path/to/data --clusters 12,11,1,3 --seu-file QC_clustering/20250624/merged_clusters.rds
"""

import argparse
import os
import subprocess
import sys
from datetime import datetime


def main():
    parser = argparse.ArgumentParser(description="Run subset reanalysis")
    parser.add_argument("--data-dir", required=True, help="Main data directory")
    parser.add_argument("--seu-file", required=True, help="Seurat object path (relative to data-dir)")
    parser.add_argument("--clusters", required=True, help="Comma-separated cluster IDs")
    parser.add_argument("--subset-name", default="subset_analysis", help="Name for this subset")
    parser.add_argument("--out-dir", default="", help="Output directory (default: data-dir/subset_reanalysis/TIMESTAMP)")
    parser.add_argument("--conda-env", default="", help="Conda environment to use")
    parser.add_argument("--r-script", default="", help="R script to run (default: subset_reanalysis_vars.R)")
    args = parser.parse_args()

    # Set environment variables
    env = os.environ.copy()
    env["DATA_DIR"] = args.data_dir
    env["SEU_FILE"] = args.seu_file
    env["CLUSTERS"] = args.clusters
    env["SUBSET_NAME"] = args.subset_name
    
    if args.out_dir:
        env["OUT_DIR"] = args.out_dir
    
    # Find R script
    if args.r_script:
        r_script = args.r_script
    else:
        # Look in same directory as this script
        script_dir = os.path.dirname(os.path.abspath(__file__))
        r_script = os.path.join(script_dir, "subset_reanalysis_vars.R")
    
    if not os.path.exists(r_script):
        print(f"Error: R script not found: {r_script}")
        sys.exit(1)
    
    # Build command
    if args.conda_env:
        cmd = f"source ~/.bashrc && conda activate {args.conda_env} && Rscript {r_script}"
    else:
        cmd = f"Rscript {r_script}"
    
    print(f"Running subset reanalysis:")
    print(f"  DATA_DIR: {args.data_dir}")
    print(f"  SEU_FILE: {args.seu_file}")
    print(f"  CLUSTERS: {args.clusters}")
    print(f"  SUBSET_NAME: {args.subset_name}")
    print(f"  Script: {r_script}")
    print()
    
    # Run
    result = subprocess.run(cmd, shell=True, env=env, executable="/bin/bash")
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
