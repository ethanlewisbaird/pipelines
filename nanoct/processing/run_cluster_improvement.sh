#!/bin/bash
# Wrapper to run nanoCT cluster improvement with conda
# Finds conda, activates nanoCT env, runs the script

# Try common conda locations
CONDA_BASE=""
for p in /opt/conda /data/ebaird/miniconda3 ~/miniconda3 /home/ebaird/miniconda3; do
    if [ -f "$p/etc/profile.d/conda.sh" ]; then
        CONDA_BASE="$p"
        break
    fi
done

if [ -n "$CONDA_BASE" ]; then
    echo "Found conda at: $CONDA_BASE"
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate nanoCT
    echo "Using: $(which python)"
    python --version
else
    echo "No conda found, using system python3"
    which python3
fi

cd /data/ebaird/scentinel/nanoCT/20260522.nanoCT || exit 1

SCRIPT="/data/ebaird/pipelines/nanoct/processing/nanoct_cluster_improvement.py"
echo "Running: $SCRIPT"
python "$SCRIPT" 2>&1
EXIT_CODE=$?
echo "Exit code: $EXIT_CODE"
exit $EXIT_CODE
