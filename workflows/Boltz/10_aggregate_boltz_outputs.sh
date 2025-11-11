#!/bin/bash

# --- Configuration ---
# Set the top-level directory containing all the prediction subdirectories.
# This script will iterate through each subdirectory inside this path.
# e.g. 
PRED_DIR=$1

# Set the directory where the output CSV files will be saved.
# The script will create this directory if it doesn't exist.
# e.g.
OUTPUT_DIR=$2

## Paths
## Use the CARPNN python
CARPNN_PYTHON="/data1/lareauc/users/chuh/miniconda3/envs/carpnn/bin/python"
CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
PYTHON_SCRIPT_PATH="${CARPNN_DIR}/workflows/Boltz/aggregate_boltz_outputs.py"

# --- Script Logic ---
# Check if the main prediction directory exists
if [ ! -d "$PRED_DIR" ]; then
    echo "Error: Prediction directory '$PRED_DIR' not found."
    exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Find all subdirectories within the main predictions directory
# The -maxdepth 1 ensures we only look one level deep
# The -type d finds only directories
for subdir in $(find "$PRED_DIR" -maxdepth 1 -type d -printf '%P\n' | grep -v '^\.$'); do
    # Get the full path of the subdirectory
    FULL_SUBDIR_PATH="$PRED_DIR/$subdir"
    
    # Define the output file name for the CSV based on the subdirectory name
    OUTPUT_CSV="$OUTPUT_DIR/${subdir}.csv"
    
    # Submit a Slurm job for each subdirectory
    sbatch --job-name="load_${subdir}" \
           --output="${OUTPUT_DIR}/slurm_output_${subdir}.log" \
           --error="${OUTPUT_DIR}/slurm_error_${subdir}.log" \
           --partition=lareauc_cpu,cpu \
           --time=1:00:00 \
           --mem=16G \
           --wrap="python ${PYTHON_SCRIPT_PATH} ${FULL_SUBDIR_PATH} --output_csv ${OUTPUT_CSV}"
    
    echo "Submitted Slurm job for directory: $FULL_SUBDIR_PATH"
done
