#!/bin/sh
#SBATCH --job-name=PAE
#SBATCH --time=1:00:00
#SBATCH --error=%A.err
#SBATCH --output=%A.out
#SBATCH --mem=4G
#SBATCH --partition=lareauc_cpu,cpu
#SBATCH --output=pae_output.out
#SBATCH --error=pae_error.err

## User-provided
INPUT=$1  # Input directory 

## Paths
## ---- Paths ----
## A path set below (uncomment and edit its line) or already exported in your environment is never
## overwritten. Anything still empty is filled in from config.sh, written by ./configure.sh.
## Previous hardcoded paths:
## Replace with the python used by the CAR-PNN
# CARPNN_PYTHON="/data1/lareauc/users/chuh/miniconda3/envs/carpnn/bin/python"
# CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
CARPNN_CONFIG="${CARPNN_CONFIG:-$HOME/.config/carpnn/config.sh}"
if [ -f "$CARPNN_CONFIG" ]; then . "$CARPNN_CONFIG"; fi
: "${CARPNN_DIR:?not set - hardcode it above or run ./configure.sh}"
: "${CARPNN_PYTHON:?not set - hardcode it above or run ./configure.sh}"
PAE_UTIL_SCRIPT="${CARPNN_DIR}/workflows/Boltz/calculate_pae.py"

# Loop through all subdirectories in the specified directory
for sub_dir in "$INPUT"/*/; do
    echo "doing ${sub_dir}"
    ${CARPNN_PYTHON} ${PAE_UTIL_SCRIPT} \
    -dir ${sub_dir} 
done



