#!/bin/bash
PARENT_DIR=$1

## ---- Paths ----
## A path set below (uncomment and edit its line) or already exported in your environment is never
## overwritten. Anything still empty is filled in from config.sh, written by ./configure.sh.
## Previous hardcoded paths:
# CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
CARPNN_CONFIG="${CARPNN_CONFIG:-$HOME/.config/carpnn/config.sh}"
if [ -f "$CARPNN_CONFIG" ]; then . "$CARPNN_CONFIG"; fi
: "${CARPNN_DIR:?not set - hardcode it above or run ./configure.sh}"
IPSAE_SCRIPT="${CARPNN_DIR}/workflows/Boltz/06_run_ipsae.sh"

PAE_CUTOFF=${2:-10}       # Default to 10 if not provided
DIST_CUTOFF=${3:-10}      # Default to 10 if not provided
MODEL_TYPE=${4:-boltz1}    # Default to 'boltz1' if not provided
OUTPUT_TYPE=${5:-pdb}    # Default to 'pdb' if not provided

## Find all subdirectories named "predictions" within the parent directory
## Note that it will find all subdirectories named "predictions" regardless of their directory depth
find "$PARENT_DIR" -type d -name "predictions" | while IFS= read -r pred_dir; do
  echo "Launching $IPSAE_SCRIPT for: $pred_dir"
  sbatch "$IPSAE_SCRIPT" "$pred_dir" ${PAE_CUTOFF} ${DIST_CUTOFF} ${MODEL_TYPE} ${OUTPUT_TYPE}
done

