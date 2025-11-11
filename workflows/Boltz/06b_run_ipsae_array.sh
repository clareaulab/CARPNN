#!/bin/bash
PARENT_DIR=$1

CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
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

