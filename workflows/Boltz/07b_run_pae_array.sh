#!/bin/sh
PARENT_DIR=$1

CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
PAE_SCRIPT="${CARPNN_DIR}/workflows/Boltz/07_run_pae.sh"

## Find all subdirectories named "predictions" within the parent directory
## Note that it will find all subdirectories named "predictions" regardless of their directory depth
find "$PARENT_DIR" -type d -name "predictions" | while IFS= read -r pred_dir; do
  echo "Launching $PAE_SCRIPT for: $pred_dir"
  sbatch "$PAE_SCRIPT" "$pred_dir"
done

