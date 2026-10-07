#!/bin/sh
PARENT_DIR=$1

## ---- Paths ----
## A path set below (uncomment and edit its line) or already exported in your environment is never
## overwritten. Anything still empty is filled in from config.sh, written by ./configure.sh.
## Previous hardcoded paths:
# CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
CARPNN_CONFIG="${CARPNN_CONFIG:-$HOME/.config/carpnn/config.sh}"
if [ -f "$CARPNN_CONFIG" ]; then . "$CARPNN_CONFIG"; fi
: "${CARPNN_DIR:?not set - hardcode it above or run ./configure.sh}"
PAE_SCRIPT="${CARPNN_DIR}/workflows/Boltz/07_run_pae.sh"

## Find all subdirectories named "predictions" within the parent directory
## Note that it will find all subdirectories named "predictions" regardless of their directory depth
find "$PARENT_DIR" -type d -name "predictions" | while IFS= read -r pred_dir; do
  echo "Launching $PAE_SCRIPT for: $pred_dir"
  sbatch "$PAE_SCRIPT" "$pred_dir"
done

