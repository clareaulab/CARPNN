#!/bin/bash

# Usage: ./run_af2_wrapper.sh /path/to/input_fasta_dir /path/to/output_dir
INPUT_DIR=$1
OUTPUT_DIR=$2

# Path to SLURM script (assumed to be in the same directory as this wrapper)
## ---- Paths ----
## A path set below (uncomment and edit its line) or already exported in your environment is never
## overwritten. Anything still empty is filled in from config.sh, written by ./configure.sh.
## Previous hardcoded paths:
## Replace with the path to this repo
# CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
CARPNN_CONFIG="${CARPNN_CONFIG:-$HOME/.config/carpnn/config.sh}"
if [ -f "$CARPNN_CONFIG" ]; then . "$CARPNN_CONFIG"; fi
: "${CARPNN_DIR:?not set - hardcode it above or run ./configure.sh}"
AF2_SLURM_SCRIPT="${CARPNN_DIR}/workflows/AlphaFold/run_alphafold_monomer.sh"

# Loop through all .fasta files in the input directory
for fasta_file in "${INPUT_DIR}"/*.fasta; do
  [ -e "$fasta_file" ] || continue  # Skip if no FASTA files
  fasta_base=$(basename "$fasta_file" .fasta)
  output_subdir="${OUTPUT_DIR}/${fasta_base}"

  echo "Submitting job for $fasta_file -> $output_subdir"
  mkdir -p "$output_subdir"
  sbatch "$AF2_SLURM_SCRIPT" "$fasta_file" "$output_subdir"
done