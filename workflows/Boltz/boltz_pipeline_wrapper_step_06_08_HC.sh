#!/bin/bash
#SBATCH --job-name=BoltzPipe
#SBATCH --output=BoltzPipe.out
#SBATCH --error=BoltzPipe.err

## Given prediction outputs from boltz, this wrapper calculates additional metrics such as:
## pae_interaction, 

PREDICTION_OUTPUT_DIRECTORY=$1

## ---- Paths ----
## A path set below (uncomment and edit its line) or already exported in your environment is never
## overwritten. Anything still empty is filled in from config.sh, written by ./configure.sh.
## Previous hardcoded paths:
# CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
CARPNN_CONFIG="${CARPNN_CONFIG:-$HOME/.config/carpnn/config.sh}"
if [ -f "$CARPNN_CONFIG" ]; then . "$CARPNN_CONFIG"; fi
: "${CARPNN_DIR:?not set - hardcode it above or run ./configure.sh}"
SCRIPT_DIR="${CARPNN_DIR}/workflows/Boltz" # Hardcoded, actually runs boltz2

IPSAE_SCRIPT="${SCRIPT_DIR}/06b_run_ipsae_array.sh"
PAE_SCRIPT="${SCRIPT_DIR}/07b_run_pae_array.sh"
ROSETTA_SCRIPT="${SCRIPT_DIR}/08b_run_rosetta_array.sh"

# -- step 6: calculate ipSAE across directories
sbatch ${IPSAE_SCRIPT} ${PREDICTION_OUTPUT_DIRECTORY}

# -- step 7: calculate pae interaction across directories
sbatch ${PAE_SCRIPT} ${PREDICTION_OUTPUT_DIRECTORY}

# -- step 8: Run BindCraft style interface calculation
sbatch ${ROSETTA_SCRIPT} ${PREDICTION_OUTPUT_DIRECTORY} -no_relax
#sbatch ${ROSETTA_SCRIPT} ${PREDICTION_OUTPUT_DIRECTORY}
