#!/bin/bash
#SBATCH --job-name=BoltzPipe
#SBATCH --output=BoltzPipe.out
#SBATCH --error=BoltzPipe.err

## Given prediction outputs from boltz, this wrapper calculates additional metrics such as:
## pae_interaction, 

PREDICTION_OUTPUT_DIRECTORY=$1

CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
SCRIPT_DIR="${CARPNN_DIR}/workflows/Boltz" # Hardcoded, actually runs boltz2

IPSAE_SCRIPT="${SCRIPT_DIR}/06b_run_ipsae_array.sh"
PAE_SCRIPT="${SCRIPT_DIR}/07b_run_pae_array.sh"
ROSETTA_SCRIPT="${SCRIPT_DIR}/08b_run_rosetta_array.sh"

# -- step 6: calculate ipSAE across directories
sbatch ${IPSAE_SCRIPT} ${PREDICTION_OUTPUT_DIRECTORY}

# -- step 7: calculate pae interaction across directories
sbatch ${PAE_SCRIPT} ${PREDICTION_OUTPUT_DIRECTORY}

# -- step 8: Run BindCraft style interface calculation
sbatch ${ROSETTA_SCRIPT} ${PREDICTION_OUTPUT_DIRECTORY}
