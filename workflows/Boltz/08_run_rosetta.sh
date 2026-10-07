#!/bin/sh
#SBATCH --job-name=BCinterface
#SBATCH --time=20:00:00
#SBATCH --error=%A.err
#SBATCH --output=%A.out
#SBATCH --mem=16G
#SBATCH --partition=lareauc_cpu,cpu
#SBATCH --output=rosetta.out
#SBATCH --error=rosetta.err

## User-provided
INPUT=$1  # Input directory 
BINDER_CHAIN=${2:-"A"}  # Default to "A" if not provided

## ---- Paths ----
## A path set below (uncomment and edit its line) or already exported in your environment is never
## overwritten. Anything still empty is filled in from config.sh, written by ./configure.sh.
## Previous hardcoded paths:
## Replace path with BindCraft python
# BINDCRAFT_PYTHON="/data1/lareauc/users/chuh/miniconda3/envs/BindCraft/bin/python"
# CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
## DSSP / DAlphaBall used by bindcraft_utils.py (both ship in BindCraft's functions/ directory)
# DSSP_PATH="/path/to/BindCraft/functions/dssp"
# DALPHABALL_PATH="/path/to/BindCraft/functions/DAlphaBall.gcc"
CARPNN_CONFIG="${CARPNN_CONFIG:-$HOME/.config/carpnn/config.sh}"
if [ -f "$CARPNN_CONFIG" ]; then . "$CARPNN_CONFIG"; fi
export DSSP_PATH DALPHABALL_PATH   # read by bindcraft_utils.py
: "${CARPNN_DIR:?not set - hardcode it above or run ./configure.sh}"
: "${BINDCRAFT_PYTHON:?not set - hardcode it above or run ./configure.sh}"
: "${DSSP_PATH:?not set - hardcode it above or run ./configure.sh}"
: "${DALPHABALL_PATH:?not set - hardcode it above or run ./configure.sh}"
BINDCRAFT_UTIL_SCRIPT="${CARPNN_DIR}/workflows/Boltz/bindcraft_utils.py"

# Check if the optional "-no_relax" flag is present in arguments 3 and beyond
RELAX_FLAG=""
for arg in "$@"; do
    if [ "$arg" = "-no_relax" ]; then
        RELAX_FLAG="-no_relax"
        break
    fi
done

## Runs BindCraft-based Rosetta Interface Scoring Outputs (bindcraft.util) for each subdirectory in INPUT
for SUBDIR in ${INPUT}/*/; do
    # if [ -f "${SUBDIR}/binding_interface.tsv" ]; then
    #     echo "Skipping ${SUBDIR} (binding_interface.tsv already exists)"
    #     continue
    # fi
    
    ${BINDCRAFT_PYTHON} ${BINDCRAFT_UTIL_SCRIPT} \
     -pdbdir ${SUBDIR} \
     -binder_chain ${BINDER_CHAIN} \
     -out ${SUBDIR}/binding_interface.tsv \
     ${RELAX_FLAG}
done


