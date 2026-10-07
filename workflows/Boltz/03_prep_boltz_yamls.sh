#!/bin/bash
#SBATCH --job-name=yaml_prep
#SBATCH --output=yaml_prep.log
#SBATCH --error=yaml_prep.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --partition=cpu

## Given colabseatch MSA directory and colabsearch input, create a directory of yaml files with properly mapped inputs
A3M_DIR=$1
COLABSEARCH_INPUT=$2
OUT_DIR=$3

## ---- Paths ----
## A path set below (uncomment and edit its line) or already exported in your environment is never
## overwritten. Anything still empty is filled in from config.sh, written by ./configure.sh.
## Previous hardcoded paths:
# CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
## Replace with the default CARPNN python
# PYTHON_PATH="/data1/lareauc/users/chuh/miniconda3/envs/carpnn/bin/python" # Any environment with Biopython will do
CARPNN_CONFIG="${CARPNN_CONFIG:-$HOME/.config/carpnn/config.sh}"
if [ -f "$CARPNN_CONFIG" ]; then . "$CARPNN_CONFIG"; fi
PYTHON_PATH="${PYTHON_PATH:-$CARPNN_PYTHON}"   # a PYTHON_PATH set above wins
: "${CARPNN_DIR:?not set - hardcode it above or run ./configure.sh}"
: "${PYTHON_PATH:?not set - hardcode it above or run ./configure.sh}"
BOLTZ_YAML_PY="${CARPNN_DIR}/workflows/Boltz/generate_yaml.py"

mkdir -p ${OUT_DIR}

## Given colabfold input, msa outdir, create a directory of yaml files with properly mapped inputs
${PYTHON_PATH} ${BOLTZ_YAML_PY} -i ${COLABSEARCH_INPUT} -msa ${A3M_DIR} -out ${OUT_DIR}
