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
OUT_DIR=$2

CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
## Replace with the default CARPNN python
PYTHON_PATH="/data1/lareauc/users/chuh/miniconda3/envs/carpnn/bin/python" # Any environment with Biopython will do
BOLTZ_YAML_PY="${CARPNN_DIR}/workflows/Boltz/generate_yaml_multimer.py"

mkdir -p ${OUT_DIR}

## Given msa outdir, create a directory of yaml files with properly mapped inputs
${PYTHON_PATH} ${BOLTZ_YAML_PY} ${A3M_DIR} -o ${OUT_DIR}
