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
COLABFOLD_CSV=$3 

CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
## Replace with the default CARPNN python
PYTHON_PATH="/data1/lareauc/users/chuh/miniconda3/envs/carpnn/bin/python" # Any environment with Biopython will do
BOLTZ_YAML_PY="${CARPNN_DIR}/workflows/Boltz/generate_yaml_multimer.py"

mkdir -p ${OUT_DIR}

# Construct the base command
# Given msa outdir, create a directory of yaml files with properly mapped inputs
CMD="${PYTHON_PATH} ${BOLTZ_YAML_PY} ${A3M_DIR} -o ${OUT_DIR}"
# Append colabfold flag if the third argument is provided
if [ -n "${COLABFOLD_CSV}" ]; then
    echo "Running in ColabFold mode with metadata: ${COLABFOLD_CSV}"
    CMD="${CMD} -colabfold ${COLABFOLD_CSV}"
else
    echo "Running in standard mode (no renaming)."
fi
## Execute the constructed command
echo "Executing: ${CMD}"
${CMD}