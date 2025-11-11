#!/bin/sh
#SBATCH --job-name=AFMonomer
#SBATCH --time=20:00:00
#SBATCH --error=AF2_out_%A.err
#SBATCH --output=AF2_out_%A.out
#SBATCH --mem=32G
#SBATCH --partition=lareauc_gpu,gpu
#SBATCH --gres=gpu:1

################################ Bash set up ##################################

# I/O Arguments
INPUT_FASTA=$1
OUTPUT_DIR=$2

## Environment Configurations
## Replace with the path to this repo
CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"

## Replace the RFDAA_AF2_PYTHON with the python in the mlfold environment
RFDAA_AF2_PYTHON="/data1/lareauc/users/chuh/miniconda3/envs/mlfold/bin/python"
RFDAA_AF2_PY="${CARPNN_DIR}/public/heme_binder_diffusion/scripts/af2/af2.py"

# Input parameters
AF2_MODELS="4"
AF2_RECYCLES=10

# Run AF2 prediction
mkdir -p ${OUTPUT_DIR}
(
  cd ${OUTPUT_DIR}/;
  ${RFDAA_AF2_PYTHON} ${RFDAA_AF2_PY} \
    --af-nrecycles ${AF2_RECYCLES} \
    --af-models ${AF2_MODELS} \
    --fasta ${INPUT_FASTA} \
    --scorefile ${OUTPUT_DIR}/AF2_score.csv
)


