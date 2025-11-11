#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --partition=lareauc_cpu
#SBATCH --time 5:00:00
#SBATCH --output=LigandMPNN_%A.log
#SBATCH --output=LigandMPNN_%A.err

INPUT_PDB=$1      # Input PDB file
INTERFACE_RESIDUES=$2 # Path to residues to redesign
SUFFIX_NAME=$3 # Usuall "_interface_redesigned" or "non_interface_redesigned"
OUTPUT_DIR=$4     # Output directory

BINDER_CHAIN="A"
TARGET_CHAIN="B"

## Total number of sequences = BATCH_SIZE * NUM_BATCH
BATCH_SIZE=5
NUM_BATCH=1
TEMPERATURE=0.4 # For large number of sequences we recommend increasing the temperature

## Replace this with the directory and python to your own LigandMPNN environment
LIGANDMPNN_DIR="/data1/lareauc/users/chuh/softwares/LigandMPNN"
LIGANDMPNN_PYTHON="/data1/lareauc/users/chuh/miniconda3/envs/ligandmpnn_env/bin/python"
## Replace this with the model param to your own LigandMPNN environment
CHECKPOINT_PATH="/data1/lareauc/users/chuh/softwares/LigandMPNN/model_params/solublempnn_v_48_020.pt"

# Create output and input directories
mkdir -p "${OUTPUT_DIR}/input"

# Copy this script to the output directory for reproducibility
cp "$0" "${OUTPUT_DIR}/$(basename "$0")"

# Copy input files
cp "${INPUT_PDB}" "${OUTPUT_DIR}/input/$(basename "${INPUT_PDB}")"
cp "${INTERFACE_RESIDUES}" "${OUTPUT_DIR}/input/$(basename "${INTERFACE_RESIDUES}")"

# Log the arguments to a file
{
  echo "INPUT_PDB=${INPUT_PDB}"
  echo "INTERFACE_RESIDUES=${INTERFACE_RESIDUES}"
  echo "OUTPUT_DIR=${OUTPUT_DIR}"
  echo "BINDER_CHAIN=${BINDER_CHAIN}"
  echo "TARGET_CHAIN=${TARGET_CHAIN}"
  echo "SUFFIX_NAME=${SUFFIX_NAME}"
} > "${OUTPUT_DIR}/run_arguments.txt"

# Run SoluableMPNN with interface residues fixed
${LIGANDMPNN_PYTHON} ${LIGANDMPNN_DIR}/run.py \
    --seed 1 \
    --model_type "soluble_mpnn" \
    --checkpoint_soluble_mpnn ${CHECKPOINT_PATH} \
    --pdb_path ${INPUT_PDB} \
    --out_folder ${OUTPUT_DIR} \
    --fixed_residues_multi ${INTERFACE_RESIDUES} \
    --temperature ${TEMPERATURE} \
    --save_stats 1 \
    --omit_AA "C" \
    --chains_to_design "${BINDER_CHAIN}" \
    --batch_size ${BATCH_SIZE} \
    --number_of_batches ${NUM_BATCH} \
    --file_ending ${SUFFIX_NAME}

