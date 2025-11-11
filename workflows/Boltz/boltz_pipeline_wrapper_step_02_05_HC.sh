#!/bin/bash
#SBATCH --job-name=BoltzPipe
#SBATCH --output=BoltzPipe.out
#SBATCH --error=BoltzPipe.err

# ==============================================================================
# Integrated Protein Interaction Prediction and Analysis Pipeline
#
# This script orchestrates a multi-step workflow for large-scale protein
# interaction prediction using the Boltz tool on an HPC cluster. It chains
# together several scripts using SBATCH dependencies to ensure correct order
# of operations.
#
# Usage:
#   sbatch this_script.sh <path_to_input_csv> <master_output_directory> <number_of_chunks>
#
# Arguments:
#   $1: Path to the input CSV file for ColabSearch.
#   $2: Master output directory where all results will be stored.
#   $3: Number of chunks to split the YAML files into for parallel processing.
# ==============================================================================

# --- User-defined variables ---

# --- Pipeline inputs from command line ---
INPUT_CSV=$1
MASTER_OUTPUT_DIR=$2
N_CHUNKS=$3

# Replace this by your own CARPNN project directory
CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"

# These are the paths to the other pipeline scripts. Assumes they are in the same directory.
SCRIPT_DIR="${CARPNN_DIR}/workflows/Boltz" 

GENERATE_MSA_SCRIPT="${SCRIPT_DIR}/02_generate_msa.sh"
PREP_YAMLS_SCRIPT="${SCRIPT_DIR}/03_prep_boltz_yamls.sh"
SPLIT_YAMLS_SCRIPT="${SCRIPT_DIR}/04_split_boltz_yamls.sh"
RUN_BOLTZ_ARRAY_SCRIPT="${SCRIPT_DIR}/05b_run_boltz_with_local_msa_array.sh"
RUN_IPSAE_ARRAY_SCRIPT="${SCRIPT_DIR}/06b_run_ipsae_array.sh"
RUN_PAE_ARRAY_SCRIPT="${SCRIPT_DIR}/07b_run_pae_array.sh"

# --- Derived directories ---
MSA_OUTPUT_DIR="${MASTER_OUTPUT_DIR}/colabsearch_output"
YAML_OUTPUT_DIR="${MASTER_OUTPUT_DIR}/boltz_inputs"
SPLIT_YAML_OUTPUT_DIR="${YAML_OUTPUT_DIR}/input_yamls_chunked"
PREDICTION_OUTPUT_DIR="${MASTER_OUTPUT_DIR}/boltz_predictions"

# --- Step 1: Generate MSAs using ColabSearch ---
echo "--- Step 1: Submitting MSA generation job..."
JOB_ID_MSA=$(sbatch --parsable "${GENERATE_MSA_SCRIPT}" "${INPUT_CSV}" "${MSA_OUTPUT_DIR}")
echo "MSA generation job ID: ${JOB_ID_MSA}"

# --- Step 2: Prepare YAML files for Boltz ---
# This step depends on the successful completion of the MSA generation job.
echo "--- Step 2: Submitting YAML preparation job..."
JOB_ID_YAML=$(sbatch --parsable --dependency=afterok:${JOB_ID_MSA} \
    "${PREP_YAMLS_SCRIPT}" "${MSA_OUTPUT_DIR}" "${INPUT_CSV}" "${YAML_OUTPUT_DIR}")
echo "YAML preparation job ID: ${JOB_ID_YAML}"

# --- Step 3: Split YAML files into chunks for array job ---
# This step depends on the successful completion of the YAML preparation.
echo "--- Step 3: Submitting YAML splitting job..."
JOB_ID_SPLIT=$(sbatch --parsable --dependency=afterok:${JOB_ID_YAML} \
    "${SPLIT_YAMLS_SCRIPT}" "${YAML_OUTPUT_DIR}/input_yamls" "${SPLIT_YAML_OUTPUT_DIR}" "${N_CHUNKS}")
echo "YAML splitting job ID: ${JOB_ID_SPLIT}"

# --- Step 4: Run Boltz predictions in parallel via an array job ---
# This job depends on the successful splitting of YAML files.
echo "--- Step 4: Submitting Boltz prediction array job..."
# Note: --wait ensures this script pauses until the array job is complete.
JOB_ID_PREDICTION=$(sbatch --wait --dependency=afterok:${JOB_ID_SPLIT} \
    "${RUN_BOLTZ_ARRAY_SCRIPT}" "${SPLIT_YAML_OUTPUT_DIR}" "${PREDICTION_OUTPUT_DIR}")
echo "Boltz prediction job ID: ${JOB_ID_PREDICTION}"
