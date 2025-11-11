#!/bin/sh
#SBATCH --job-name=ColabSearch
#SBATCH --time=40:00:00
#SBATCH --mem=512G
#SBATCH --partition=lareauc_cpu,cpu
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --output=colabsearch.out
#SBATCH --error=colabsearch.err

## User-provided
INPUT=$1  # Input file (colabsearch-like csv)
OUTDIR=$2  # Output directory for all outputs

## Replace with you path to the CARPNN directory
CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
## Replace with you path to ColabFold singularity container
SIF_PATH="/data1/lareauc/shared_resources/colabfold_sif/colabfold_1.5.5-cuda12.2.2.sif"
## Replace with you path to the ColabFold database
DB_FOLDER="/data1/lareauc/shared_resources/colabfold_db"

INPUT_DIR=$(dirname ${INPUT})
INPUT_BASENAME=$(basename ${INPUT})

# Create output directory for this job
mkdir -p ${OUTDIR}

# Run singularity container
singularity run \
    -B ${INPUT_DIR}:/input_dir \
    -B $(pwd):/work \
    -B ${OUTDIR}:/outputs \
    -B ${DB_FOLDER}:/db_folder \
    ${SIF_PATH} \
    colabfold_search /input_dir/${INPUT_BASENAME} /db_folder /outputs






