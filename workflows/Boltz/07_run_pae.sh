#!/bin/sh
#SBATCH --job-name=PAE
#SBATCH --time=1:00:00
#SBATCH --error=%A.err
#SBATCH --output=%A.out
#SBATCH --mem=4G
#SBATCH --partition=lareauc_cpu,cpu
#SBATCH --output=pae_output.out
#SBATCH --error=pae_error.err

## User-provided
INPUT=$1  # Input directory 

## Paths
## Replace with the python used by the CAR-PNN
CARPNN_PYTHON="/data1/lareauc/users/chuh/miniconda3/envs/carpnn/bin/python"
CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
PAE_UTIL_SCRIPT="${CARPNN_DIR}/workflows/Boltz/calculate_pae.py"

# Loop through all subdirectories in the specified directory
for sub_dir in "$INPUT"/*/; do
    echo "doing ${sub_dir}"
    ${CARPNN_PYTHON} ${PAE_UTIL_SCRIPT} \
    -dir ${sub_dir} 
done



