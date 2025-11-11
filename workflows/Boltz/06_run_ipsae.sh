#!/bin/bash
#SBATCH --job-name=ipsae
#SBATCH --output=ipsae.log
#SBATCH --error=ipsae.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --partition=lareauc_cpu,cpu

# Set variables with defaults
DIR=$1
PAE_CUTOFF=${2:-10}       # Default to 10 if not provided
DIST_CUTOFF=${3:-10}      # Default to 10 if not provided
MODEL_TYPE=${4:-boltz1}    # Default to 'boltz1' if not provided
OUTPUT_TYPE=${5:-pdb}    # Default to 'pdb' if not provided

## Replace with your CAR-PNN directory
CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
## Replace with CAR-PNN python
PYTHON_PATH="/data1/lareauc/users/chuh/miniconda3/envs/carpnn/bin/python" # Any python with numpy will do

IPSAE_SCIRPT_PATH="${CARPNN_DIR}/public/ipSAE/ipsae_modified.py"

# Loop through all subdirectories in the specified directory
for sub_dir in "$DIR"/*/; do
    # Loop through all PDB files in the current subdirectory
    for pdb_file in "$sub_dir"*.pdb; do
        # Extract the basename of the PDB file and prepend 'pae_' while changing the extension to .npz
        npz_file="$sub_dir/pae_$(basename "${pdb_file%.pdb}.npz")"
        # Check if the corresponding PAE file exists
        if [[ -f "$npz_file" ]]; then
            echo "Processing $npz_file and $pdb_file with pae_cutoff=$PAE_CUTOFF and dist_cutoff=$DIST_CUTOFF"
            "$PYTHON_PATH" "$IPSAE_SCIRPT_PATH" "$npz_file" "$pdb_file" "$PAE_CUTOFF" "$DIST_CUTOFF" "$MODEL_TYPE" "$OUTPUT_TYPE"
        else
            echo "No corresponding PAE file for $pdb_file in $sub_dir. Skipping..."
        fi
    done
done


