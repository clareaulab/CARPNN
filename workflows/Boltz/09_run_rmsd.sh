#!/bin/sh
#SBATCH --job-name=RMSD
#SBATCH --time=20:00:00
#SBATCH --error=%A.err
#SBATCH --output=%A.out
#SBATCH --mem=16G
#SBATCH --partition=lareauc_cpu,cpu
#SBATCH --output=rmsd_output.out
#SBATCH --error=rmsd_error.err

## User-provided
INPUT=$1  # Csv file containing paths to monomer pdbs / parental pdb
OUTPUT=$2 #EX: /data1/lareauc/users/chuh/ProteinDesign/ActiveCampaigns/BCMABinder/ProteinMPNN_Optimization/boltz_data/tables/rmsd_monomer.csv or rmsd_parental.csv

## Paths
## Use the CARPNN python
CARPNN_PYTHON="/data1/lareauc/users/chuh/miniconda3/envs/carpnn/bin/python"
CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
RMSD_UTIL_SCRIPT="${CARPNN_DIR}/workflows/Boltz/calculate_rmsd.py"

# This blocks calculates the monomer RMSD only
${CARPNN_PYTHON} ${RMSD_UTIL_SCRIPT} \
 -csv ${INPUT} \
 -pred_col "binder_path" \
 -monomer_col "monomer_pdb" \
 --compute_monomer_rmsd \
 --binder_chain "A" \
 --target_chain "B" \
 -out ${OUTPUT} 


#  #This block calculates the full thing (monomer rmsd, rmsd from parental, hotspot rmsd) for (BindCraft + MPNN)
# Assumeing the parental binder has binder on chain A and target on chain B
# ${CARPNN_PYTHON} ${RMSD_UTIL_SCRIPT} \
#  -csv ${INPUT} \
#  -pred_col "binder_path" \
#  -monomer_col "monomer_pdb" \
#  -backbone_col "parental_pdb" \
#  --binder_chain "A" \
#  --target_chain "B" \
#  --all \
#  -out ${OUTPUT} 

# # This block calculates the full thing (monomer rmsd, rmsd from parental, hotspot rmsd) (RFDiffusion)
# ${CARPNN_PYTHON} ${RMSD_UTIL_SCRIPT} \
#  -csv ${INPUT} \
#  -pred_col "boltz1_pred_pdb_path" \
#  -monomer_col "af2_monomer_pdb_path" \
#  -backbone_col "rfd_backbone_path" \
#  --binder_chain "A" \
#  --target_chain "B" \
#  --all \
#  -out ${OUTPUT} 

# # This block calculates the parental RMSD
# ${CARPNN_PYTHON} ${RMSD_UTIL_SCRIPT} \
#  -csv ${INPUT} \
#  -pred_col "binder_path" \
#  -monomer_col "parental_pdb" \
#  --compute_monomer_rmsd \
#  --binder_chain "A" \
#  --target_chain "B" \
#  -out ${OUTPUT}
