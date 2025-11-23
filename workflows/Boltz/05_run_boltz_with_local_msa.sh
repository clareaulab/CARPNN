#!/bin/bash
#SBATCH --job-name=boltz
#SBATCH --output=run_boltz.out
#SBATCH --error=run_boltz.err
#SBATCH --partition=lareauc_gpu,gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=END,FAIL

INPUT_PATH=$1
OUTPUT_DIR=$2

DIFFUSION_SAMPLES=5
RECYCLING_STEPS=10

## Replace with path to the boltz environment
BOLTZ_PATH="/data1/lareauc/users/chuh/miniconda3/envs/boltz2/bin/boltz"

## To perform local MSA, we need to provide additional msa directory
${BOLTZ_PATH} predict ${INPUT_PATH} --out_dir ${OUTPUT_DIR} --recycling_steps ${RECYCLING_STEPS} --diffusion_samples ${DIFFUSION_SAMPLES} --output_format "pdb" --write_full_pae	
