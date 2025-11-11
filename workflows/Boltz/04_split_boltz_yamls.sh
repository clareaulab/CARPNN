#!/bin/bash
INPUT_DIR=$1
OUTPUT_DIR=$2
N_CHUNKS=$3

CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
## Replace with the default CARPNN python
PYTHON_PATH="/data1/lareauc/users/chuh/miniconda3/envs/carpnn/bin/python" # Any environment with Biopython will do
SPLIT_DIR_PY="${CARPNN_DIR}/workflows/Boltz/split_dir.py"

${PYTHON_PATH} ${SPLIT_DIR_PY} -src ${INPUT_DIR} -dest ${OUTPUT_DIR} -n ${N_CHUNKS} -csv ${OUTPUT_DIR}/split.csv