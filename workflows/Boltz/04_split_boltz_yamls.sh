#!/bin/bash
INPUT_DIR=$1
OUTPUT_DIR=$2
N_CHUNKS=$3

## ---- Paths ----
## A path set below (uncomment and edit its line) or already exported in your environment is never
## overwritten. Anything still empty is filled in from config.sh, written by ./configure.sh.
## Previous hardcoded paths:
# CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
## Replace with the default CARPNN python
# PYTHON_PATH="/data1/lareauc/users/chuh/miniconda3/envs/carpnn/bin/python" # Any environment with Biopython will do
CARPNN_CONFIG="${CARPNN_CONFIG:-$HOME/.config/carpnn/config.sh}"
if [ -f "$CARPNN_CONFIG" ]; then . "$CARPNN_CONFIG"; fi
PYTHON_PATH="${PYTHON_PATH:-$CARPNN_PYTHON}"   # a PYTHON_PATH set above wins
: "${CARPNN_DIR:?not set - hardcode it above or run ./configure.sh}"
: "${PYTHON_PATH:?not set - hardcode it above or run ./configure.sh}"
SPLIT_DIR_PY="${CARPNN_DIR}/workflows/Boltz/split_dir.py"

${PYTHON_PATH} ${SPLIT_DIR_PY} -src ${INPUT_DIR} -dest ${OUTPUT_DIR} -n ${N_CHUNKS} -csv ${OUTPUT_DIR}/split.csv