#!/bin/bash
INPUT_PARENT_DIR=$1    # e.g. /path/to/inputs
OUTPUT_PARENT_DIR=$2   # e.g. /path/to/outputs

## ---- Paths ----
## A path set below (uncomment and edit its line) or already exported in your environment is never
## overwritten. Anything still empty is filled in from config.sh, written by ./configure.sh.
## Previous hardcoded paths:
## Replace with you path to the CARPNN directory
# CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
CARPNN_CONFIG="${CARPNN_CONFIG:-$HOME/.config/carpnn/config.sh}"
if [ -f "$CARPNN_CONFIG" ]; then . "$CARPNN_CONFIG"; fi
: "${CARPNN_DIR:?not set - hardcode it above or run ./configure.sh}"
BOLTZ_RUN_SCRIPT="${CARPNN_DIR}/workflows/Boltz/05_run_boltz_with_local_msa.sh"

for SUBDIR in "${INPUT_PARENT_DIR}"/*; do
    if [ -d "$SUBDIR" ]; then
        BASENAME=$(basename "$SUBDIR")
        INPUT_PATH="${SUBDIR}"
        OUTPUT_PATH="${OUTPUT_PARENT_DIR}/${BASENAME}"

        # Check if a non-empty 'predictions' directory exists anywhere inside OUTPUT_PATH
        FOUND_NONEMPTY_PRED=$(find "$OUTPUT_PATH" -type d -name "predictions" -exec bash -c '[ "$(find "{}" -mindepth 1 | wc -l)" -gt 0 ]' \; -print -quit)

        if [ -n "$FOUND_NONEMPTY_PRED" ]; then
            echo "⚠️  Skipping ${BASENAME}: found non-empty 'predictions' directory at ${FOUND_NONEMPTY_PRED}"
            continue
        fi

        echo "🚀 Submitting sbatch ${BOLTZ_RUN_SCRIPT} ${INPUT_PATH} ${OUTPUT_PATH}"
        sbatch ${BOLTZ_RUN_SCRIPT} "${INPUT_PATH}" "${OUTPUT_PATH}"
    fi
done