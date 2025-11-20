#!/bin/sh
PARENT_DIR=$1
BINDER_CHAIN=${2:-"A"} # Default to "A" if not provided
# Check if the optional "-no_relax" flag is present in arguments 3 and beyond
RELAX_FLAG=""
for arg in "$@"; do
    if [ "$arg" = "-no_relax" ]; then
        RELAX_FLAG="-no_relax"
        break
    fi
done

CARPNN_DIR="/data1/lareauc/users/chuh/softwares/CARPNN"
ROSETTA_SCRIPT="${CARPNN_DIR}/workflows/Boltz/08_run_rosetta.sh"

## Find all subdirectories named "predictions" within the parent directory
## Note that it will find all subdirectories named "predictions" regardless of their directory depth
find "$PARENT_DIR" -type d -name "predictions" | while IFS= read -r pred_dir; do
  echo "Launching $ROSETTA_SCRIPT for: $pred_dir"
  sbatch "$ROSETTA_SCRIPT" "$pred_dir" ${BINDER_CHAIN} ${RELAX_FLAG}
done
