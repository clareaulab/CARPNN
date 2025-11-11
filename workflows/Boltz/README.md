## Large-scale Boltz Refolding Workflow
Warning: Running Boltz locally require downloading the entirely of colabfold's database which is ~2TB in size and could take substantial amount of time

## tl;dr
1. The `boltz_pipeline_wrapper_step_02_05_HC.sh` script takes in colabfold-style csv file and runs step 02-05 to predict all of the structures.
2. The `boltz_pipeline_wrapper_step_06_08_HC.sh` script runs step 06-08 to extract additional metrics from the outputs from the previous step.
3. The `10_aggregate_boltz_outputs.sh` script aggregates scores from previous outputs for downstream analysis

For example of what the output structure should look like, see

## Script setup and Descriptions

### All scripts

The `CARPNN_DIR` variable should be updated to the path to this CARPNN installation directory.

The sbatch headers configurations (ex: `#SBATCH --partition=cpu`) should be updated to fit your HPC system. 

### `02_generate_msa.sh`

This scripts takes in a colabfold-style csv and performs multi-sequence alignment (paired and unpaired) against a local database.
Follow the instruction from https://github.com/sokrypton/ColabFold/wiki/Running-ColabFold-in-Docker to set up the singularity container and the instructions from https://colabfold.mmseqs.com/ to set up the colabfold database. In particular, you should update these two variables.

- The `SIF_PATH` variable should be updated to the path to the singularity container downloaded
  - `singularity pull docker://ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2`
- The `DB_FOLDER` variable should be updated to the parent directory where the databases are hosted
  - `wget https://raw.githubusercontent.com/sokrypton/ColabFold/main/setup_databases.sh`
  - `chmod +x setup_databases.sh`
  - `./setup_databases.sh database/`

Currently the input csv must have two columns named `id` and `sequence`. The id must be formatted by the binder name and target name delimitted by a double underscored. For example: `binder_1__BCMA`. The sequence must contain exactly two sequence delimited by `:`. For example: `BINDERSEQ:TARGETSEQ`

Setting up a local MSA search also require massive amount of memory so make sure your sbatch configuration is very high in memory (at least 512GB) for efficient search.

### `03_prep_boltz_yamls`

This scripts takes in the colabsearch alignment outputs (a3m files) from `02_generate_msa.sh` and converts them into yaml files that can be processed by Boltz.

- The `CARPNN_DIR` variable should be updated to the path to this CARPNN installation directory
- The `PYTHON_PATH` should be updated to the CAR-PNN environment python

### `04_split_boltz_yamls`

This scripts takes in the yaml file generated from `03_prep_boltz_yamls.sh` and split them into N chunks of directories for downstream parallelized processing

- The `CARPNN_DIR` variable should be updated to the path to this CARPNN installation directory
- The `PYTHON_PATH` should be updated to the CAR-PNN environment python

### `05_run_boltz_with_local_msa.sh`

This is the main script that calls Boltz to predict structures. By default the number diffusion samples is set to 5 and recycle is set to 10.

- The `CARPNN_DIR` variable should be updated to the path to this CARPNN installation directory
- The `BOLTZ_PATH` variable should be updated to the Boltz executable in the boltz environment

### `05b_run_boltz_with_local_msa_array.sh`

A convenience script that runs `05_run_boltz_with_local_msa.sh` over a list of directories. In the event that some of the predictions failed (ex: due to CUDA ECC error) one can relaunch this script over the same directories and it will relaunch for those where the output directory is empty.

- The `CARPNN_DIR` variable should be updated to the path to this CARPNN installation directory

### `06_run_ipsae.sh`

This script calculates ipSAE outputs based on the provided boltz prediction directory.

- The `CARPNN_DIR` variable should be updated to the path to this CARPNN installation directory

### `06b_run_ipsae_array.sh`

A convenience script that runs `06_run_ipsae.sh` over a list of directories. 

- The `CARPNN_DIR` variable should be updated to the path to this CARPNN installation directory
- The `PYTHON_PATH` should be updated to the CAR-PNN environment python

### `07_run_pae.sh`

This script calculates pAE interaction based on the provided boltz prediction directory.

- The `CARPNN_DIR` variable should be updated to the path to this CARPNN installation directory

### `07b_run_pae_array.sh`

A convenience script that runs `07_run_pae.sh` over a list of directories. 

- The `CARPNN_DIR` variable should be updated to the path to this CARPNN installation directory

### `08_run_rosetta.sh`

This script calculates the same set of biophyisical metrics calculated by BindCraft on the provided boltz prediction directory. Following the instruction at https://github.com/martinpacesa/BindCraft/tree/main to install BindCraft. Alternatively, you could install scipy and pyrosetta and use the CARPNN python environment instead

- The `BINDCRAFT_PYTHON` variable should be updated to the python used in the BindCraft environment
- The `CARPNN_DIR` variable should be updated to the path to this CARPNN installation directory

Note that this sccipt can also be ran with an addition flag `-no_relax`. If specified, this will skip the relaxation step of the interface biophyisical metric calculation which could be very slow for larger complexes. However this will result in predicted stuctures with more clashes which will impact dG calculations and the filters used downtream should be adjusted accordingly.

### `08b_run_rosetta_array.sh`

A convenience script that runs `08_run_rosetta.sh` over a list of directories. The `-no_relax` flag can also be specified here.

- The `CARPNN_DIR` variable should be updated to the path to this CARPNN installation directory

### `09_run_rmsd.sh`

This script compares the RMSD between two structures for all rows listed in a csv. Note that this table is not created by any of the previous steps in this workflow and will need to be created separately. More description on how this table is created cna be found in the example notebook.

- The `CARPNN_DIR` variable should be updated to the path to this CARPNN installation directory
- The `CARPNN_PYTHON` should be updated to the CAR-PNN environment python

### `10_aggregate_boltz_outputs.sh`

This script aggregates the predictions files in each directory into a csv that can be easily processed by downstream analysis.

- The `CARPNN_DIR` variable should be updated to the path to this CARPNN installation directory
- The `CARPNN_PYTHON` should be updated to the CAR-PNN environment python


If you use the code here. Please also cite the studies where the code is based on:
- Boltz1/2
- BindCraft
- ColabFold