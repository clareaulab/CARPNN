## MPNN-related workflow (ProteinMPNN/SoluableMPNN)

This directory contains workflows to run SoluableMPNN (LigandMPNN is just the name of the repository)

To setup, follow the instruction from the LigandMPNN repo and replace the paths indicated in the scripts accordingly

## Script setup and Descriptions

### All scripts

The sbatch headers configurations (ex: `#SBATCH --partition=cpu`) should be updated to fit your HPC system. 

### `swap_chain.py`
Utility script. Use it to swap chain A and chain B of a pdb structure.

Example: `python swap_chain.py examples/02_lead_candidate_mutagenized/input_structure/BCMA_l54_s662851_mpnn3_model2.pdb examples/02_lead_candidate_mutagenized/input_structure/BCMA_l54_s662851_mpnn3_model2_chain_swapped.pdb`

### `get_interface_residues.py`
Utility script. Use it to obtain a json listing interface residues in a given pdb file.

Example: `python get_interface_residues.py /data1/lareauc/users/chuh/softwares/CARPNN/examples/02_lead_candidate_mutagenized/input_structure/BCMA_l54_s662851_mpnn3_model2_chain_swapped.pdb A B --cutoff 4 --output_json /data1/lareauc/users/chuh/softwares/CARPNN/examples/02_lead_candidate_mutagenized/input_structure/BCMA_l54_interface_residues.jsonl`

### `SoluableMPNN_noninterface_mutagenesis.sh`
This script runs SoluableMPNN on the non-interface residues of the provided structure based on the interface residues json provided (can be obtained using `get_interface_residues.py`). The total number of sequence produced = BATCH_SIZE * NUM_BATCH so modify it to suit your needs.

- The `LIGANDMPNN_DIR` variable should be updated to your installation directory of ligandMPNN
- The `LIGANDMPNN_PYTHON` variable should be updated to the path to python in your LigandMPNN conda environments
- The `CHECKPOINT_PATH` variable should be updated to your downloaded weights of LigandMPNN

### `SoluableMPNN_interface_mutagenesis.sh`
This script runs SoluableMPNN on the interface residues of the provided structure based on the interface residues json provided (can be obtained using `get_interface_residues.py`). The total number of sequence produced = BATCH_SIZE * NUM_BATCH so modify it to suit your needs.

- The `LIGANDMPNN_DIR` variable should be updated to your installation directory of ligandMPNN
- The `LIGANDMPNN_PYTHON` variable should be updated to the path to python in your LigandMPNN conda environments
- The `CHECKPOINT_PATH` variable should be updated to your downloaded weights of LigandMPNN

### `SoluableMPNN_all_mutagenesis.sh`
This script runs SoluableMPNN on all residues of the provided structure. Note that the interface json argument is kept for consistency with other scripts but is not used to actually run the tool. The total number of sequence produced = BATCH_SIZE * NUM_BATCH so modify it to suit your needs.

- The `LIGANDMPNN_DIR` variable should be updated to your installation directory of ligandMPNN
- The `LIGANDMPNN_PYTHON` variable should be updated to the path to python in your LigandMPNN conda environments
- The `CHECKPOINT_PATH` variable should be updated to your downloaded weights of LigandMPNN

## Others

If you use this code, please cite the orginal LigandMPNN tool

https://github.com/dauparas/LigandMPNN
