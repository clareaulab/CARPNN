## AlphaFold Single-sequence Mode Refolding Workflow

This directory contains the workflow to run AlphaFold in single-sequence mode over large batch of monomer sequences. It is optimized for folding up sequences with identical length (for example MPNN draws)

## Workflow setup

To set up the workflow, run:

`bash set_up_AlphaFold_workflow.sh`

This will download the [heme_binder_diffusion](https://github.com/ikalvet/heme_binder_diffusion) repository from the [RFDiffusion-AllAtom](https://github.com/baker-laboratory/rf_diffusion_all_atom) paper (Krishna R et al. Science 2024) (with the LigandMPNN submodule removed) and place it repository under the `public/` repository. It will also download the alphafold weights.

After the repository has been downloaded, you would also need to create a `mlfold` conda environment using the yaml file provided in `public/heme_binder_diffusion/envs/mlfold` 

`conda env create -f public/heme_binder_diffusion/envs/mlfold.yml`

Once the environment has been created, note the python path to the environment 
```
conda activate mlfold
which python # <- remember this path
```

## Script set up & descriptions

### All scripts

The sbatch headers configurations (ex: `#SBATCH --partition=cpu`) should be updated to fit your HPC system. 

### `run_alphafold_monomer.sh`
This scripts takes in a fasta file of interest and will fold up all sequences in a fasta file as monomers and write out the lDDT scores as well as the pdb files in the output directory. Example commands can be found both in the `notebooks/bindcraft_to_boltz` and `notebooks/lead_optimization.ipynb` notebooks.

- The `CARPNN_DIR` variable should be updated to the path to this CARPNN installation directory.
- The `RFDAA_AF2_PYTHON` variable should be updated to the python path used by the `mlfold` environment

By default we run model 4 for 10 recycles. Adjust these parameters as you see fit.

### `run_alphafold_monomer_array.sh`
This is a convenience script that will run `run_alphafold_monomer.sh` over a directory containing a list of fastas, launching one job per fasta.

## Others
If you use the code here, please also cite the study where the code based on:

https://github.com/ikalvet/heme_binder_diffusion

