## AlphaFold Single-sequence Mode Refolding Workflow
This runs AlphaFold in single-sequence mode over large batch of sequences.

## Script setup and Descriptions

To set up the workflow, run:
`bash set_up_AlphaFold_workflow.sh`
which will download the heme_binder_diffusion repository from et al. (removing the LigandMPNN submodule) and place it repository under the `public/` repository as well as downloading the alphafold weights.

You would also need to create a `mlfold` conda environment using the yaml file provided in `envs/mlfold` and replace the scripts python path to the one used by the environment.

`run_alphafold_monomer.sh`
This scripts will fold up all sequences in a fasta file as monomers and write out the lDDT scores as well as the pdb files in the output directory.

`run_alphafold_monomer_array.sh`
This is a convenience script that will run `run_alphafold_monomer.sh` over a directory containing a list of fastas, lauch one job per fasta.

If you use the code here, please also cite the study where the code based on:
heme diffusion

