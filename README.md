# CAR-PNN: Chimeric-Antigen-Receptor ProteinMPNN 

![project_logo](assets/logo.png)

The **CAR-PNN (Chimeric-Antigen-Receptor ProteinMPNN)** workflow is a general purposed, HPC-ready workflow that allows you to perfom various protein design / protein interaction tasks such as
- Large-scale local folding of protein sequences via Boltz
- Implement custom binder filtering configuration similar to the one used by BindCraft
- Lead optimization of de novo binder via SoluableMPNN for downstream applications such as CAR

## Installation

### Step 1: Install environment required for running scripts in this tool

```
conda env create -f carpnn.yml
conda activate carpnn

## Install ipykernel for jupyter notebook
python -m ipykernel install --user --name=carpnn

## Keep a note of the path printed out by this command, this path will be needed to set up some of the following workflows
which python
```
The environment for the tool itself is relatively lightweight. See Step 2 for setting up individual workflows.

### Step 2: Install environments required for individual workflows

If you already have these tools installed else where in your HPC environment you can just replace the paths to their corresponding installation directory without reinstalling. Some of the github repository required include:

- LigandMPNN (For running Protein/SoluableMPNN, https://github.com/dauparas/LigandMPNN)
- Boltz2 (For structure prediction and filtering, https://github.com/jwohlwend/boltz)
- ColabFold singularity container (For using ColabSearch module to perform MSA using local database, https://github.com/YoshitakaMo/localcolabfold)
- BindCraft (For calculation of interface metrics, https://github.com/martinpacesa/BindCraft)
- ipSAE (Already provided in this repo, for interaction scoring, https://github.com/DunbrackLab/IPSAE/tree/main)

For details, read about how to set up the individual workflows in the README files in the individual workflow directorys under `workflows/`. The workflows are decoupled into individual scripts such that if you already have another workflow you are using to create a certain component (ex: binder generation / MSA) they should be easily swappable to your needs.

## How to Run

In short, the binder sequences to lead optimization workflow consist of 4 major steps:

- **Step 1**: Refold candidate binders via Boltz against target antigen(ex: binder sequences from tools like RFDiffusion, BindCraft, BoltzGen, etc.)
- **Step 2**: Identify lead de novo binder (experimentally validated or have high ipSAE score > 0.85)
- **Step 3**: Mutagensis of (non-interface/interface) residues of the lead binder and perform sequence level filter (ex: charge / immunogenicity)
- **Step 4**: Refold mutagenized sequences via Boltz and identify top binders for experimental testing

To perform Step 1 and Step 2: see an example in the `notebooks/bindcraft_to_boltz.ipynb` notebook where we identify the most promising in silico candidates from BindCraft outputs by refolding the sequence against target antigen using Boltz and checking for ipSAE score

To perform Step 3 and Step 4: see an example in the `notebooks/lead_optimization.ipynb` notebook where took the lead binder found from step 2 and used SoluableMPNN to perform mutagenesis on surface and interface residues of the lead binder to further diversify and improve its attributes.

All the inputs required and outputs from the tool are provided in the `examples/` directory.

## Directory Structure
- `examples` - Required inputs and workflow outputs
- `filters` - BindCraft-style json for Binder Filtering
- `notebooks` - Notebooks outlining the pipeline
- `public` - directory where other tools can be installed
- `workflows` - parental directory containing scripts for individual workflows

## Common Questions

### Can this be used to lead optimize variants of antibody formats (scFVs/VHHs)?

As the tool was developed with de novo binders in mind, we have not tested this tool on scFVs and do not recommend using this workflow as the individual tools (e.g. ProteinMPNN/Boltz/ipSAE) have not been validated for this purpose. However as new tools come up may get included in this workflow in the future.

### Can this be used to fold up my binder against a multi-chain protein?
Unfortunately this feature is not implemented at the moment. As a get around, you can add glycines linkers (ex: 20Gs) between the multichain protein and treat them as monomer

### I am bumping into JAX / deep learning library related errors, what should I do?
It is likely that you have compiled the library on a cpu-node without gpu access. This can usually be solved by requesting a gpu interactive session and reinstalling the libraries.

Have more questions / comments? Reach out to Hoyin at chuh@mskcc.org!

## License 
MIT License. See LICENSE file 

**Citation**:
If you are using CARPNN in your research, please cite:
```
@article {Chow2025.12.12.694033,
	author = {Chow, Arthur and Chu, Hoyin and Li, Ruofan and Nalbant, Benan and Dozic, Abdul and Kida, Laura and Lareau, Caleb},
	title = {Sequence and structural determinants of efficacious de novo chimeric antigen receptors},
	year = {2025},
	doi = {10.64898/2025.12.12.694033},
	journal = {bioRxiv}
}
```