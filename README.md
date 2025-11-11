# CAR-PNN workflow: Chimeric-Antigen-Receptor ProteinMPNN Workflow

![project_logo](assets/logo.png)

CAR-PNN is a HPC-ready workflow that performs leads optimization of de novo binder into CAR candidates using variety of protein design tools. 

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
### Step 2: Install environments required for individual workflows

If you already have these tools installed else where in your HPC environment you can just replace the paths to their corresponding installation directory without reinstalling. Some of the github repository required include:

- LigandMPNN (For running Protein/SoluableMPNN, https://github.com/dauparas/LigandMPNN)
- Boltz2 (For structure prediction and filtering, https://github.com/jwohlwend/boltz)
- ColabFold container (For using ColabSearch module to perform MSA using local database, https://github.com/YoshitakaMo/localcolabfold)
- BindCraft (For calculation of interface metrics, https://github.com/martinpacesa/BindCraft)
- ipSAE (Already provided in this repo, for interaction scoring, https://github.com/DunbrackLab/IPSAE/tree/main)

For details, read about how to set up the individual workflows in the README files in the individual workflow directorys under `workflows/`. The workflows are decoupled into individual scripts such that if you already have another workflow you are using to create a certain component (ex: binder generation / MSA) they should be easily swappable to your needs.

## How to Run

In short, the workflow consist of 4 major steps:

- Step 1: Refold candidate binders via Boltz against target antigen(ex: binder sequences from tools like RFDiffusion, BindCraft, BoltzGen, ProteinHunter etc.)
- Step 2: Identify lead de novo binder (experimentally validated or have high ipSAE score > 0.8)
- Step 3: Mutagensis (non-interface/interface) residues of the lead binder and perform sequence level filter (ex: charge / immunogenicity)
- Step 4: Refold mutagenized sequences via Boltz and identify top binders for experimental testing

Workflow:
### Step 1 (Not included in this repo): Generate protein binders against your target of interest
- Some recommended tools:
    - RFDiffusion, BindCraft, BoltzGen, ProteinHunter

### Step 2: Identify lead binder candidate from protein design tools. We recommend this should be a binder with ipSAE of at least 0.8 (higher the better)

In the `notebooks/bindcraft_to_boltz.ipynb` notebook, we provide an example workflow that identifies most promising in silico candidates from BindCraft outputs by refolding the sequence against target antigen using Boltz and checking for ipSAE score

### Step 3: Generate candidate sequences to evolve to using SoluableMPNN. 
In the first part of the `notebooks/lead_optimization.ipynb` notebook, we take the best candidate identified from step 2 and used SoluableMPNN to generate mutagenesis and performed basic sequence filtering

### Step 4: Refold the candidate sequences against the target antigen
In the second part of the `notebooks/lead_optimization.ipynb` notebook, we refold the mutagenized sequences and calculated additional interface metrics such as those used by BindCraft and ipSAE

### Step 5: Calculate additional structural metrics based on output and aggregate the outputs across splits
In the last part of the `notebooks/lead_optimization.ipynb` notebook, we apply a BindCraft style filter to eliminate sequences and ranked the final outputs by ipSAE (high to low)

Example input and outputs of the workflow can be found in the `examples` directory.

Please cite this work and the work where each workflows are based on (listed in individual workflow README) if you found this work to be helpful.


License: MIT 