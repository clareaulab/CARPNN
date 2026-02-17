### This script is copied and modified from multiple util functions in the BindCraft github ###
### https://github.com/martinpacesa/BindCraft/blob/main/functions/  ###

####################################
################ BioPython functions
####################################
import os
import math
import numpy as np
from collections import defaultdict
from scipy.spatial import cKDTree
from Bio import BiopythonWarning
from Bio.PDB import PDBParser, DSSP, Selection, Polypeptide, PDBIO, Select, Chain, Superimposer
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.Polypeptide import is_aa

####################################
################ PyRosetta functions
####################################
### Import dependencies
#import os
import pyrosetta as pr
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.protocols.simple_moves import AlignChainMover
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.io import pose_from_pose
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

## These imports are replaced by direct function definition
# from .generic_utils import clean_pdb
def clean_pdb(pdb_file):
    # Read the pdb file and filter relevant lines
    with open(pdb_file, 'r') as f_in:
        relevant_lines = [line for line in f_in if line.startswith(('ATOM', 'HETATM', 'MODEL', 'TER', 'END'))]

    # Write the cleaned lines back to the original pdb file
    with open(pdb_file, 'w') as f_out:
        f_out.writelines(relevant_lines)

# analyze sequence composition of design
def validate_design_sequence(sequence, num_clashes, advanced_settings):
    note_array = []

    # Check if protein contains clashes after relaxation
    if num_clashes > 0:
        note_array.append('Relaxed structure contains clashes.')

    # Check if the sequence contains disallowed amino acids
    if advanced_settings["omit_AAs"]:
        restricted_AAs = advanced_settings["omit_AAs"].split(',')
        for restricted_AA in restricted_AAs:
            if restricted_AA in sequence:
                note_array.append('Contains: '+restricted_AA+'!')

    # Analyze the protein
    analysis = ProteinAnalysis(sequence)

    # Calculate the reduced extinction coefficient per 1% solution
    extinction_coefficient_reduced = analysis.molar_extinction_coefficient()[0]
    molecular_weight = round(analysis.molecular_weight() / 1000, 2)
    extinction_coefficient_reduced_1 = round(extinction_coefficient_reduced / molecular_weight * 0.01, 2)

    # Check if the absorption is high enough
    if extinction_coefficient_reduced_1 <= 2:
        note_array.append(f'Absorption value is {extinction_coefficient_reduced_1}, consider adding tryptophane to design.')

    # Join the notes into a single string
    notes = ' '.join(note_array)

    return notes

# temporary function, calculate RMSD of input PDB and trajectory target
def target_pdb_rmsd(trajectory_pdb, starting_pdb, chain_ids_string):
    # Parse the PDB files
    parser = PDBParser(QUIET=True)
    structure_trajectory = parser.get_structure('trajectory', trajectory_pdb)
    structure_starting = parser.get_structure('starting', starting_pdb)
    
    # Extract chain A from trajectory_pdb
    chain_trajectory = structure_trajectory[0]['A']
    
    # Extract the specified chains from starting_pdb
    chain_ids = chain_ids_string.split(',')
    residues_starting = []
    for chain_id in chain_ids:
        chain_id = chain_id.strip()
        chain = structure_starting[0][chain_id]
        for residue in chain:
            if is_aa(residue, standard=True):
                residues_starting.append(residue)
    
    # Extract residues from chain A in trajectory_pdb
    residues_trajectory = [residue for residue in chain_trajectory if is_aa(residue, standard=True)]
    
    # Ensure that both structures have the same number of residues
    min_length = min(len(residues_starting), len(residues_trajectory))
    residues_starting = residues_starting[:min_length]
    residues_trajectory = residues_trajectory[:min_length]
    
    # Collect CA atoms from the two sets of residues
    atoms_starting = [residue['CA'] for residue in residues_starting if 'CA' in residue]
    atoms_trajectory = [residue['CA'] for residue in residues_trajectory if 'CA' in residue]
    
    # Calculate RMSD using structural alignment
    sup = Superimposer()
    sup.set_atoms(atoms_starting, atoms_trajectory)
    rmsd = sup.rms
    
    return round(rmsd, 2)

# detect C alpha clashes for deformed trajectories
def calculate_clash_score(pdb_file, threshold=2.4, only_ca=False):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    atoms = []
    atom_info = []  # Detailed atom info for debugging and processing

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == 'H':  # Skip hydrogen atoms
                        continue
                    if only_ca and atom.get_name() != 'CA':
                        continue
                    atoms.append(atom.coord)
                    atom_info.append((chain.id, residue.id[1], atom.get_name(), atom.coord))

    tree = cKDTree(atoms)
    pairs = tree.query_pairs(threshold)

    valid_pairs = set()
    for (i, j) in pairs:
        chain_i, res_i, name_i, coord_i = atom_info[i]
        chain_j, res_j, name_j, coord_j = atom_info[j]

        # Exclude clashes within the same residue
        if chain_i == chain_j and res_i == res_j:
            continue

        # Exclude directly sequential residues in the same chain for all atoms
        if chain_i == chain_j and abs(res_i - res_j) == 1:
            continue

        # If calculating sidechain clashes, only consider clashes between different chains
        if not only_ca and chain_i == chain_j:
            continue

        valid_pairs.add((i, j))

    return len(valid_pairs)

three_to_one_map = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

# identify interacting residues at the binder interface
def hotspot_residues(trajectory_pdb, binder_chain="B", atom_distance_cutoff=4.0):
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", trajectory_pdb)

    # Get the specified chain
    binder_atoms = Selection.unfold_entities(structure[0][binder_chain], 'A')
    binder_coords = np.array([atom.coord for atom in binder_atoms])

    # Get atoms and coords for the target chain
    target_chain = "A" if binder_chain == "B" else "B"
    target_atoms = Selection.unfold_entities(structure[0][target_chain], 'A')
    target_coords = np.array([atom.coord for atom in target_atoms])

    # Build KD trees for both chains
    binder_tree = cKDTree(binder_coords)
    target_tree = cKDTree(target_coords)

    # Prepare to collect interacting residues
    interacting_residues = {}

    # Query the tree for pairs of atoms within the distance cutoff
    pairs = binder_tree.query_ball_tree(target_tree, atom_distance_cutoff)

    # Process each binder atom's interactions
    for binder_idx, close_indices in enumerate(pairs):
        binder_residue = binder_atoms[binder_idx].get_parent()
        binder_resname = binder_residue.get_resname()

        # Convert three-letter code to single-letter code using the manual dictionary
        if binder_resname in three_to_one_map:
            aa_single_letter = three_to_one_map[binder_resname]
            for close_idx in close_indices:
                target_residue = target_atoms[close_idx].get_parent()
                interacting_residues[binder_residue.id[1]] = aa_single_letter

    return interacting_residues

# calculate secondary structure percentage of design
def calc_ss_percentage(pdb_file, advanced_settings, chain_id="B", atom_distance_cutoff=4.0):
    # Parse the structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    model = structure[0]  # Consider only the first model in the structure

    # Calculate DSSP for the model
    dssp = DSSP(model, pdb_file, dssp=advanced_settings["dssp_path"])

    # Prepare to count residues
    ss_counts = defaultdict(int)
    ss_interface_counts = defaultdict(int)
    plddts_interface = []
    plddts_ss = []

    # Get chain and interacting residues once
    chain = model[chain_id]
    interacting_residues = set(hotspot_residues(pdb_file, chain_id, atom_distance_cutoff).keys())

    for residue in chain:
        residue_id = residue.id[1]
        if (chain_id, residue_id) in dssp:
            ss = dssp[(chain_id, residue_id)][2]  # Get the secondary structure
            ss_type = 'loop'
            if ss in ['H', 'G', 'I']:
                ss_type = 'helix'
            elif ss == 'E':
                ss_type = 'sheet'

            ss_counts[ss_type] += 1

            if ss_type != 'loop':
                # calculate secondary structure normalised pLDDT
                avg_plddt_ss = sum(atom.bfactor for atom in residue) / len(residue)
                plddts_ss.append(avg_plddt_ss)

            if residue_id in interacting_residues:
                ss_interface_counts[ss_type] += 1

                # calculate interface pLDDT
                avg_plddt_residue = sum(atom.bfactor for atom in residue) / len(residue)
                plddts_interface.append(avg_plddt_residue)

    # Calculate percentages
    total_residues = sum(ss_counts.values())
    total_interface_residues = sum(ss_interface_counts.values())

    percentages = calculate_percentages(total_residues, ss_counts['helix'], ss_counts['sheet'])
    interface_percentages = calculate_percentages(total_interface_residues, ss_interface_counts['helix'], ss_interface_counts['sheet'])

    i_plddt = round(sum(plddts_interface) / len(plddts_interface) / 100, 2) if plddts_interface else 0
    ss_plddt = round(sum(plddts_ss) / len(plddts_ss) / 100, 2) if plddts_ss else 0

    return (*percentages, *interface_percentages, i_plddt, ss_plddt)

def calculate_percentages(total, helix, sheet):
    helix_percentage = round((helix / total) * 100,2) if total > 0 else 0
    sheet_percentage = round((sheet / total) * 100,2) if total > 0 else 0
    loop_percentage = round(((total - helix - sheet) / total) * 100,2) if total > 0 else 0

    return helix_percentage, sheet_percentage, loop_percentage


####################################
################ PyRosetta section
####################################

# # Rosetta interface scores
# def score_interface(pdb_file, binder_chain="B"):
#     # load pose
#     pose = pr.pose_from_pdb(pdb_file)

#     # --- ELEGANT MULTI-CHAIN LOGIC START ---
#     # 1. Identify all chain IDs in the pose
#     all_chain_ids = [pose.pdb_info().chain(pose.conformation().chain_begin(i)) 
#                      for i in range(1, pose.num_chains() + 1)]
    
#     # 2. Group all other chains as the target
#     target_chains = "".join([c for c in all_chain_ids if c != binder_chain])
    
#     # 3. Define the interface string (e.g., "B_AC")
#     # This tells Rosetta: Partner 1 is binder_chain, Partner 2 is the rest
#     interface_string = f"{binder_chain}_{target_chains}"
#     # --- ELEGANT MULTI-CHAIN LOGIC END ---

#     # analyze interface statistics
#     iam = InterfaceAnalyzerMover()
    
#     # Apply the logical binary interface
#     iam.set_interface(interface_string)
    
#     scorefxn = pr.get_fa_scorefxn()
#     iam.set_scorefunction(scorefxn)
#     iam.set_compute_packstat(True)
#     iam.set_compute_interface_energy(True)
#     iam.set_calc_dSASA(True)
#     iam.set_calc_hbond_sasaE(True)
#     iam.set_compute_interface_sc(True)
#     iam.set_pack_separated(True)
#     iam.apply(pose)
#     # # load pose
#     # pose = pr.pose_from_pdb(pdb_file)

#     # # analyze interface statistics
#     # iam = InterfaceAnalyzerMover()
#     # scorefxn = pr.get_fa_scorefxn()
#     # iam.set_scorefunction(scorefxn)
#     # iam.set_compute_packstat(True)
#     # iam.set_compute_interface_energy(True)
#     # iam.set_calc_dSASA(True)
#     # iam.set_calc_hbond_sasaE(True)
#     # iam.set_compute_interface_sc(True)
#     # iam.set_pack_separated(True)
#     # iam.apply(pose)

#     # Initialize dictionary with all amino acids
#     interface_AA = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}

#     # Initialize list to store PDB residue IDs at the interface
#     interface_residues_set = hotspot_residues(pdb_file, binder_chain)
#     interface_residues_pdb_ids = []
    
#     # Iterate over the interface residues
#     for pdb_res_num, aa_type in interface_residues_set.items():
#         # Increase the count for this amino acid type
#         interface_AA[aa_type] += 1

#         # Append the binder_chain and the PDB residue number to the list
#         interface_residues_pdb_ids.append(f"{binder_chain}{pdb_res_num}")

#     # count interface residues
#     interface_nres = len(interface_residues_pdb_ids)

#     # Convert the list into a comma-separated string
#     interface_residues_pdb_ids_str = ','.join(interface_residues_pdb_ids)

#     # Calculate the percentage of hydrophobic residues at the interface of the binder
#     hydrophobic_aa = set('ACFILMPVWY')
#     hydrophobic_count = sum(interface_AA[aa] for aa in hydrophobic_aa)
#     if interface_nres != 0:
#         interface_hydrophobicity = (hydrophobic_count / interface_nres) * 100
#     else:
#         interface_hydrophobicity = 0

#     # retrieve statistics
#     interfacescore = iam.get_all_data()
#     interface_sc = interfacescore.sc_value # shape complementarity
#     interface_interface_hbonds = interfacescore.interface_hbonds # number of interface H-bonds
#     interface_dG = iam.get_interface_dG() # interface dG
#     interface_dSASA = iam.get_interface_delta_sasa() # interface dSASA (interface surface area)
#     interface_packstat = iam.get_interface_packstat() # interface pack stat score
#     interface_dG_SASA_ratio = interfacescore.dG_dSASA_ratio * 100 # ratio of dG/dSASA (normalised energy for interface area size)
#     buns_filter = XmlObjects.static_get_filter('<BuriedUnsatHbonds report_all_heavy_atom_unsats="true" scorefxn="scorefxn" ignore_surface_res="false" use_ddG_style="true" dalphaball_sasa="1" probe_radius="1.1" burial_cutoff_apo="0.2" confidence="0" />')
#     interface_delta_unsat_hbonds = buns_filter.report_sm(pose)

#     if interface_nres != 0:
#         interface_hbond_percentage = (interface_interface_hbonds / interface_nres) * 100 # Hbonds per interface size percentage
#         interface_bunsch_percentage = (interface_delta_unsat_hbonds / interface_nres) * 100 # Unsaturated H-bonds per percentage
#     else:
#         interface_hbond_percentage = None
#         interface_bunsch_percentage = None

#     # calculate binder energy score
#     chain_design = ChainSelector(binder_chain)
#     tem = pr.rosetta.core.simple_metrics.metrics.TotalEnergyMetric()
#     tem.set_scorefunction(scorefxn)
#     tem.set_residue_selector(chain_design)
#     binder_score = tem.calculate(pose)

#     # calculate binder SASA fraction
#     bsasa = pr.rosetta.core.simple_metrics.metrics.SasaMetric()
#     bsasa.set_residue_selector(chain_design)
#     binder_sasa = bsasa.calculate(pose)

#     if binder_sasa > 0:
#         interface_binder_fraction = (interface_dSASA / binder_sasa) * 100
#     else:
#         interface_binder_fraction = 0

#     # calculate surface hydrophobicity
#     binder_pose = {pose.pdb_info().chain(pose.conformation().chain_begin(i)): p for i, p in zip(range(1, pose.num_chains()+1), pose.split_by_chain())}[binder_chain]

#     layer_sel = pr.rosetta.core.select.residue_selector.LayerSelector()
#     layer_sel.set_layers(pick_core = False, pick_boundary = False, pick_surface = True)
#     surface_res = layer_sel.apply(binder_pose)

#     exp_apol_count = 0
#     total_count = 0 
    
#     # count apolar and aromatic residues at the surface
#     for i in range(1, len(surface_res) + 1):
#         if surface_res[i] == True:
#             res = binder_pose.residue(i)

#             # count apolar and aromatic residues as hydrophobic
#             if res.is_apolar() == True or res.name() == 'PHE' or res.name() == 'TRP' or res.name() == 'TYR':
#                 exp_apol_count += 1
#             total_count += 1

#     surface_hydrophobicity = exp_apol_count/total_count

#     # output interface score array and amino acid counts at the interface
#     interface_scores = {
#     'binder_score': binder_score,
#     'surface_hydrophobicity': surface_hydrophobicity,
#     'interface_sc': interface_sc,
#     'interface_packstat': interface_packstat,
#     'interface_dG': interface_dG,
#     'interface_dSASA': interface_dSASA,
#     'interface_dG_SASA_ratio': interface_dG_SASA_ratio,
#     'interface_fraction': interface_binder_fraction,
#     'interface_hydrophobicity': interface_hydrophobicity,
#     'interface_nres': interface_nres,
#     'interface_interface_hbonds': interface_interface_hbonds,
#     'interface_hbond_percentage': interface_hbond_percentage,
#     'interface_delta_unsat_hbonds': interface_delta_unsat_hbonds,
#     'interface_delta_unsat_hbonds_percentage': interface_bunsch_percentage
#     }

#     # round to two decimal places
#     interface_scores = {k: round(v, 2) if isinstance(v, float) else v for k, v in interface_scores.items()}

#     return interface_scores, interface_AA, interface_residues_pdb_ids_str

def score_interface(pdb_file, binder_chain="B"):
    # load pose
    pose = pr.pose_from_pdb(pdb_file)

    # --- ELEGANT MULTI-CHAIN LOGIC START ---
    all_chain_ids = [pose.pdb_info().chain(pose.conformation().chain_begin(i)) 
                     for i in range(1, pose.num_chains() + 1)]
    target_chains = "".join([c for c in all_chain_ids if c != binder_chain])
    interface_string = f"{binder_chain}_{target_chains}"
    # --- ELEGANT MULTI-CHAIN LOGIC END ---

    iam = InterfaceAnalyzerMover()
    iam.set_interface(interface_string)
    
    scorefxn = pr.get_fa_scorefxn()
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True)
    iam.set_calc_hbond_sasaE(True)
    iam.set_compute_interface_sc(True)
    iam.set_pack_separated(True)
    iam.apply(pose)

    # Initialize interface analysis
    interface_AA = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}
    interface_residues_set = hotspot_residues(pdb_file, binder_chain)
    interface_residues_pdb_ids = []
    
    for pdb_res_num, aa_type in interface_residues_set.items():
        interface_AA[aa_type] += 1
        interface_residues_pdb_ids.append(f"{binder_chain}{pdb_res_num}")

    interface_nres = len(interface_residues_pdb_ids)
    interface_residues_pdb_ids_str = ','.join(interface_residues_pdb_ids)

    # Hydrophobicity
    hydrophobic_aa = set('ACFILMPVWY')
    hydrophobic_count = sum(interface_AA[aa] for aa in hydrophobic_aa)
    interface_hydrophobicity = (hydrophobic_count / interface_nres) * 100 if interface_nres != 0 else 0

    # Retrieve statistics from IAM
    interfacescore = iam.get_all_data()
    interface_sc = interfacescore.sc_value
    interface_interface_hbonds = interfacescore.interface_hbonds
    interface_dG = iam.get_interface_dG()
    interface_dSASA = iam.get_interface_delta_sasa()
    interface_packstat = iam.get_interface_packstat()
    interface_dG_SASA_ratio = interfacescore.dG_dSASA_ratio * 100

    # --- CONDITIONAL BURIED UNSAT CALCULATION ---
    # Rosetta's BuriedUnsatHbondFilter with use_ddG_style crashes if num_chains > 3
    if pose.num_chains() < 4:
        buns_filter = XmlObjects.static_get_filter('<BuriedUnsatHbonds report_all_heavy_atom_unsats="true" scorefxn="scorefxn" ignore_surface_res="false" use_ddG_style="true" dalphaball_sasa="1" probe_radius="1.1" burial_cutoff_apo="0.2" confidence="0" />')
        interface_delta_unsat_hbonds = buns_filter.report_sm(pose)
    else:
        # Gracefully skip if too many chains
        interface_delta_unsat_hbonds = np.NaN
        print(f"Skipping BuriedUnsatHbond calculation for {pdb_file}: too many chains ({pose.num_chains()})")

    if interface_nres != 0:
        interface_hbond_percentage = (interface_interface_hbonds / interface_nres) * 100
        interface_bunsch_percentage = (interface_delta_unsat_hbonds / interface_nres) * 100
    else:
        interface_hbond_percentage = 0.0
        interface_bunsch_percentage = 0.0

    # --- REMAINING METRICS (Binder Score/SASA/Hydrophobicity) ---
    chain_design = ChainSelector(binder_chain)
    tem = pr.rosetta.core.simple_metrics.metrics.TotalEnergyMetric()
    tem.set_scorefunction(scorefxn)
    tem.set_residue_selector(chain_design)
    binder_score = tem.calculate(pose)

    bsasa = pr.rosetta.core.simple_metrics.metrics.SasaMetric()
    bsasa.set_residue_selector(chain_design)
    binder_sasa = bsasa.calculate(pose)
    interface_binder_fraction = (interface_dSASA / binder_sasa) * 100 if binder_sasa > 0 else 0

    # Surface Hydrophobicity
    binder_pose = {pose.pdb_info().chain(pose.conformation().chain_begin(i)): p for i, p in zip(range(1, pose.num_chains()+1), pose.split_by_chain())}[binder_chain]
    layer_sel = pr.rosetta.core.select.residue_selector.LayerSelector()
    layer_sel.set_layers(pick_core = False, pick_boundary = False, pick_surface = True)
    surface_res = layer_sel.apply(binder_pose)

    exp_apol_count = 0
    total_count = 0 
    for i in range(1, len(surface_res) + 1):
        if surface_res[i]:
            res = binder_pose.residue(i)
            if res.is_apolar() or res.name() in ['PHE', 'TRP', 'TYR']:
                exp_apol_count += 1
            total_count += 1
    surface_hydrophobicity = exp_apol_count/total_count if total_count > 0 else 0

    interface_scores = {
        'binder_score': binder_score,
        'surface_hydrophobicity': surface_hydrophobicity,
        'interface_sc': interface_sc,
        'interface_packstat': interface_packstat,
        'interface_dG': interface_dG,
        'interface_dSASA': interface_dSASA,
        'interface_dG_SASA_ratio': interface_dG_SASA_ratio,
        'interface_fraction': interface_binder_fraction,
        'interface_hydrophobicity': interface_hydrophobicity,
        'interface_nres': interface_nres,
        'interface_interface_hbonds': interface_interface_hbonds,
        'interface_hbond_percentage': interface_hbond_percentage,
        'interface_delta_unsat_hbonds': interface_delta_unsat_hbonds,
        'interface_delta_unsat_hbonds_percentage': interface_bunsch_percentage
    }

    interface_scores = {k: round(v, 2) if isinstance(v, float) else v for k, v in interface_scores.items()}

    return interface_scores, interface_AA, interface_residues_pdb_ids_str

# align pdbs to have same orientation
def align_pdbs(reference_pdb, align_pdb, reference_chain_id, align_chain_id):
    # initiate poses
    reference_pose = pr.pose_from_pdb(reference_pdb)
    align_pose = pr.pose_from_pdb(align_pdb)

    align = AlignChainMover()
    align.pose(reference_pose)

    # If the chain IDs contain commas, split them and only take the first value
    reference_chain_id = reference_chain_id.split(',')[0]
    align_chain_id = align_chain_id.split(',')[0]

    # Get the chain number corresponding to the chain ID in the poses
    reference_chain = pr.rosetta.core.pose.get_chain_id_from_chain(reference_chain_id, reference_pose)
    align_chain = pr.rosetta.core.pose.get_chain_id_from_chain(align_chain_id, align_pose)

    align.source_chain(align_chain)
    align.target_chain(reference_chain)
    align.apply(align_pose)

    # Overwrite aligned pdb
    align_pose.dump_pdb(align_pdb)
    clean_pdb(align_pdb)

# calculate the rmsd without alignment
def unaligned_rmsd(reference_pdb, align_pdb, reference_chain_id, align_chain_id):
    reference_pose = pr.pose_from_pdb(reference_pdb)
    align_pose = pr.pose_from_pdb(align_pdb)

    # Define chain selectors for the reference and align chains
    reference_chain_selector = ChainSelector(reference_chain_id)
    align_chain_selector = ChainSelector(align_chain_id)

    # Apply selectors to get residue subsets
    reference_chain_subset = reference_chain_selector.apply(reference_pose)
    align_chain_subset = align_chain_selector.apply(align_pose)

    # Convert subsets to residue index vectors
    reference_residue_indices = get_residues_from_subset(reference_chain_subset)
    align_residue_indices = get_residues_from_subset(align_chain_subset)

    # Create empty subposes
    reference_chain_pose = pr.Pose()
    align_chain_pose = pr.Pose()

    # Fill subposes
    pose_from_pose(reference_chain_pose, reference_pose, reference_residue_indices)
    pose_from_pose(align_chain_pose, align_pose, align_residue_indices)

    # Calculate RMSD using the RMSDMetric
    rmsd_metric = RMSDMetric()
    rmsd_metric.set_comparison_pose(reference_chain_pose)
    rmsd = rmsd_metric.calculate(align_chain_pose)

    return round(rmsd, 2)

# Relax designed structure
def pr_relax(pdb_file, relaxed_pdb_path):
    if not os.path.exists(relaxed_pdb_path):
        # Generate pose
        pose = pr.pose_from_pdb(pdb_file)
        start_pose = pose.clone()

        ### Generate movemaps
        mmf = MoveMap()
        mmf.set_chi(True) # enable sidechain movement
        mmf.set_bb(True) # enable backbone movement, can be disabled to increase speed by 30% but makes metrics look worse on average
        mmf.set_jump(False) # disable whole chain movement

        # Run FastRelax
        fastrelax = FastRelax()
        scorefxn = pr.get_fa_scorefxn()
        fastrelax.set_scorefxn(scorefxn)
        fastrelax.set_movemap(mmf) # set MoveMap
        fastrelax.max_iter(200) # default iterations is 2500
        fastrelax.min_type("lbfgs_armijo_nonmonotone")
        fastrelax.constrain_relax_to_start_coords(True)
        fastrelax.apply(pose)

        # Align relaxed structure to original trajectory
        align = AlignChainMover()
        align.source_chain(0)
        align.target_chain(0)
        align.pose(start_pose)
        align.apply(pose)

        # Copy B factors from start_pose to pose
        for resid in range(1, pose.total_residue() + 1):
            if pose.residue(resid).is_protein():
                # Get the B factor of the first heavy atom in the residue
                bfactor = start_pose.pdb_info().bfactor(resid, 1)
                for atom_id in range(1, pose.residue(resid).natoms() + 1):
                    pose.pdb_info().bfactor(resid, atom_id, bfactor)

        # output relaxed and aligned PDB
        pose.dump_pdb(relaxed_pdb_path)
        clean_pdb(relaxed_pdb_path)

### Starting from this section this is custom
import glob
import argparse
from tqdm import tqdm
import pandas as pd

def make_relaxed_pdbdir(pdbdir,outdir=None):
    """
    For a directory of pdbs, make a new directory of relaxed pdbs with the same name at the new directory
    """
    pdb_paths = glob.glob(f"{pdbdir}/*.pdb")
    outdir = f"{pdbdir}/relaxed" if outdir is None else outdir
    os.makedirs(outdir,exist_ok=True)
    for pdb in tqdm(pdb_paths):
        pdb_basename = os.path.basename(pdb)
        pdb_relaxed = f"{outdir}/{pdb_basename}"
        pr_relax(pdb, pdb_relaxed)
    print(f"Finished writing relaxed pdbs at {outdir}")
    return outdir
    

def interface_score_dataframe(pdbdir,binder_chain="A",relaxed_pdbdir=None,no_relax=False,regex="*.pdb",is_colabfold=False):
    """
    Given a directory of pdbs, create a relaxed pdb and score according pyrosetta metrics as calculated by BindCraft
    """
    pdb_paths = glob.glob(f"{pdbdir}/{regex}")
    # Creates a partial bindcraft parameters file that only contains the path to the DSSP software
    advanced_settings = {"dssp_path":DSSP_PATH}
    # Make relaxed version of the pdbs if not given
    if (not no_relax) and (not is_colabfold):
        if relaxed_pdbdir is None:
            print("Relaxed pdb not detected. Creating relaxed pdbs")
            relaxed_pdbdir = make_relaxed_pdbdir(pdbdir)
    else:
        print("Calculating BindCraft interface scores without relaxing.")
        relaxed_pdbdir = pdbdir
    rows = []

    for pdb in tqdm(pdb_paths):
        ## Obtain the relaxed version of the pdb
        pdb_basename = os.path.basename(pdb)
        relaxed_pdb = f"{relaxed_pdbdir}/{pdb_basename}"
        if is_colabfold:
            if "_relaxed_" in pdb_basename:
                # Skip to avoid calculating twice
                continue
            else:
                relaxed_pdb = f"{relaxed_pdbdir}/{pdb_basename}".replace("_unrelaxed_","_relaxed_")
        # Calculate clashes before and after relaxation
        num_clashes_mpnn = calculate_clash_score(pdb)
        num_clashes_mpnn_relaxed = calculate_clash_score(relaxed_pdb)

        ## Score the interface (of the relaxed pdb)
        mpnn_interface_scores, mpnn_interface_AA, mpnn_interface_residues = score_interface(relaxed_pdb,binder_chain=binder_chain)

        # secondary structure content of the binder
        mpnn_alpha, mpnn_beta, mpnn_loops, mpnn_alpha_interface, mpnn_beta_interface, mpnn_loops_interface, mpnn_i_plddt, mpnn_ss_plddt = calc_ss_percentage(pdb, advanced_settings, chain_id=binder_chain)

        ## Statistic to track
        statistic_dict= {
            'i_pLDDT': mpnn_i_plddt,
            'ss_pLDDT': mpnn_ss_plddt,
            'Unrelaxed_Clashes': num_clashes_mpnn,
            'Relaxed_Clashes': num_clashes_mpnn_relaxed,
            'Binder_Energy_Score': mpnn_interface_scores['binder_score'],
            'Surface_Hydrophobicity': mpnn_interface_scores['surface_hydrophobicity'],
            'ShapeComplementarity': mpnn_interface_scores['interface_sc'],
            'PackStat': mpnn_interface_scores['interface_packstat'],
            'dG': mpnn_interface_scores['interface_dG'],
            'dSASA': mpnn_interface_scores['interface_dSASA'], 
            'dG/dSASA': mpnn_interface_scores['interface_dG_SASA_ratio'],
            'Interface_SASA_%': mpnn_interface_scores['interface_fraction'],
            'Interface_Hydrophobicity': mpnn_interface_scores['interface_hydrophobicity'],
            'n_InterfaceResidues': mpnn_interface_scores['interface_nres'],
            'n_InterfaceHbonds': mpnn_interface_scores['interface_interface_hbonds'],
            'InterfaceHbondsPercentage': mpnn_interface_scores['interface_hbond_percentage'],
            'n_InterfaceUnsatHbonds': mpnn_interface_scores['interface_delta_unsat_hbonds'],
            'InterfaceUnsatHbondsPercentage': mpnn_interface_scores['interface_delta_unsat_hbonds_percentage'],
            'InterfaceAAs': mpnn_interface_AA,
            'Interface_Helix%': mpnn_alpha_interface,
            'Interface_BetaSheet%': mpnn_beta_interface,
            'Interface_Loop%': mpnn_loops_interface,
            'Binder_Helix%': mpnn_alpha,
            'Binder_BetaSheet%': mpnn_beta,
            'Binder_Loop%': mpnn_loops,
        }

        ## Non-statistic values
        statistic_dict["binder_path"] = pdb
        statistic_dict["binder_relaxed_path"] = relaxed_pdb
        rows.append(statistic_dict)
        
    interface_df = pd.DataFrame(rows)
    return interface_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    DALPHABALL_PATH="/data1/lareauc/users/chuh/softwares/BindCraft/functions/DAlphaBall.gcc"
    DSSP_PATH="/data1/lareauc/users/chuh/softwares/BindCraft/functions/dssp"

    # I/O Arguments
    parser.add_argument("-pdbdir", type=str,required=True,help='The directory containing the pdbs to be scored')
    parser.add_argument("-relaxed_pdbdir", type=str,default=None,help='The directory containing the relaxed pdbs. Must have the same number and same names. Creates a new one at the specified directory if not specified')
    parser.add_argument("-regex", type=str,default="*.pdb",help='Glob string to be used to match pdbs in the given directory. Default is *.pdb. For processing batched best ranked ColabFold output use *_relaxed_rank_001_*.pdb')
    parser.add_argument("-binder_chain", type=str, default="A",help='The chain of the binder in each pdb file. Default A which is different than the BindCraft one')
    parser.add_argument("-no_relax", action="store_true",default=False,help='If specified, do not relax the input pdbs and directly calculate interface-related metrics by assuming they are already relaxed')
    parser.add_argument("-is_colabfold", action="store_true",default=False,help='If specified, assume ColabFold naming convention between unrelaxed vs. relaxed pdb. (aka we will search for _unrelaxed_ and _relaxed_ pair)')
    parser.add_argument("-out", type=str, default="./binding_interface.tsv", help="Path to the output score file")
    parser.add_argument("-seed", type=int, default=42, help="Seed to use for Pyrosetta")
    
    # Add a special argument that specifically alows for mapping between relaxed and unrelaxed pdbs based on ColabFold outputs
    args = parser.parse_args()

    ## Initialize pyrosetta
    #"-run:constant_seed 1 -jran 10"
    pr_init_string = f"-ignore_unrecognized_res -ignore_zero_occupancy -mute all -holes:dalphaball {DALPHABALL_PATH} -corrections::beta_nov16 true -relax:default_repeats 1 -run:constant_seed 1 -jran {args.seed}"
    pr.init(pr_init_string)

    interface_df = interface_score_dataframe(args.pdbdir,binder_chain=args.binder_chain,relaxed_pdbdir=args.relaxed_pdbdir,no_relax=args.no_relax,regex=args.regex,is_colabfold=args.is_colabfold)
    interface_df.to_csv(args.out,sep="\t",index=False)
    print("Done writing interface df")
