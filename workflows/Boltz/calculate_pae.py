#!/usr/bin/env python3

import argparse
import os
import json
import numpy as np
import glob
from Bio.PDB import PDBParser
from itertools import combinations
import pandas as pd

def get_chain_absolute_indices(pdb_path, chain_id):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb_path)
    all_residues = []
    target_indices = []
    for model in structure:
        for chain in model:
            for res in chain.get_residues():
                if res.id[0] != ' ':  # skip hetero/water
                    continue
                all_residues.append((chain.id, res))
    for i, (cid, _) in enumerate(all_residues):
        if cid == chain_id:
            target_indices.append(i)
    return target_indices[0], target_indices[-1]

def extract_boltz_pae(pdb_path, pae_path):
    # Load the PAE matrix from Boltz output
    pae_file = np.load(pae_path)["pae"]
    # 1. Identify all unique chains in the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb_path)
    # Extract IDs for all chains present in the first model
    chains = sorted([chain.id for chain in structure[0]])
    pae_res = {
        "pae_all": np.mean(pae_file)
    }
    # Pre-calculate indices for all chains to avoid repeated PDB parsing
    chain_idx_map = {}
    for chain_id in chains:
        start, end = get_chain_absolute_indices(pdb_path, chain_id)
        chain_idx_map[chain_id] = (start, end)
    # 2. Calculate Intra-chain PAE: pae_chain_<chain>
    for chain_id in chains:
        start, end = chain_idx_map[chain_id]
        # Slice the diagonal block for the specific chain
        intra_slice = pae_file[start:end+1, start:end+1]
        pae_res[f"pae_chain_{chain_id}"] = np.mean(intra_slice)
    # 3. Calculate Inter-chain PAE: pae_interaction_<chain1>_<chain2>
    # combinations(chains, 2) ensures we calculate each pair only once (e.g., A-B)
    for c1, c2 in combinations(chains, 2):
        s1, e1 = chain_idx_map[c1]
        s2, e2 = chain_idx_map[c2]        
        # Interaction 1: Chain 1 relative to Chain 2
        pae_inter1 = np.mean(pae_file[s1:e1+1, s2:e2+1])
        # Interaction 2: Chain 2 relative to Chain 1
        pae_inter2 = np.mean(pae_file[s2:e2+1, s1:e1+1])
        # Calculate the symmetric mean interaction
        pae_interaction = (pae_inter1 + pae_inter2) / 2
        pae_res[f"pae_interaction_{c1}_{c2}"] = pae_interaction

    # For backwards compatibility, make pae_interaction the value between chain A and B
    pae_res[f"pae_interaction"] = pae_res[f"pae_interaction_A_B"]
    pae_res_serialized = {
        k: round(float(v), 4) for k, v in pae_res.items()
    }
    return pae_res_serialized

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-dir", help="Input directory")
    args = parser.parse_args()

    results = {}
    for pdb_path in glob.glob(os.path.join(args.dir, "*.pdb")):
        base = os.path.splitext(os.path.basename(pdb_path))[0]
        pae_path = os.path.join(args.dir, f"pae_{base}.npz")
        pae_out_path = os.path.join(args.dir, f"pae_{base}.json")
        if not os.path.exists(pae_path):
            print(f"Skipping {base}: missing {pae_path}")
            continue
        try:
            results = extract_boltz_pae(pdb_path, pae_path)
            #print(results)
            with open(pae_out_path, "w") as f:
                json.dump(results, f, indent=2)
        except Exception as e:
            print(f"Error with {base}: {e}")
        
if __name__ == "__main__":
    main()