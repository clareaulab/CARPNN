#!/usr/bin/env python3

import argparse
import os
import json
import numpy as np
import glob
from Bio.PDB import PDBParser

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
    pae_file = np.load(pae_path)["pae"]
    chain_a_start, chain_a_end = get_chain_absolute_indices(pdb_path, "A")
    chain_b_start, chain_b_end = get_chain_absolute_indices(pdb_path, "B")
    pae_all = np.mean(pae_file)
    pae_interaction1 = np.mean(pae_file[chain_a_start:chain_a_end+1, chain_b_start:chain_b_end+1])
    pae_interaction2 = np.mean(pae_file[chain_b_start:chain_b_end+1, chain_a_start:chain_a_end+1])
    pae_binder = np.mean(pae_file[chain_a_start:chain_a_end+1, chain_a_start:chain_a_end+1])
    pae_target = np.mean(pae_file[chain_b_start:chain_b_end+1, chain_b_start:chain_b_end+1])
    pae_interaction = (pae_interaction1 + pae_interaction2) / 2
    ## Turning to float so it's serializable
    return {
        "pae_all": float(pae_all),
        "pae_chain_a": float(pae_binder),
        "pae_chain_b": float(pae_target),
        "pae_interaction": float(pae_interaction)
    }

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