#!/usr/bin/env python
import sys
import os
import json
import argparse
from Bio.PDB import PDBParser, NeighborSearch, Selection
from Bio.PDB.Polypeptide import is_aa  # ✅ Correct import

def get_interface_residues(pdb_file, binder_chain_id, target_chain_id, distance_cutoff, invert):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    model = structure[0]

    try:
        binder_chain = model[binder_chain_id]
        target_chain = model[target_chain_id]
    except KeyError as e:
        raise ValueError(f"Chain ID not found: {e}")

    target_atoms = Selection.unfold_entities(target_chain, 'A')
    ns = NeighborSearch(target_atoms)

    interface_residues = set()
    all_binder_residues = set()

    for residue in binder_chain:
        if not is_aa(residue, standard=True):
            continue
        res_id = f"{binder_chain_id}{residue.id[1]}"
        all_binder_residues.add(res_id)
        for atom in residue:
            neighbors = ns.search(atom.coord, distance_cutoff, level='R')
            if any(n.get_parent().id == target_chain_id for n in neighbors):
                interface_residues.add(res_id)
                break

    if invert:
        selected_residues = all_binder_residues - interface_residues
    else:
        selected_residues = interface_residues

    return sorted(selected_residues, key=lambda x: int(x[1:]))

def main():
    parser = argparse.ArgumentParser(description="Extract interface (or non-interface) residues from a PDB file.")
    parser.add_argument("pdb_file", help="Path to the PDB file")
    parser.add_argument("binder_chain", help="Binder chain ID (e.g., B)")
    parser.add_argument("target_chain", help="Target chain ID (e.g., A)")
    parser.add_argument("--cutoff", type=float, default=4.0, help="Distance cutoff in Å (default: 4.0)")
    parser.add_argument("--invert", action="store_true", help="Return non-interface residues instead")
    parser.add_argument("--output_json", default=None, help="Optional path to save the output JSON")

    args = parser.parse_args()

    try:
        residues = get_interface_residues(
            args.pdb_file,
            args.binder_chain,
            args.target_chain,
            args.cutoff,
            args.invert
        )
        output = {args.pdb_file: " ".join(residues)}

        # Determine output path if not provided
        if args.output_json is None:
            pdb_base = os.path.splitext(os.path.basename(args.pdb_file))[0]
            output_path = os.path.join(os.path.dirname(args.pdb_file), f"{pdb_base}_interface.jsonl")
        else:
            output_path = args.output_json

        with open(output_path, "w") as f:
            f.write(json.dumps(output) + "\n")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        exit(1)

if __name__ == "__main__":
    main()