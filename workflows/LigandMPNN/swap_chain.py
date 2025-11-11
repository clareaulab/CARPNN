import argparse
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

def copy_and_renumber_chain(chain, new_id, renumber=False):
    new_chain = Chain(new_id)
    new_resid = 1

    for residue in chain:
        if renumber:
            hetfield, _, icode = residue.id
            new_residue_id = (hetfield, new_resid, icode)
            new_residue = Residue(new_residue_id, residue.resname, residue.segid)
            for atom in residue:
                new_residue.add(atom.copy())
            new_chain.add(new_residue)
            new_resid += 1
        else:
            new_chain.add(residue.copy())

    return new_chain

def swap_and_reorder_chains(pdb_input, pdb_output, 
                            chain_a_id='B', chain_b_id='A',
                            renumber=False):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_input)
    model = structure[0]  # Assume one model

    try:
        chain_a = model[chain_a_id]
        chain_b = model[chain_b_id]
    except KeyError:
        print(f"Error: PDB must contain both chain '{chain_a_id}' and chain '{chain_b_id}'.")
        return

    # Swap and optionally renumber
    new_chain_a = copy_and_renumber_chain(chain_a, "A", renumber)
    new_chain_b = copy_and_renumber_chain(chain_b, "B", renumber)

    new_model = Model(0)
    new_model.add(new_chain_a)
    new_model.add(new_chain_b)

    new_structure = Structure("swapped_structure")
    new_structure.add(new_model)

    io = PDBIO()
    io.set_structure(new_structure)
    io.save(pdb_output)

    print(f"Saved swapped and reordered chains to '{pdb_output}'")
    print(f"Original chain '{chain_a_id}' is now 'A' (appears second)")
    print(f"Original chain '{chain_b_id}' is now 'B' (appears first)")
    if renumber:
        print("Residues were renumbered contiguously starting from 1.")

def main():
    parser = argparse.ArgumentParser(description="Swap and reorder two chains in a PDB file.")
    parser.add_argument("input", help="Input PDB file")
    parser.add_argument("output", help="Output PDB file")
    parser.add_argument("--renumber", action="store_true", help="Renumber residues contiguously starting from 1")
    parser.add_argument("--chain-a", default="B", help="Chain ID in original pdb file to turn into chain A (default: B)")
    parser.add_argument("--chain-b", default="A", help="Chain ID in original pdb file to turn into chain B (default: A)")

    args = parser.parse_args()

    swap_and_reorder_chains(
        pdb_input=args.input,
        pdb_output=args.output,
        chain_a_id=args.chain_a,
        chain_b_id=args.chain_b,
        renumber=args.renumber
    )

if __name__ == "__main__":
    main()