import pandas as pd
import numpy as np
import os
import tqdm

from Bio.PDB import PDBParser, MMCIFParser, Superimposer
from Bio import pairwise2

import argparse

def three_to_one(resname):
    """Convert 3-letter amino acid code to 1-letter."""
    AA_3TO1 = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLU": "E",
        "GLN": "Q",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V"
    }
    return AA_3TO1[resname]


def get_sequence_and_residues(chain):
    seq = []
    residues = []
    for res in chain:
        if 'CA' in res:
            try:
                aa = three_to_one(res.get_resname())
                seq.append(aa)
                residues.append(res)
            except KeyError:
                continue  # skip non-standard residues
    return seq, residues

def aligned_chain_rmsd_over_other_chain(
    reference_pdb,
    align_pdb,
    reference_align_chain_id,
    align_align_chain_id,
    reference_rmsd_chain_id,
    align_rmsd_chain_id,
    seq_align_first=True,
    seq_align_second=False
):
    """
    A pymol-style alignment where sequence alignment is first computed in specific chains.
    Then the CA of the sequence-aligned residues between the two PDBs are structurally aligned.
    RMSD is then computed over the specified chains.
    
    Args:
        reference_pdb: path to reference PDB file.
        align_pdb: path to structure to be aligned.
        reference_align_chain_id: chain ID in reference structure used for alignment.
        align_align_chain_id: chain ID in align structure used for alignment.
        reference_rmsd_chain_id: chain ID in reference used for RMSD calculation.
        align_rmsd_chain_id: chain ID in aligned structure used for RMSD calculation.
        seq_align_first: whether to sequence-align before structural alignment.
        seq_align_second: whether to compute RMSD only over sequence-aligned residues in the RMSD chains.

    Returns:
        RMSD (float)
    """
    if reference_pdb.endswith(".pdb"):
        ref_structure = PDBParser(QUIET=True).get_structure("ref", reference_pdb)
    elif reference_pdb.endswith(".cif"):
        ref_structure = MMCIFParser(QUIET=True).get_structure("ref", reference_pdb)

    if align_pdb.endswith(".pdb"):
        align_structure = PDBParser(QUIET=True).get_structure("align", align_pdb)
    elif align_pdb.endswith(".cif"):
        align_structure = MMCIFParser(QUIET=True).get_structure("align", align_pdb)

    # Alignment chains (used for alignment and superposition)
    ref_align_chain = ref_structure[0][reference_align_chain_id]
    align_align_chain = align_structure[0][align_align_chain_id]

    # RMSD chains (used for final RMSD calculation)
    ref_rmsd_chain = ref_structure[0][reference_rmsd_chain_id]
    align_rmsd_chain = align_structure[0][align_rmsd_chain_id]

    # Sequence and residue extraction from alignment chains
    ref_seq, ref_residues = get_sequence_and_residues(ref_align_chain)
    align_seq, align_residues = get_sequence_and_residues(align_align_chain)

    ref_seq_str = ''.join(ref_seq)
    align_seq_str = ''.join(align_seq)

    atoms_ref = []
    atoms_align = []

    # First sequence alignment (used for superposition)
    if seq_align_first:
        alignments = pairwise2.align.globalxx(ref_seq_str, align_seq_str)
        aligned1, aligned2 = alignments[0].seqA, alignments[0].seqB

        idx_ref = idx_align = 0
        aligned_count = 0

        for a1, a2 in zip(aligned1, aligned2):
            if a1 != '-' and a2 != '-':
                atoms_ref.append(ref_residues[idx_ref]['CA'])
                atoms_align.append(align_residues[idx_align]['CA'])
                aligned_count += 1
            if a1 != '-':
                idx_ref += 1
            if a2 != '-':
                idx_align += 1

        #print(f"[Alignment 1: chain {align_align_chain_id} to chain {reference_align_chain_id}] Aligned residues: {aligned_count}")
        if aligned_count < len(ref_seq)*0.8:
            print("Warning: The total sequence-aligned residues are less than 80% of the reference sequence. This result is likely inacccurate")
    else:
        # Directly match residues by index without sequence alignment
        for ref_res, align_res in zip(ref_residues, align_residues):
            if 'CA' in ref_res and 'CA' in align_res:
                atoms_ref.append(ref_res['CA'])
                atoms_align.append(align_res['CA'])
        #print(f"[Alignment 1 chain {align_align_chain_id} to chain {reference_align_chain_id}] Aligned residues (by index): {len(atoms_ref)}")

    # Superimpose align_structure onto ref_structure
    #print(f"atoms ref: {len(atoms_ref)}")
    #print(f"atoms align: {len(atoms_align)}")

    sup = Superimposer()
    sup.set_atoms(atoms_ref, atoms_align)
    sup.apply(align_structure.get_atoms())

    #print("Superimposition done.")

    # Extract coordinates for RMSD computation
    if seq_align_second:
        # Sequence-align RMSD chains and compute RMSD only on aligned residues
        ref_seq2, ref_res2 = get_sequence_and_residues(ref_rmsd_chain)
        align_seq2, align_res2 = get_sequence_and_residues(align_rmsd_chain)

        ref_seq2_str = ''.join(ref_seq2)
        align_seq2_str = ''.join(align_seq2)

        alignments2 = pairwise2.align.globalxx(ref_seq2_str, align_seq2_str)
        aligned1, aligned2 = alignments2[0].seqA, alignments2[0].seqB

        idx_ref = idx_align = 0
        ref_coords = []
        align_coords = []
        aligned_count2 = 0

        for a1, a2 in zip(aligned1, aligned2):
            if a1 != '-' and a2 != '-':
                if 'CA' in ref_res2[idx_ref] and 'CA' in align_res2[idx_align]:
                    ref_coords.append(ref_res2[idx_ref]['CA'].get_coord())
                    align_coords.append(align_res2[idx_align]['CA'].get_coord())
                    aligned_count2 += 1
            if a1 != '-':
                idx_ref += 1
            if a2 != '-':
                idx_align += 1

        #print(f"[Alignment 2: chain {align_rmsd_chain_id} to chain {reference_rmsd_chain_id}] Aligned residues: {aligned_count2}")
        if aligned_count2 < len(ref_seq2)*0.8:
            print("Warning: The total sequence-aligned residues are less than 80% of the reference sequence. This result is likely inacccurate")
    else:
        ref_coords = [res['CA'].get_coord() for res in ref_rmsd_chain if 'CA' in res]
        align_coords = [res['CA'].get_coord() for res in align_rmsd_chain if 'CA' in res]
        #print(f"[Alignment 2: chain {align_rmsd_chain_id} to chain {reference_rmsd_chain_id}] Aligned residues (by index): {len(ref_coords)}")

    #print(f"ref coords (chain {ref_rmsd_chain}): {len(ref_coords)}")
    #print(f"align coords (chain {align_rmsd_chain}): {len(align_coords)}")

    if len(ref_coords) != len(align_coords):
        raise ValueError("Number of C-alpha atoms in RMSD chains must match after filtering.")

    ref_coords = np.array(ref_coords)
    align_coords = np.array(align_coords)


    diff = ref_coords - align_coords
    rmsd = np.sqrt((diff ** 2).sum() / len(ref_coords))
    return round(rmsd, 3)

def compute_rmsd_between_pdbs(binder_complex_pdb, binder_backbone_pdb, binder_monomer_pdb,
                               compute_monomer_rmsd=True,
                               compute_target_rmsd=True,
                               compute_hotspot_rmsd=True,
                               binder_chain="A",
                               target_chain="B"):
    res_dict = {}

    if compute_monomer_rmsd:
        binder_complex_monomer_rmsd = aligned_chain_rmsd_over_other_chain(
            reference_pdb=binder_complex_pdb,
            align_pdb=binder_monomer_pdb,
            reference_align_chain_id=binder_chain,
            align_align_chain_id=binder_chain,
            reference_rmsd_chain_id=binder_chain,
            align_rmsd_chain_id=binder_chain,
            seq_align_first=False,
            seq_align_second=False
        )
        res_dict["binder_complex_monomer_rmsd"] = binder_complex_monomer_rmsd

    if compute_target_rmsd:
        target_complex_backbone_rmsd = aligned_chain_rmsd_over_other_chain(
            reference_pdb=binder_backbone_pdb,
            align_pdb=binder_complex_pdb,
            reference_align_chain_id=target_chain,
            align_align_chain_id=target_chain,
            reference_rmsd_chain_id=target_chain,
            align_rmsd_chain_id=target_chain,
            seq_align_first=True,
            seq_align_second=True
        )
        res_dict["target_complex_backbone_rmsd"] = target_complex_backbone_rmsd

    if compute_hotspot_rmsd:
        binder_complex_backbone_target_aligned_rmsd = aligned_chain_rmsd_over_other_chain(
            reference_pdb=binder_backbone_pdb,
            align_pdb=binder_complex_pdb,
            reference_align_chain_id=target_chain,
            align_align_chain_id=target_chain,
            reference_rmsd_chain_id=binder_chain,
            align_rmsd_chain_id=binder_chain,
            seq_align_first=True,
            seq_align_second=False
        )
        res_dict["binder_complex_backbone_target_aligned_rmsd"] = binder_complex_backbone_target_aligned_rmsd

    return res_dict

# def compute_rmsd_table(df, pred_col, monomer_col, backbone_col,
#                        compute_monomer_rmsd=True,
#                        compute_target_rmsd=True,
#                        compute_hotspot_rmsd=True,
#                        binder_chain="A",
#                        target_chain="B"):
#     tqdm.tqdm.pandas()

#     def row_fn(x):
#         pred_pdb = x[pred_col]
#         backbone_pdb = x[backbone_col] if compute_target_rmsd or compute_hotspot_rmsd else None
#         monomer_pdb = x[monomer_col] if compute_monomer_rmsd else None
#         return compute_rmsd_between_pdbs(
#             pred_pdb, backbone_pdb, monomer_pdb,
#             compute_monomer_rmsd=compute_monomer_rmsd,
#             compute_target_rmsd=compute_target_rmsd,
#             compute_hotspot_rmsd=compute_hotspot_rmsd,
#             binder_chain=binder_chain,
#             target_chain=target_chain
#         )

#     rmsd_rows = df.progress_apply(row_fn, axis=1).tolist()
#     rmsd_df = pd.DataFrame(rmsd_rows)
#     merged_df = pd.concat([df, rmsd_df], axis=1)
#     return merged_df

def compute_rmsd_table(df, pred_col, monomer_col, backbone_col,
                       compute_monomer_rmsd=True,
                       compute_target_rmsd=True,
                       compute_hotspot_rmsd=True,
                       binder_chain="A",
                       target_chain="B"):
    tqdm.tqdm.pandas()

    # Define the keys we expect to see in the output to handle failures gracefully
    expected_keys = []
    if compute_monomer_rmsd: expected_keys.append("monomer_rmsd")
    if compute_target_rmsd: expected_keys.append("target_rmsd")
    if compute_hotspot_rmsd: expected_keys.append("hotspot_rmsd")

    def row_fn(x):
        try:
            pred_pdb = x[pred_col]
            # Handle potentially missing file paths or NaN values in the dataframe
            if pd.isna(pred_pdb):
                raise ValueError(f"Prediction PDB path is NaN for row {x.name}")

            backbone_pdb = x[backbone_col] if (compute_target_rmsd or compute_hotspot_rmsd) else None
            monomer_pdb = x[monomer_col] if compute_monomer_rmsd else None

            return compute_rmsd_between_pdbs(
                pred_pdb, backbone_pdb, monomer_pdb,
                compute_monomer_rmsd=compute_monomer_rmsd,
                compute_target_rmsd=compute_target_rmsd,
                compute_hotspot_rmsd=compute_hotspot_rmsd,
                binder_chain=binder_chain,
                target_chain=target_chain
            )
        except Exception as e:
            # Print the error for debugging, but don't stop the whole table
            print(f"Error processing row {x.name} ({x.get('id', 'unknown')}): {e}")
            # Return a dictionary of NaNs so the DataFrame remains the same length
            return {key: np.nan for key in expected_keys}

    # Apply the function across the rows
    rmsd_rows = df.progress_apply(row_fn, axis=1).tolist()
    
    # Create the RMSD dataframe
    rmsd_df = pd.DataFrame(rmsd_rows)
    
    # Ensure indices match before concat (crucial if df has been filtered previously)
    rmsd_df.index = df.index
    
    merged_df = pd.concat([df, rmsd_df], axis=1)
    return merged_df

def main():
    parser = argparse.ArgumentParser(description="Given a csv containing the binder and target PDBs, compute the RMSD of the binder and target backbone.")
    parser.add_argument("-csv", help="Path to the csv file containing the binder and target PDBs.")
    parser.add_argument("-pred_col", default="boltz1_pred_pdb_path", help="Column name for the prediction PDB (binder A + target B) file.")
    parser.add_argument("-monomer_col", default="af2_monomer_pdb_path", help="Column name for the monomer binder (binder A) file.")
    parser.add_argument("-backbone_col", default="rfd_backbone_path", help="Column name for the diffused backbone (binder A + target B) PDB file.")
    parser.add_argument("--binder_chain", default="A", help="Chain ID for the binder (default: A)")
    parser.add_argument("--target_chain", default="B", help="Chain ID for the target (default: B)")
    ## Options to specifcy which RMSD to compute
    parser.add_argument("--compute_monomer_rmsd", action="store_true", default=None,
                        help="Compute binder monomer RMSD (binder A vs monomer A).")
    parser.add_argument("--compute_target_rmsd", action="store_true", default=None,
                        help="Compute target backbone RMSD (target B vs complex B).")
    parser.add_argument("--compute_hotspot_rmsd", action="store_true", default=None,
                        help="Compute hotspot RMSD (binder A vs complex A after aligning target B).")
    parser.add_argument("--all", action="store_true", default=True,
                        help="Compute all types of RMSD. Ignored if any specific --compute_* flag is set.")

    parser.add_argument("-out", help="Path to the output CSV file. If not specified, output will be written to the same directory as the input CSV.")

    args = parser.parse_args()
    input_df = pd.read_csv(args.csv)

    if args.out:
        output_path = args.out
    else:
        input_dir = os.path.dirname(args.csv)
        output_path = os.path.join(input_dir, "rmsd_results.csv")

    # If any individual compute flag is explicitly set, disable --all
    specific_flags = [
        args.compute_monomer_rmsd,
        args.compute_target_rmsd,
        args.compute_hotspot_rmsd
    ]
    

    if any(flag is True for flag in specific_flags):
        args.all = False

    # If --all is still True, turn on all compute flags
    if args.all:
        args.compute_monomer_rmsd = True
        args.compute_target_rmsd = True
        args.compute_hotspot_rmsd = True

    print("calculating RMSD...")
    merged_df = compute_rmsd_table(
        input_df,
        pred_col=args.pred_col,
        monomer_col=args.monomer_col,
        backbone_col=args.backbone_col,
        compute_monomer_rmsd=args.compute_monomer_rmsd,
        compute_target_rmsd=args.compute_target_rmsd,
        compute_hotspot_rmsd=args.compute_hotspot_rmsd,
        binder_chain=args.binder_chain,
        target_chain=args.target_chain
    )
    merged_df.to_csv(output_path, index=False)
    print(f"RMSD table saved to {output_path}")

if __name__ == "__main__":
    main()