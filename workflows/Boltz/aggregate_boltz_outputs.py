import os
import glob
import json
import re
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, PPBuilder
import tqdm
import argparse

def load_boltz(pred_dir,colab_input=None,output_format="pdb"):
    """
    Given the boltz prediction output directory, load all the availalbe metrics including:
    - metrics in the confidence json file (ex: iptm, plddt)
    - if generated, metrics from the ipsae file
    - if generated, rosetta metrics from running bindcraft_util
    """
    prediction_dir = pred_dir
    all_pdb_dirs = os.listdir(prediction_dir)
    results_df_rows = []
    confidence_rows = []
    pae_rows = []
    for pdb_dir_name in all_pdb_dirs:
        pdb_dir_path = os.path.join(prediction_dir,pdb_dir_name)
        all_pdbs_paths = glob.glob(os.path.join(pdb_dir_path,f"{pdb_dir_name}_model*.{output_format}"))
        ## Infer the number of samples from the number of pdbs
        for pdb_path in all_pdbs_paths:
            pdb_name = os.path.basename(pdb_path).replace(".pdb","")
            confidence_path = os.path.join(pdb_dir_path,f"confidence_{pdb_name}.json")
            pae_path = os.path.join(pdb_dir_path,f"pae_{pdb_name}.npz")
            pae_json_path = os.path.join(pdb_dir_path,f"pae_{pdb_name}.json") # This assumes the pae json was pre-extracted
            plddt_path = os.path.join(pdb_dir_path,f"plddt_{pdb_name}.npz")
            ## Extract sequences
            ## So slow
            #chain_a_seq,chain_b_seq = extract_boltz_seq(pdb_path)
            res_json  = {
                "prediction_name": pdb_dir_name,
                "pdb_model_name": pdb_name,
                "pdb_path": pdb_path,
                "confidence_path": confidence_path,
                "pae_path": pae_path,
                "plddt_path": plddt_path
                #"binder_seq":chain_a_seq,
                #"target_seq":chain_b_seq
            }
            results_df_rows.append(res_json)

            ## Add other metric data like confidence and pae
            extracted_confidence = extract_boltz_confidence(confidence_path)
            confidence_rows.append(extracted_confidence)
            if os.path.exists(pae_json_path):
                ## If pae is pre-extracted as json this is much faster
                extracted_pae = extract_boltz_pae(pae_json_path)
                pae_rows.append(extracted_pae)
            else:
                print("slow pae")
                ## Slow
                extracted_pae = extract_boltz_pae_npz(pdb_path,pae_path)
                pae_rows.append(extracted_pae)

    results_df = pd.DataFrame(results_df_rows)
    if len(confidence_rows)==0:
        print(f"Warning: No results found in directory {pred_dir}. Returning empty dataframe")
        return pd.DataFrame()
        
    confidence_df = pd.concat(confidence_rows,ignore_index=True)
    pae_df = pd.concat(pae_rows,ignore_index=True)
    merged_df = pd.concat([results_df,confidence_df,pae_df],axis=1)

    if colab_input:
        colab_input_df = pd.read_csv(colab_input)
        colab_input_df["chain_a_seq"] = colab_input_df["sequence"].apply(lambda x: x.split(":")[0])
        colab_input_df["chain_b_seq"] = colab_input_df["sequence"].apply(lambda x: x.split(":")[1])
        merged_df = merged_df.merge(colab_input_df[["id","chain_a_seq","chain_b_seq"]],left_on="prediction_name",right_on="id",how="left")
    
    ## Load ipsae outputs
    ipsae_df = load_ipsae_dir(prediction_dir)
    if len(ipsae_df) > 0:
        ipsae_df["pdb_model_name"] = ipsae_df["Model"].apply(lambda x: os.path.basename(x))
        merged_df = merged_df.merge(ipsae_df,on='pdb_model_name',how="left")
        
    ## Load rosetta/bindcraft util outputs
    combined_interface_df = load_interface_dir(prediction_dir)
    if len(combined_interface_df) > 0:
        combined_interface_df["pdb_model_name"] = combined_interface_df["binder_path"].apply(lambda x: os.path.basename(x).replace(".pdb",""))
        merged_df = merged_df.merge(combined_interface_df,on="pdb_model_name",how="left")
    return merged_df

def load_pdb_ipsae(pdb,max_only=True):
    """
    Given PDB filename, load its ipsae metrics assuming the ipsae files are available in the same directory
    """
    # If ipSAE was also run, get the related paths too
    pdb_name = os.path.basename(pdb)
    pdb_basename = pdb_name.replace(".pdb","")
    pdb_dir = os.path.dirname(pdb)
    # Search pattern for the corresponding ipsae score file with two variable numbers
    byres_search_pattern = os.path.join(pdb_dir, pdb_basename + "_*_*_byres.txt")
    matching_files = glob.glob(byres_search_pattern)
    regex_pattern = re.compile(rf"{re.escape(pdb_basename)}_(\d+)_(\d+)_byres\.txt")
    for f in matching_files:
        filename = os.path.basename(f)
        match = regex_pattern.match(filename)
        if match:
            pae_cutoff, dist_cutoff = match.groups()
            break
    else:
        print(f"No matching ipSAE byres file found for {pdb_name}. Returning empty df")
        return pd.DataFrame()
    ## Now infer the other two ipsae file name
    ipsae_pml = os.path.join(pdb_dir,f"{pdb_basename}_{pae_cutoff}_{dist_cutoff}.pml")
    ipsae_txt = os.path.join(pdb_dir,f"{pdb_basename}_{pae_cutoff}_{dist_cutoff}.txt")
    ## TODO: do some additional logic with the other two files
    ## For now just return the plain interchain ipsae metric
    ipsae_df = pd.read_csv(ipsae_txt,sep="\s+")
    if max_only:
        # Identify all unique chain pairs present in the file
        pairs = ipsae_df[ipsae_df["Type"] == "asym"][["Chn1", "Chn2"]].drop_duplicates()
        for _, row in pairs.iterrows():
            c1, c2 = row["Chn1"], row["Chn2"]
            # Extract Asymmetric ipSAE (Directional: Chn1 -> Chn2)
            asym_val = ipsae_df[
                (ipsae_df["Chn1"] == c1) & 
                (ipsae_df["Chn2"] == c2) & 
                (ipsae_df["Type"] == "asym")
            ]["ipSAE"]
            if not asym_val.empty:
                col_name_asym = f"Chn1_{c1}_to_Chn2_{c2}_ipSAE"
                ipsae_df[col_name_asym] = asym_val.item()
            # Extract Max ipSAE (Interaction between Chn1 and Chn2)
            max_val = ipsae_df[
                (ipsae_df["Chn1"] == c1) & 
                (ipsae_df["Chn2"] == c2) & 
                (ipsae_df["Type"] == "max")
            ]["ipSAE"]
            if not max_val.empty:
                col_name_max = f"Chn1_{c1}_to_Chn2_{c2}_ipSAE_max"
                ipsae_df[col_name_max] = max_val.item()
        # Collapse the dataframe to a single row containing all new columns
        # We keep only the global 'max' type row and drop duplicates
        ipsae_df = ipsae_df[ipsae_df["Type"] == "max"].copy().drop_duplicates("Type")
    # if max_only:
    #     # This collapse the ipsae metrics to a single number per prediction
    #     #ipsae_max_df = ipsae_df[ipsae_df["Type"]=="max"].copy()
    #     print(ipsae_df.columns)
    #     print(ipsae_df["Chn1"])
    #     print(ipsae_df["Chn2"])
    #     ## 2025-12-14 update:
    #     ## Also add the A->B and B->A ipSAE
    #     a_to_b_ipSAE = ipsae_df[(ipsae_df["Chn1"]=="A")&(ipsae_df["Chn2"]=="B")&(ipsae_df["Type"]=="asym")]["ipSAE"].item()
    #     b_to_a_ipSAE = ipsae_df[(ipsae_df["Chn1"]=="B")&(ipsae_df["Chn2"]=="A")&(ipsae_df["Type"]=="asym")]["ipSAE"].item()
    #     a_to_b_ipSAE_max = ipsae_df[(ipsae_df["Chn1"]=="A")&(ipsae_df["Chn2"]=="B")&(ipsae_df["Type"]=="max")]["ipSAE"].item()
    #     ipsae_df["Chn1_A_to_Chn2_B_ipSAE"] = a_to_b_ipSAE
    #     ipsae_df["Chn1_B_to_Chn2_A_ipSAE"] = b_to_a_ipSAE
    #     ipsae_df["Chn1_A_to_Chn2_B_ipSAE_max"] = a_to_b_ipSAE_max
    #     ## If we are dealing with >2 chains take the ipSAE outputs for the first pair 
    #     ipsae_df = ipsae_df[ipsae_df["Type"]=="max"].copy().drop_duplicates("max")
    return ipsae_df

def load_ipsae_dir(dir):
    """
    For each pdb in the given directory, combine the ipsae outputs
    """
    dir_pdbs = glob.glob(f"{dir}/**/*.pdb")
    combined_ipsae_df = pd.concat([load_pdb_ipsae(pdb,max_only=True) for pdb in dir_pdbs],ignore_index=True)
    return combined_ipsae_df

def load_interface_dir(dir):
    interface_dfs = glob.glob(f"{dir}/**/*.tsv")
    if len(interface_dfs) > 0:
        combined_bindcraft_interface_df = pd.concat([pd.read_csv(df,sep="\t") for df in interface_dfs],ignore_index=True)
        return combined_bindcraft_interface_df
    else:
        return pd.DataFrame()

def extract_boltz_confidence(path):
    with open(path) as f:
        boltz_confidence_json = json.load(f)
        return pd.json_normalize(boltz_confidence_json,sep="_")

def extract_boltz_pae(path):
    with open(path) as f:
        boltz_pae_json = json.load(f)
        return pd.json_normalize(boltz_pae_json,sep="_")

def extract_boltz_seq(pdb_path, chains=("A", "B")):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_path)
    ppb = PPBuilder()
    result = {}
    for model in structure:
        for chain in model:
            if chain.id in chains:
                peptides = ppb.build_peptides(chain)
                sequence = ''.join(str(pp.get_sequence()) for pp in peptides)
                result[chain.id] = sequence

    return result
    

def extract_boltz_pae_npz(pdb_path, pae_path):
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

def find_predictions_dirs(parent_dir):
    predictions_dirs = []
    # Walk through all subdirectories
    for root, dirs, files in os.walk(parent_dir):
        for d in dirs:
            if d == "predictions":
                predictions_dirs.append(os.path.join(root, d))
    return predictions_dirs

def load_all_boltz_dfs(pred_dir,colab_input=None):
    all_dirs_to_load = find_predictions_dirs(pred_dir)
    #print(all_dirs_to_load)
    all_dfs = [
        load_boltz(subdir,colab_input) for subdir in tqdm.tqdm(all_dirs_to_load)
    ]
    combined_df = pd.concat(all_dfs)
    return combined_df

def main():
    """
    Main function to parse command-line arguments and run the data loading.
    """
    parser = argparse.ArgumentParser(description="Load Boltz prediction metrics from a directory.")
    parser.add_argument("pred_dir", type=str,
                        help="The top-level directory containing the Boltz predictions.")
    parser.add_argument("--colab_input", type=str, default=None,
                        help="Path to the input CSV file from a Colab notebook.")
    parser.add_argument("--output_csv", type=str, default=None,
                        help="Path to save the output CSV file.")
    
    args = parser.parse_args()

    print(f"Loading data from directory: {args.pred_dir}")
    merged_df = load_all_boltz_dfs(args.pred_dir, args.colab_input)
    
    if not merged_df.empty:
        if args.output_csv:
            merged_df.to_csv(args.output_csv, index=False)
            print(f"Data saved to {args.output_csv}")
        else:
            print("Data loaded successfully. Here is the head of the resulting DataFrame:")
            print(merged_df.head())
    else:
        print("No data was loaded. The resulting DataFrame is empty.")

if __name__ == "__main__":
    main()
