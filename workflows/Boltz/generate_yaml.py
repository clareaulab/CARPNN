import argparse
import sys
import os
import pandas as pd
import yaml

from Bio import PDB
from dataclasses import dataclass
from pathlib import Path

import io
import json
import shutil
from collections import defaultdict

@dataclass
class Seq:
    name: str
    sequence: str

MAX_MSA_SEQS = 16384
MAX_PAIRED_SEQS = 8192

## This is modified version of the msatojson.py script from Yoshitaka Moriwaki
## original script: 
## https://github.com/cddlab/alphafold3_tools/blob/main/alphafold3tools/msatojson.py

def get_msa_json(inputmsafile):
    """
    This function and its child functions are adopted from Yoshitaka Moriwaki's script 
    to extract MSA from af3 input
    https://github.com/cddlab/alphafold3_tools/blob/main/alphafold3tools/msatojson.py
    Write AlphaFold3 input JSON file from a3m-format MSA file.

    Args:
        inputmsafile (str): Input MSA file path.
        cardinality (int): The number of distinct polypeptide chains.
        stoichiometries (list[int]): Stoichiometries of each polypeptide chain.
        pairedmsas (list[list[Seq]]): Paired MSAs.
        unpairedmsas (list[list[Seq]]): Unpaired MSAs.
        outputfile (str): Output file path.
    """
    name = ""
    with open(inputmsafile, "r") as f:
        lines = f.readlines()
    residue_lens, stoichiometries = get_residuelens_stoichiometries(lines)
    if len(residue_lens) != len(stoichiometries):
        raise ValueError("Length of residue_lens and stoichiometries must be the same.")
    cardinality = len(residue_lens)
    print(
        f"The input MSA file contains {cardinality} distinct polypeptide chains."
    )
    print(f"Residue lengths: {residue_lens}")
    print(f"Stoichiometries: {stoichiometries}")
    pairedmsas, unpairedmsas = get_paired_and_unpaired_msa(
        lines, residue_lens, cardinality
    )
    content = generate_input_json_content(
        name=f"{name}",
        cardinality=cardinality,
        stoichiometries=stoichiometries,
        pairedmsas=pairedmsas,
        unpairedmsas=unpairedmsas,
    )
    return content
    #return pairedmsas, unpairedmsas

def int_id_to_str_id(num: int) -> str:
    """Encodes a number as a string, using reverse spreadsheet style naming.
    This block is cited from
    https://github.com/google-deepmind/alphafold3/blob/main/src/alphafold3/structure/mmcif.py#L40

    Args:
      num: A positive integer.

    Returns:
      A string that encodes the positive integer using reverse spreadsheet style,
      naming e.g. 1 = A, 2 = B, ..., 27 = AA, 28 = BA, 29 = CA, ... This is the
      usual way to encode chain IDs in mmCIF files.
    """
    if num <= 0:
        raise ValueError(f"Only positive integers allowed, got {num}.")

    num = num - 1  # 1-based indexing.
    output = []
    while num >= 0:
        output.append(chr(num % 26 + ord("A")))
        num = num // 26 - 1
    return "".join(output)

def convert_msas_to_str(msas):
    """convert MSAs to str format for AlphaFold3 input JSON file."""
    if msas == []:
        return ""
    else:
        return "\n".join(f"{seq.name}{seq.sequence}" for seq in msas) + "\n"

def split_a3msequences(residue_lens, line) -> list[str]:
    """Split a3m sequences into a list of a3m sequences.
    Note: The a3m-format MSA file represents inserted residues with lowercase.
    The first line (starting with '#') of the MSA file contains residue lengths
    and stoichiometries of each polypeptide chain.
    From the second line, the first sequence is the query.
    After this, the paired MSA blocks are followed by the unpaired MSA.
    Args:
        residue_lens: list[int]
            Residue lengths of each polypeptide chain
        line: str
            A3M sequences
    Returns:
        a3msequences: list[str]
            A3M sequences, len(a3msequences) should be the same as len(residue_lens).
    """
    a3msequences = [""] * len(residue_lens)
    i = 0
    count = 0
    current_residue = []

    for char in line:
        current_residue.append(char)
        if char == "-" or char.isupper():
            count += 1
        if count == residue_lens[i]:
            a3msequences[i] = "".join(current_residue)
            current_residue = []
            count = 0
            i += 1
            if i == len(residue_lens):
                break

    if current_residue and i < len(residue_lens):
        a3msequences[i] = "".join(current_residue)

    return a3msequences

def get_residuelens_stoichiometries(lines) -> tuple[list[int], list[int]]:
    """Get residue lengths and stoichiometries from msa file.
    Args:
        lines: list[str]
            Lines of input msa file
    Returns:
        residue_lens: list[int]
            Residue lengths of each polypeptide chain
        stoichiometries: list[int]
            Stoichiomerties of each polypeptide chain
    """
    if lines[0].startswith("#"):
        residue_lens_, stoichiometries_ = lines[0].split("\t")
        residue_lens = list(map(int, residue_lens_.lstrip("#").split(",")))
        stoichiometries = list(map(int, stoichiometries_.split(",")))
    else:
        # If the first line does not start with '#',
        # get the residue length from the first sequence.
        # Always assume a monomer prediction.
        if not lines[0].startswith(">"):
            raise ValueError(
                "The first line of the input MSA file must start with '#' or '>'."
            )
        residue_lens = [len(lines[1].strip())]
        stoichiometries = [1]
    return residue_lens, stoichiometries

def get_paired_and_unpaired_msa(
    lines: list[str], residue_lens: list[int], cardinality: int
) -> tuple[list[list[Seq]], list[list[Seq]]]:
    """Get paired and unpaired MSAs from input MSA file.
    Args:
        lines: list[str]
            Lines of input MSA file
        residue_lens: list[int]
            Residue lengths of each polypeptide chain
        cardinality: int
            Number of polypeptide chains
        query_seqnames: list[int]
            Query sequence names
    Returns:
        pairedmsas: list[list[Seq]]
            Paired MSAs, len(pairedmsa) should be the cardinality.
            If cardinality is 1, pairedmsas returns [[Seq("", "")]].
        unpairedmsas: list[list[Seq]]
            Unpaired MSAs, len(unpairedmsa) should be the cardinality.
    """
    pairedmsas: list[list[Seq]] = [[] for _ in range(cardinality)]
    unpairedmsas: list[list[Seq]] = [[] for _ in range(cardinality)]
    pairedflag = False
    unpairedflag = False
    seen = False
    seqnames_seen = []
    query_seqnames = [int(101 + i) for i in range(cardinality)]
    chain = -1
    start = 1 if lines[0].startswith("#") else 0
    for line in lines[start:]:
        if line.startswith(">"):
            if line not in seqnames_seen:
                seqnames_seen.append(line)
            else:
                seen = True
                continue
            if cardinality > 1 and line.startswith(
                ">" + "\t".join(map(str, query_seqnames)) + "\n"
            ):
                pairedflag = True
                unpairedflag = False
            elif any(line.startswith(f">{seq}\n") for seq in query_seqnames):
                pairedflag = False
                unpairedflag = True
                chain += 1
            seqname = line
        else:
            if seen:
                seen = False
                continue
            if pairedflag:
                a3mseqs = split_a3msequences(residue_lens, line)
                for i in range(cardinality):
                    pairedmsas[i].append(Seq(seqname, a3mseqs[i]))

            elif unpairedflag:
                a3mseqs = split_a3msequences(residue_lens, line)
                for i in range(cardinality):
                    # Remove all-gapped sequences
                    if a3mseqs[i] == "-" * residue_lens[i]:
                        continue
                    unpairedmsas[i].append(Seq(seqname, a3mseqs[i]))
            else:
                raise ValueError("Flag must be either paired or unpaired.")
    return pairedmsas, unpairedmsas

def generate_input_json_content(
    name: str,
    cardinality: int,
    stoichiometries: list[int],
    pairedmsas: list[list[Seq]],
    unpairedmsas: list[list[Seq]],
) -> str:
    """generate AlphaFold3 input JSON file.

    Args:
        name (str): Name of the protein complex.
                    Used for the name field in the JSON file.
        cardinality (int): The number of distinct polypeptide chains.
        stoichiometries (list[int]): Stoichiometries of each polypeptide chain.
        pairedmsas (list[list[Seq]]): Paired MSAs.
        unpairedmsas (list[list[Seq]]): Unpaired MSAs.
    Returns:
        str: JSON string for AlphaFold3 input file.
    """
    sequences: list[dict] = []
    chain_id_count = 0
    null = None
    for i in range(cardinality):
        # unpairedmsa[i][0] is more appropriate than pairedmsa[i][0].
        query_seq = unpairedmsas[i][0].sequence
        chain_ids = [
            int_id_to_str_id(chain_id_count + j + 1) for j in range(stoichiometries[i])
        ]
        chain_id_count += stoichiometries[i]
        sequences.append(
            {
                "protein": {
                    "id": chain_ids,
                    "sequence": query_seq,
                    "modifications": [],
                    "unpairedMsa": convert_msas_to_str(unpairedmsas[i]),
                    "pairedMsa": convert_msas_to_str(pairedmsas[i]),
                    "templates": [],
                }
            }
        )
    content = json.dumps(
        {
            "dialect": "alphafold3",
            "version": 1,
            "name": f"{name}",
            "sequences": sequences,
            "modelSeeds": [1],
            "bondedAtomPairs": null,
            "userCCD": null,
        },
        indent=4,
    )
    return content

def gather_msa_lines(x,a3m_files):
    '''
    Lifted a part of boltz1 processing logic to gather MSA lines from a3m files
    '''
    # process input x
    seqs = [x] if isinstance(x, str) else x
    N, REDO = 101, True
    # deduplicate and keep track of order
    seqs_unique = []
    [seqs_unique.append(x) for x in seqs if x not in seqs_unique]
    Ms = [N + seqs_unique.index(seq) for seq in seqs]
    
    # gather a3m lines
    a3m_lines = {}
    for a3m_file in a3m_files:
        update_M, M = True, None
        for line in open(a3m_file, "r"):
            if len(line) > 0:
                if "\x00" in line:
                    line = line.replace("\x00", "")
                    update_M = True
                if line.startswith(">") and update_M:
                    M = int(line[1:].rstrip())
                    update_M = False
                    if M not in a3m_lines:
                        a3m_lines[M] = []
                a3m_lines[M].append(line)

    a3m_lines = ["".join(a3m_lines[n]) for n in Ms]
    return a3m_lines


def convert_paired_msa(input_path, output_path):
    """
    Convert the paired MSA output from colabsearch into boltz friendly format
    """
    line_counts = 0
    idx_to_keep = 0
    pair_name = None
    entry_ids = None
    empty_char = "\x00"
    with open(input_path, 'r') as infile, open(output_path, 'a') as outfile:
        for line in infile:
            line = line.rstrip()
            ## For the first line get the protein ids
            if line_counts == 0:
                pair_name = line
                entry_ids = line[1:].split("\t")
                line_to_write = f">{entry_ids[idx_to_keep]}\n"
                outfile.write(line_to_write)
                line_counts += 1
            else:
                ## If we bump into another row that looks the same as the first row
                ## This is the section for a new protein
                ## increasement the index write the new id instead
                ## Also add an empty character to the line because boltz likes it for some reason
                if line == pair_name:
                    idx_to_keep += 1
                    outfile.write(f"{empty_char}>{entry_ids[idx_to_keep]}\n")
                    line_counts += 1
                else:
                    ## Even lines are uniref id names
                    if line_counts % 2 == 0:
                        uniref_id_to_use = line[1:].split("\t")[idx_to_keep]
                        line_to_write = f">{uniref_id_to_use}\n"
                        outfile.write(line_to_write)
                        line_counts += 1
                    else:
                        ## Odd lines are the actual alignment, no changes needed
                        outfile.write(f"{line}\n")
                        line_counts += 1
        outfile.write(empty_char)
        line_counts += 1
    print("Done converting paired MSA")

def convert_unpaired_msa(input_path, output_path):
    """
    Convert the unpaired MSA output from colabsearch into boltz friendly format
    """
    empty_char = "\x00"
    with open(input_path, 'r') as infile, open(output_path, 'a') as outfile:
        for line in infile:
            ## If the line starts with > and does not have \t it's a new protein
            ## and we need to add that weird little empty character
            if line.startswith(">") and "\t" not in line:
                outfile.write(empty_char+line)
            else:
                outfile.write(line)
        outfile.write(empty_char)
    print("Done converting paired MSA")


def convert_msa_to_boltz_csv(data, paired_msas, unpaired_msa, msa_dir):
    """
    Dump MSA files into boltz1 comptabitle csv formats
    data: a prediction name sequence pair
    e.g.{'EGF_EGFR_8hgs__AF-P06213-F1-model_v4_0': 'ECPLSHDGYCLHDGVCMYIEALDKYACNCVVGYIGERCQYRDLKWWE'}
    
    """
    for idx, name in enumerate(data):
        # Get paired sequences
        paired = paired_msas[idx].strip().splitlines()
        paired = paired[1::2]  # ignore headers
        paired = paired[: MAX_PAIRED_SEQS]

        # Set key per row and remove empty sequences
        keys = [idx for idx, s in enumerate(paired) if s != "-" * len(s)]
        paired = [s for s in paired if s != "-" * len(s)]

        # Combine paired-unpaired sequences
        unpaired = unpaired_msa[idx].strip().splitlines()
        unpaired = unpaired[1::2]
        unpaired = unpaired[: (MAX_MSA_SEQS - len(paired))]
        if paired:
            unpaired = unpaired[1:]  # ignore query is already present

        # Combine
        seqs = paired + unpaired
        keys = keys + [-1] * len(unpaired)

        # Dump MSA
        csv_str = ["key,sequence"] + [f"{key},{seq}" for key, seq in zip(keys, seqs)]
        msa_path = os.path.join(msa_dir, f"{name}.csv")
        with open(msa_path, "w") as f:
            f.write("\n".join(csv_str))
        #msa_path = os.path.join(msa_dir,f"{name}.csv")
        # msa_path = msa_dir / f"{name}.csv"
        # with msa_path.open("w") as f:
        #     f.write("\n".join(csv_str))

def generate_msa_dir(prediction_json,msa_path,pred_name,outdir):
    '''
    Prase the MSA file given and dump the outputs into the given directory 
    in a structured format that can be used by boltz1.
    '''
    ## Generate a directory to store the processed MSA files
    paired_dir = os.path.join(outdir, f"{pred_name}_paired")
    unpaired_dir = os.path.join(outdir, f"{pred_name}_unpaired")

    ## Remove the directory if alresady exists
    if os.path.exists(paired_dir):
        shutil.rmtree(paired_dir)
    if os.path.exists(unpaired_dir):
        shutil.rmtree(unpaired_dir)
    os.makedirs(paired_dir, exist_ok=True)
    os.makedirs(unpaired_dir, exist_ok=True)

    ## Extract information from MSA file
    content = json.loads(get_msa_json(msa_path))
    sequences = content["sequences"]

    ## Dump the information in the corresponding directory
    for item in sequences:
        unpaired_msa = item["protein"]["unpairedMsa"]
        paired_msa = item["protein"]["pairedMsa"]
        with open(os.path.join(unpaired_dir, f"unpaired.a3m"), 'a') as f:
            f.write(unpaired_msa)
        with open(os.path.join(paired_dir, f"paired.a3m"), 'a') as f:
            f.write(paired_msa)

    ## Do a little modification to paired MSA to make it match the output from mmseq server
    convert_unpaired_msa(os.path.join(unpaired_dir, f"unpaired.a3m"), os.path.join(unpaired_dir, f"unpaired_modified.a3m"))
    convert_paired_msa(os.path.join(paired_dir, f"paired.a3m"), os.path.join(paired_dir, f"paired_modified.a3m"))
    
    ## Conver the pair/unpaired a3m files into csv files readable by boltz
    unpaired_a3m_files = [os.path.join(unpaired_dir, f"unpaired_modified.a3m")]
    paired_a3m_files = [os.path.join(paired_dir, f"paired_modified.a3m")]
    unpaired_msas = gather_msa_lines(list(prediction_json.values()),unpaired_a3m_files)
    paired_msas = gather_msa_lines(list(prediction_json.values()),paired_a3m_files)
    convert_msa_to_boltz_csv(prediction_json, paired_msas, unpaired_msas, outdir)
    print("Conversion to MSA done")

## This is way how colabfold/search sanitized file names
def safe_filename(file: str) -> str:
    return "".join([c if c.isalnum() or c in ["_", ".", "-"] else "_" for c in file])

def colabsearch_to_boltz_yaml(colabsearch_df, msa_dir, outdir, behavior="strict"):
    '''
    Convert the filtered Foldseek results into a yaml file. Returns a dataframe
    containing the mapping between the binder and target analogs and the yaml file.
    '''
    yaml_df_rows = []
    for i, row in colabsearch_df.iterrows():
        row_id = row["id"]
        row_id_sanitized = safe_filename(row_id) # Colabsearch adjusts the input name 
        sequence = row["sequence"]
        ## Assuming __ as the separator for the binder and target names
        target_name = f"{row_id.split('__')[-1]}"
        binder_name = row_id.replace("__"+target_name,"")
        #binder_name = f"{row_id.split('__')[0]}"
        
        binder_seq = f"{sequence.split(':')[0]}"
        target_seq = f"{sequence.split(':')[1]}"
        msa_path = os.path.join(msa_dir, f"{row_id_sanitized}.a3m")
        row_json = {
            f"{row_id_sanitized}_0": binder_seq,
            f"{row_id_sanitized}_1": target_seq
        }

        ## check if the msa file exists
        if not os.path.exists(msa_path):
            raise FileNotFoundError(f"MSA file {msa_path} not found.")

        ## Process the corresponding MSA file
        generate_msa_dir(row_json, msa_path, row_id, msa_dir)

        # Construct the json that will be converted into boltz compatible yaml
        boltz_json = {
            "version": 1
        }
        binder_json = {
            "protein":{
                "id": "A",
                "sequence": binder_seq,
                "msa": f"{msa_dir}/{row_id_sanitized}_0.csv"
            }
        }
        target_json = {
            "protein":{
                "id": "B",
                "sequence": target_seq,
                "msa": f"{msa_dir}/{row_id_sanitized}_1.csv"
            }
        }
        boltz_json["sequences"] = [binder_json, target_json]
        ## Write the yaml file
        yaml_path = os.path.join(outdir, f"{row_id}.yaml")
        with open(yaml_path, 'w') as outfile:
            yaml.dump(boltz_json, outfile, default_flow_style=False)

        yaml_df_rows.append({
            "prediction_name": row_id,
            "prediction_name_sanitized": row_id_sanitized,
            "binder_name": binder_name,
            "binder_seq": binder_seq,
            "target_name": target_name,
            "target_seq": target_seq,
            "msa_path": msa_path,
            "yaml_path": yaml_path
        })
    yaml_df = pd.DataFrame(yaml_df_rows)
    return yaml_df


def colabsearch_to_boltz_yaml_monomer(colabsearch_df, msa_dir, outdir, behavior="strict"):
    '''
    Convert the colabfold input into a yaml file (monomer tasks only). Returns a dataframe
    containing the mapping between the binder and target analogs and the yaml file.
    '''
    yaml_df_rows = []
    for i, row in colabsearch_df.iterrows():
        row_id = row["id"]
        row_id_sanitized = safe_filename(row_id) # Colabsearch adjusts the input name 
        sequence = row["sequence"]
        msa_path = os.path.join(msa_dir, f"{i}.a3m")

        ## check if the msa file exists
        if not os.path.exists(msa_path):
            raise FileNotFoundError(f"MSA file {msa_path} not found.")

        ## For monomer the a3m can be directly used without further processing
        # Construct the json that will be converted into boltz compatible yaml
        boltz_json = {
            "version": 1
        }
        monomer_json = {
            "protein":{
                "id": "A",
                "sequence": sequence,
                "msa": msa_path
            }
        }
        boltz_json["sequences"] = [monomer_json]
        ## Write the yaml file
        yaml_path = os.path.join(outdir, f"{row_id}.yaml")
        with open(yaml_path, 'w') as outfile:
            yaml.dump(boltz_json, outfile, default_flow_style=False)

        yaml_df_rows.append({
            "prediction_name": row_id,
            "prediction_name_sanitized": row_id_sanitized,
            "monomer_seq": sequence,
            "msa_path": msa_path,
            "yaml_path": yaml_path
        })
    yaml_df = pd.DataFrame(yaml_df_rows)
    return yaml_df

def main():
    parser = argparse.ArgumentParser(description="Genereate YAML inputs based on local MSA outputs")
    parser.add_argument("-input", help="Path to csv used for colabsearch.")
    parser.add_argument("-msa", help="Directory where the MSAs (a3m) are generated")
    parser.add_argument("-out", help="Path to save the YAML results.")
    parser.add_argument("-monomer", action="store_true", help="If set, treat the input as monomer tasks.")
    
    args = parser.parse_args()

    colabsearch_df = pd.read_csv(args.input)
    yaml_outdir = os.path.join(args.out, "input_yamls")
    os.makedirs(yaml_outdir, exist_ok=True)
    print(args.monomer)

    # Generate the yaml files that will be used for boltz1 downsrtream
    if args.monomer:
        print("Generating YAML for monomer tasks")
        yaml_df = colabsearch_to_boltz_yaml_monomer(colabsearch_df, args.msa, yaml_outdir)
    else:
        print("Generating YAML for multimer tasks (currently only supports 2 chains)")
        yaml_df = colabsearch_to_boltz_yaml(colabsearch_df, args.msa, yaml_outdir)
    #yaml_df = colabsearch_to_boltz_yaml(colabsearch_df, args.msa, yaml_outdir)
    yaml_df.to_csv(os.path.join(args.out, "01_yaml_mapping.csv"), index=False)
    print("Done")
    

if __name__ == "__main__":
    main()