# import yaml
# import os
# import re
# import argparse
# import pandas as pd
# from pathlib import Path
# import itertools

# # Custom representer for [B, C] style list formatting
# class FlowSeq(list): pass
# def flow_seq_representer(dumper, data):
#     return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)
# yaml.add_representer(FlowSeq, flow_seq_representer)

# def get_chain_ids(num_needed):
#     """Generates a list of chain IDs: A-Z, then AA, AB..."""
#     alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
#     generator = itertools.chain(
#         alphabet, 
#         (''.join(pair) for pair in itertools.product(alphabet, repeat=2))
#     )
#     return [next(generator) for _ in range(num_needed)]


# def split_a3m_to_csv(a3m_path, lengths):
#     """
#     Splits A3M into CSVs per monomer. 
#     Keys: 
#     - Query sequence is always 0.
#     - If multimer: Paired rows increment (1, 2...), Unpaired rows are -1.
#     - If monomer: All MSA rows are -1.
#     """
#     msa_dir = a3m_path.parent / a3m_path.stem
#     msa_dir.mkdir(parents=True, exist_ok=True)
    
#     with open(a3m_path, 'r') as f:
#         # Load and clean lines, skip comments/empty
#         lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]

#     if len(lines) < 2: return [], []

#     orig_ids = lines[0].lstrip('>').split()
#     csv_paths = []
    
#     # 1. Clean A3M insertions (lowercase) from the whole file to align width
#     query_seq = re.sub(r'[a-z]', '', lines[1])
#     msa_entries = []
#     for i in range(2, len(lines), 2):
#         clean_seq = re.sub(r'[a-z]', '', lines[i+1])
#         msa_entries.append((lines[i], clean_seq))

#     is_multimer_complex = len(lengths) > 1
#     start_pos = 0

#     for i, length in enumerate(lengths):
#         end_pos = start_pos + length
#         label = orig_ids[i] if i < len(orig_ids) else str(i+1)
#         out_file = msa_dir / f"{label}.csv"
        
#         rows = []
#         # --- QUERY RECORD (Always Key 0) ---
#         query_slice = query_seq[start_pos:end_pos].upper().replace('.', '-')
#         rows.append({"key": 0, "sequence": query_slice})
        
#         # --- MSA RECORDS ---
#         current_key = 1
#         for header, clean_seq in msa_entries:
#             sub_seq = clean_seq[start_pos:end_pos].upper().replace('.', '-')
            
#             # Skip if this monomer slice is only gaps
#             if not any(c.isalpha() for c in sub_seq):
#                 continue

#             if is_multimer_complex:
#                 # Logic for complexes:
#                 other_parts = clean_seq[:start_pos] + clean_seq[end_pos:]
#                 is_unpaired = all(c == '-' or c == '.' for c in other_parts)
                
#                 if is_unpaired:
#                     assigned_key = -1
#                 else:
#                     assigned_key = current_key
#                     current_key += 1
#             else:
#                 # Logic for Monomer: MSA entries are forced to -1
#                 assigned_key = -1
            
#             rows.append({"key": assigned_key, "sequence": sub_seq})
            
#         # Save CSV
#         df = pd.DataFrame(rows)
#         df.to_csv(out_file, index=False)
#         csv_paths.append(str(out_file.resolve()))
#         start_pos = end_pos
        
#     return csv_paths, orig_ids  


# def process_file(a3m_path, yaml_dir):
#     """Processes A3M, writes YAML, and handles missing # headers for monomers."""
#     a3m_path = Path(a3m_path)
#     with open(a3m_path, 'r') as f:
#         header_line = f.readline()
#         lines = f.readlines()
        
#     # Find the query sequence to determine length if header is missing
#     query_seq_full = ""
#     for line in lines:
#         if not line.startswith(">") and not line.startswith("#"):
#             query_seq_full = line.strip().replace("-", "")
#             break

#     # 1. Handle Header Parsing
#     match = re.search(r'#([\d,]+)\s+([\d,]+)', header_line)
    
#     if match:
#         # It's a complex (or a monomer with ColabFold-style header)
#         lengths = [int(x) for x in match.group(1).split(',')]
#         counts = [int(x) for x in match.group(2).split(',')]
#     else:
#         # It's a standard monomer without the # header
#         print(f"No # header found for {a3m_path.name}, assuming single monomer.")
#         lengths = [len(query_seq_full)]
#         counts = [1]
    
#     # 2. Split MSAs into CSVs (using the updated regex-sub logic for insertions)
#     individual_csvs, _ = split_a3m_to_csv(a3m_path, lengths)
    
#     total_chains = sum(counts)
#     all_alphabet_ids = get_chain_ids(total_chains)

#     entities = []
#     chain_ptr, seq_ptr = 0, 0
#     for i, length in enumerate(lengths):
#         num_copies = counts[i]
#         ids = all_alphabet_ids[chain_ptr : chain_ptr + num_copies]
        
#         # Determine the sequence slice for this entity
#         entity_seq = query_seq_full[seq_ptr : seq_ptr + length]
        
#         entities.append({
#             "protein": {
#                 "id": FlowSeq(ids) if len(ids) > 1 else ids[0],
#                 "sequence": entity_seq,
#                 "msa": individual_csvs[i]
#             }
#         })
#         seq_ptr += length
#         chain_ptr += num_copies

#     # 3. Write YAML
#     yaml_path = yaml_dir / f"{a3m_path.stem}.yaml"
#     with open(yaml_path, 'w') as f:
#         yaml.dump({"sequences": entities}, f, default_flow_style=False, sort_keys=False)

#     return {
#         "prediction_name": a3m_path.stem,
#         "a3m_path": str(a3m_path.resolve()),
#         "yaml_path": str(yaml_path.resolve())
#     }

# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument("input_dir", help="Input folder with .a3m files")
#     parser.add_argument("-o", "--output_dir", required=True, help="Main output folder")
#     args = parser.parse_args()

#     base_out = Path(args.output_dir)
#     yaml_out = base_out / "input_yamls"
#     yaml_out.mkdir(parents=True, exist_ok=True)

#     prediction_data = []
#     a3m_files = list(Path(args.input_dir).glob("*.a3m"))

#     for f in a3m_files:
#         try:
#             row = process_file(f, yaml_out)
#             prediction_data.append(row)
#             print(f"Success: {f.name}")
#         except Exception as e:
#             print(f"Error {f.name}: {e}")

#     if prediction_data:
#         df = pd.DataFrame(prediction_data)
#         table_path = base_out / "prediction_metadata.tsv"
#         df.to_csv(table_path, sep="\t", index=False)
#         print(f"\nSummary table created at: {table_path}")

# if __name__ == "__main__":
#     main()

import yaml
import os
import re
import argparse
import pandas as pd
from pathlib import Path
import itertools

# Custom representer for [B, C] style list formatting
class FlowSeq(list): pass
def flow_seq_representer(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)
yaml.add_representer(FlowSeq, flow_seq_representer)

def get_chain_ids(num_needed):
    """Generates a list of chain IDs: A-Z, then AA, AB..."""
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    generator = itertools.chain(
        alphabet, 
        (''.join(pair) for pair in itertools.product(alphabet, repeat=2))
    )
    return [next(generator) for _ in range(num_needed)]

def split_a3m_to_csv(a3m_path, lengths):
    """Splits A3M into CSVs. If monomer, ALL rows (including query) are key -1."""
    msa_dir = a3m_path.parent / a3m_path.stem
    msa_dir.mkdir(parents=True, exist_ok=True)
    
    with open(a3m_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    if len(lines) < 2: return [], []

    orig_ids = lines[0].lstrip('>').split()
    csv_paths = []
    
    # 1. Clean A3M insertions (lowercase)
    query_seq = re.sub(r'[a-z]', '', lines[1])
    msa_entries = []
    for i in range(2, len(lines), 2):
        clean_seq = re.sub(r'[a-z]', '', lines[i+1])
        msa_entries.append((lines[i], clean_seq))

    is_multimer_complex = len(lengths) > 1
    start_pos = 0

    for i, length in enumerate(lengths):
        end_pos = start_pos + length
        label = orig_ids[i] if i < len(orig_ids) else str(i+1)
        out_file = msa_dir / f"{label}.csv"
        
        rows = []
        
        # --- QUERY RECORD ---
        query_slice = query_seq[start_pos:end_pos].upper().replace('.', '-')
        
        # If monomer, even the first line is -1. If multimer, it's 0.
        query_key = 0 if is_multimer_complex else -1
        rows.append({"key": query_key, "sequence": query_slice})
        
        # --- MSA RECORDS ---
        current_key = 1
        for header, clean_seq in msa_entries:
            sub_seq = clean_seq[start_pos:end_pos].upper().replace('.', '-')
            if not any(c.isalpha() for c in sub_seq):
                continue

            if is_multimer_complex:
                other_parts = clean_seq[:start_pos] + clean_seq[end_pos:]
                is_unpaired = all(c == '-' or c == '.' for c in other_parts)
                if is_unpaired:
                    assigned_key = -1
                else:
                    assigned_key = current_key
                    current_key += 1
            else:
                # Monomer MSA entries are also -1
                assigned_key = -1
            
            rows.append({"key": assigned_key, "sequence": sub_seq})
            
        pd.DataFrame(rows).to_csv(out_file, index=False)
        csv_paths.append(str(out_file.resolve()))
        start_pos = end_pos
        
    return csv_paths, orig_ids

def process_file(a3m_path, yaml_dir, custom_name=None):
    """Processes A3M, writes YAML, and returns metadata row."""
    a3m_path = Path(a3m_path)
    stem_name = custom_name if custom_name else a3m_path.stem

    with open(a3m_path, 'r') as f:
        header_line = f.readline()
        lines = f.readlines()
        
    query_seq_full = ""
    for line in lines:
        if not line.startswith(">") and not line.startswith("#"):
            query_seq_full = line.strip().replace("-", "")
            break

    match = re.search(r'#([\d,]+)\s+([\d,]+)', header_line)
    if match:
        lengths = [int(x) for x in match.group(1).split(',')]
        counts = [int(x) for x in match.group(2).split(',')]
    else:
        lengths = [len(query_seq_full)]
        counts = [1]
    
    individual_csvs, _ = split_a3m_to_csv(a3m_path, lengths)
    total_chains = sum(counts)
    all_alphabet_ids = get_chain_ids(total_chains)

    entities = []
    chain_ptr, seq_ptr = 0, 0
    for i, length in enumerate(lengths):
        num_copies = counts[i]
        ids = all_alphabet_ids[chain_ptr : chain_ptr + num_copies]
        entities.append({
            "protein": {
                "id": FlowSeq(ids) if len(ids) > 1 else ids[0],
                "sequence": query_seq_full[seq_ptr : seq_ptr + length],
                "msa": individual_csvs[i]
            }
        })
        seq_ptr += length
        chain_ptr += num_copies

    yaml_path = yaml_dir / f"{stem_name}.yaml"
    with open(yaml_path, 'w') as f:
        yaml.dump({"sequences": entities}, f, default_flow_style=False, sort_keys=False)

    return {
        "prediction_name": stem_name,
        "a3m_path": str(a3m_path.resolve()),
        "yaml_path": str(yaml_path.resolve())
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", help="Input folder with .a3m files")
    parser.add_argument("-o", "--output_dir", required=True, help="Main output folder")
    parser.add_argument("-colabfold", type=str, help="Path to ColabFold input CSV for renaming")
    args = parser.parse_args()

    base_out = Path(args.output_dir)
    yaml_out = base_out / "input_yamls"
    yaml_out.mkdir(parents=True, exist_ok=True)

    # Gather A3M files
    # We sort numerically if they are 0.a3m, 1.a3m etc.
    a3m_files = sorted(list(Path(args.input_dir).glob("*.a3m")), 
                       key=lambda x: int(x.stem) if x.stem.isdigit() else x.stem)

    # Handle ColabFold renaming logic
    id_map = None
    if args.colabfold:
        cf_df = pd.read_csv(args.colabfold)
        if 'id' not in cf_df.columns:
            raise ValueError("ColabFold CSV must contain an 'id' column.")
        
        id_map = cf_df['id'].tolist()
        
        if len(id_map) != len(a3m_files):
            raise ValueError(f"Counts mismatch! CSV has {len(id_map)} IDs, but folder has {len(a3m_files)} A3M files.")
        print(f"ColabFold mode: Renaming {len(a3m_files)} files using CSV IDs.")

    prediction_data = []
    for i, f in enumerate(a3m_files):
        try:
            custom_name = id_map[i] if id_map else None
            row = process_file(f, yaml_out, custom_name=custom_name)
            prediction_data.append(row)
            print(f"Success: {f.name} -> {row['prediction_name']}")
        except Exception as e:
            print(f"Error {f.name}: {e}")

    if prediction_data:
        df = pd.DataFrame(prediction_data)
        table_path = base_out / "prediction_metadata.tsv"
        df.to_csv(table_path, sep="\t", index=False)
        print(f"\nSummary table created at: {table_path}")

if __name__ == "__main__":
    main()