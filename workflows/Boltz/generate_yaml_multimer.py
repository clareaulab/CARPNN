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

# def split_a3m(a3m_path, lengths):
#     """Splits concatenated A3M into isolated files named by original A3M IDs."""
#     msa_dir = a3m_path.parent / a3m_path.stem
#     msa_dir.mkdir(parents=True, exist_ok=True)
    
#     with open(a3m_path, 'r') as f:
#         lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]

#     if len(lines) < 2: return []

#     orig_ids = lines[0].lstrip('>').split()
#     msa_paths = []
    
#     start_pos = 0
#     for i, length in enumerate(lengths):
#         end_pos = start_pos + length
#         label = orig_ids[i] if i < len(orig_ids) else str(i+1)
#         out_file = msa_dir / f"{label}.a3m"
        
#         with open(out_file, 'w') as out_f:
#             out_f.write(f">{label}\n{lines[1][start_pos:end_pos]}\n")
#             for j in range(2, len(lines), 2):
#                 header = lines[j]
#                 full_seq = lines[j+1]
#                 sub_seq = full_seq[start_pos:end_pos]
#                 if any(c.isalpha() for c in sub_seq):
#                     out_f.write(f"{header}\n{sub_seq}\n")
        
#         msa_paths.append(str(out_file.resolve()))
#         start_pos = end_pos
        
#     return msa_paths, orig_ids

# def process_file(a3m_path, yaml_dir):
#     """Processes A3M, writes YAML to yaml_dir, and returns metadata row."""
#     a3m_path = Path(a3m_path)
#     with open(a3m_path, 'r') as f:
#         header_line = f.readline()
#         lines = f.readlines()
        
#         query_seq_full = ""
#         for line in lines:
#             if not line.startswith(">"):
#                 query_seq_full = line.strip().replace("-", "")
#                 break

#     match = re.search(r'#([\d,]+)\s+([\d,]+)', header_line)
#     if not match: raise ValueError(f"Header error in {a3m_path}")
        
#     lengths = [int(x) for x in match.group(1).split(',')]
#     counts = [int(x) for x in match.group(2).split(',')]
    
#     individual_msas, _ = split_a3m(a3m_path, lengths)
#     total_chains = sum(counts)
#     all_alphabet_ids = get_chain_ids(total_chains)

#     entities = []
#     chain_ptr, seq_ptr = 0, 0
#     for i, length in enumerate(lengths):
#         num_copies = counts[i]
#         ids = all_alphabet_ids[chain_ptr : chain_ptr + num_copies]
#         entities.append({
#             "protein": {
#                 "id": FlowSeq(ids) if len(ids) > 1 else ids[0],
#                 "sequence": query_seq_full[seq_ptr : seq_ptr + length],
#                 "msa": individual_msas[i]
#             }
#         })
#         seq_ptr += length
#         chain_ptr += num_copies

#     # Write YAML to the specific input_yamls directory
#     yaml_path = yaml_dir / f"{a3m_path.stem}.yaml"
#     with open(yaml_path, 'w') as f:
#         yaml.dump({"sequences": entities}, f, default_flow_style=False, sort_keys=False)

#     # Return metadata for the summary table
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

#     # Setup directories
#     base_out = Path(args.output_dir)
#     yaml_out = base_out / "input_yamls"
#     yaml_out.mkdir(parents=True, exist_ok=True)

#     prediction_data = []

#     # Process files
#     a3m_files = list(Path(args.input_dir).glob("*.a3m"))
#     print(f"Found {len(a3m_files)} a3m files. Processing...")

#     for f in a3m_files:
#         try:
#             row = process_file(f, yaml_out)
#             prediction_data.append(row)
#             print(f"Success: {f.name}")
#         except Exception as e:
#             print(f"Error {f.name}: {e}")

#     # Create and save the summary table
#     if prediction_data:
#         df = pd.DataFrame(prediction_data)
#         table_path = base_out / "prediction_metadata.tsv"
#         df.to_csv(table_path, sep="\t", index=False)
#         print(f"\nSummary table created at: {table_path}")
#     else:
#         print("\nNo metadata collected; table not created.")

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
    """Splits A3M into CSVs per monomer with pairing keys."""
    msa_dir = a3m_path.parent / a3m_path.stem
    msa_dir.mkdir(parents=True, exist_ok=True)
    
    with open(a3m_path, 'r') as f:
        # Load all lines, skip comments
        lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]

    if len(lines) < 2: return [], []

    orig_ids = lines[0].lstrip('>').split()
    csv_paths = []
    
    # Extract sequences and headers into pairs for easier processing
    # lines[0:1] is the first header/query, lines[2:] are MSA entries
    query_header = lines[0]
    query_seq = lines[1]
    msa_entries = []
    for i in range(2, len(lines), 2):
        msa_entries.append((lines[i], lines[i+1]))

    start_pos = 0
    for i, length in enumerate(lengths):
        end_pos = start_pos + length
        label = orig_ids[i] if i < len(orig_ids) else str(i+1)
        out_file = msa_dir / f"{label}.csv"
        
        rows = []
        
        # 1. The Query (Always Key 0)
        query_slice = query_seq[start_pos:end_pos].upper().replace('.', '-')
        rows.append({"key": 0, "sequence": query_slice})
        
        current_key = 1
        for header, full_seq in msa_entries:
            # Clean and slice
            full_seq_upper = full_seq.upper().replace('.', '-')
            sub_seq = full_seq_upper[start_pos:end_pos]
            
            # Logic to determine if this is a paired or unpaired alignment
            # In paired A3Ms, unpaired lines are padded with gaps for other chains.
            # We check if the ENTIRE alignment row is gaps EXCEPT for this monomer.
            other_parts = full_seq_upper[:start_pos] + full_seq_upper[end_pos:]
            is_unpaired = all(c == '-' for c in other_parts)
            
            # Determine the key
            # If the current sub_seq is just gaps, it doesn't belong in this monomer's CSV
            if not any(c.isalpha() for c in sub_seq):
                continue
                
            if is_unpaired:
                assigned_key = -1
            else:
                assigned_key = current_key
                current_key += 1
            
            rows.append({"key": assigned_key, "sequence": sub_seq})
            
        # Save as CSV
        df = pd.DataFrame(rows)
        df.to_csv(out_file, index=False)
        
        csv_paths.append(str(out_file.resolve()))
        start_pos = end_pos
        
    return csv_paths, orig_ids

def process_file(a3m_path, yaml_dir):
    """Processes A3M, writes YAML to yaml_dir, and returns metadata row."""
    a3m_path = Path(a3m_path)
    with open(a3m_path, 'r') as f:
        header_line = f.readline()
        # Read lines and skip header to find query sequence
        lines = f.readlines()
        query_seq_full = ""
        for line in lines:
            if not line.startswith(">") and not line.startswith("#"):
                query_seq_full = line.strip().replace("-", "")
                break

    match = re.search(r'#([\d,]+)\s+([\d,]+)', header_line)
    if not match: raise ValueError(f"Header error in {a3m_path}")
        
    lengths = [int(x) for x in match.group(1).split(',')]
    counts = [int(x) for x in match.group(2).split(',')]
    
    # Split MSAs into CSVs
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
                "msa": individual_csvs[i] # Now pointing to .csv
            }
        })
        seq_ptr += length
        chain_ptr += num_copies

    yaml_path = yaml_dir / f"{a3m_path.stem}.yaml"
    with open(yaml_path, 'w') as f:
        yaml.dump({"sequences": entities}, f, default_flow_style=False, sort_keys=False)

    return {
        "prediction_name": a3m_path.stem,
        "a3m_path": str(a3m_path.resolve()),
        "yaml_path": str(yaml_path.resolve())
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", help="Input folder with .a3m files")
    parser.add_argument("-o", "--output_dir", required=True, help="Main output folder")
    args = parser.parse_args()

    base_out = Path(args.output_dir)
    yaml_out = base_out / "input_yamls"
    yaml_out.mkdir(parents=True, exist_ok=True)

    prediction_data = []
    a3m_files = list(Path(args.input_dir).glob("*.a3m"))

    for f in a3m_files:
        try:
            row = process_file(f, yaml_out)
            prediction_data.append(row)
            print(f"Success: {f.name}")
        except Exception as e:
            print(f"Error {f.name}: {e}")

    if prediction_data:
        df = pd.DataFrame(prediction_data)
        table_path = base_out / "prediction_metadata.tsv"
        df.to_csv(table_path, sep="\t", index=False)
        print(f"\nSummary table created at: {table_path}")

if __name__ == "__main__":
    main()