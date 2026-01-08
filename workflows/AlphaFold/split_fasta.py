import os
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter # Added for single-line writing
from pathlib import Path
import math

def setup_cli():
    parser = argparse.ArgumentParser(description="Prepare input tables and chunked monomer FASTAs.")
    parser.add_argument("-i", "--input", required=True, help="Input file (Master .fasta or ColabFold .csv)")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory (will create input_tables inside)")
    parser.add_argument("-n", "--num_chunks", type=int, default=20, help="Number of chunks to split the sequences into")
    parser.add_argument("--no_sort", action="store_false", dest="sort_by_length", help="Disable sorting monomers by length")
    parser.set_defaults(sort_by_length=True)
    return parser.parse_args()

def load_data(input_path):
    """Loads sequences and returns a list of dictionaries."""
    path = Path(input_path)
    data = []

    if path.suffix.lower() in ['.fasta', '.fa']:
        print(f"Detected FASTA input: {path.name}")
        for record in SeqIO.parse(path, "fasta"):
            data.append({"id": record.id, "full_seq": str(record.seq)})
        
    elif path.suffix.lower() == '.csv':
        print(f"Detected CSV input: {path.name}")
        df = pd.read_csv(path)
        if 'id' not in df.columns or 'sequence' not in df.columns:
            raise ValueError("CSV must contain 'id' and 'sequence' columns.")
        for _, row in df.iterrows():
            data.append({"id": str(row['id']), "full_seq": str(row['sequence'])})
    else:
        raise ValueError("Unsupported file format. Please provide a .fasta or .csv file.")
    
    return data

def write_single_line_fasta(records, output_path):
    """Helper to write FASTA records with the sequence on exactly one line."""
    with open(output_path, "w") as handle:
        writer = FastaWriter(handle, wrap=None) # wrap=None prevents line wrapping
        writer.write_file(records)

def main():
    args = setup_cli()
    
    # 1. Setup Directory Structure
    base_dir = Path(args.outdir) / "input_tables"
    chunk_dir = base_dir / "fasta_chunked"
    
    base_dir.mkdir(parents=True, exist_ok=True)
    chunk_dir.mkdir(parents=True, exist_ok=True)
    
    # 2. Load and Process Data
    raw_data = load_data(args.input)
    if not raw_data:
        print("No sequences found.")
        return

    processed_items = []
    for item in raw_data:
        # Clean the ID: Take first item before "__"
        clean_id = item['id'].split('__')[0]
        
        # Extract Monomer: Take first item before ":"
        full_seq = item['full_seq']
        monomer_seq = full_seq.split(':')[0]
        
        processed_items.append({
            "id": clean_id,
            "original_id": item['id'],
            "full_seq": full_seq,
            "monomer_seq": monomer_seq,
            "monomer_len": len(monomer_seq)
        })

    # 3. Sorting Logic (Descending by Monomer Length)
    if args.sort_by_length:
        print("Sorting sequences by binder (monomer) length...")
        processed_items.sort(key=lambda x: x['monomer_len'], reverse=True)

    monomer_records = []
    full_csv_rows = []

    for item in processed_items:
        monomer_records.append(SeqRecord(
            Seq(item['monomer_seq']),
            id=item['id'],
            description=""
        ))
        full_csv_rows.append({"id": item['id'], "sequence": item['full_seq']})

    # 4. Create colabsearch_input.csv
    pd.DataFrame(full_csv_rows).to_csv(base_dir / "colabsearch_input.csv", index=False)
    
    # 5. Create binder_monomer.fasta (Single-line sequence)
    monomer_fasta_path = base_dir / "binder_monomer.fasta"
    write_single_line_fasta(monomer_records, monomer_fasta_path)
    
    # 6. Chunking Logic
    total_seqs = len(monomer_records)
    actual_chunks = min(args.num_chunks, total_seqs)
    chunk_size = math.ceil(total_seqs / actual_chunks)
    
    print(f"Splitting {total_seqs} sequences into {actual_chunks} chunks...")
    
    for i in range(actual_chunks):
        start_idx = i * chunk_size
        end_idx = min(start_idx + chunk_size, total_seqs)
        chunk_slice = monomer_records[start_idx:end_idx]
        
        if not chunk_slice:
            break
            
        chunk_filename = chunk_dir / f"part_{i+1}.fasta"
        write_single_line_fasta(chunk_slice, chunk_filename)

    print(f"Done! Outputs deposited in {base_dir}")

if __name__ == "__main__":
    main()

# import os
# import argparse
# import pandas as pd
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from pathlib import Path
# import math

# def setup_cli():
#     parser = argparse.ArgumentParser(description="Prepare input tables and chunked monomer FASTAs.")
#     parser.add_argument("-i", "--input", required=True, help="Input file (Master .fasta or ColabFold .csv)")
#     parser.add_argument("-o", "--outdir", required=True, help="Output directory (will create input_tables inside)")
#     parser.add_argument("-n", "--num_chunks", type=int, default=20, help="Number of chunks to split the sequences into")
#     parser.add_argument("--no_sort", action="store_false", dest="sort_by_length", help="Disable sorting monomers by length")
#     parser.set_defaults(sort_by_length=True)
#     return parser.parse_args()

# def load_data(input_path):
#     """Loads sequences and returns a list of dictionaries."""
#     path = Path(input_path)
#     data = []

#     if path.suffix.lower() in ['.fasta', '.fa']:
#         print(f"Detected FASTA input: {path.name}")
#         for record in SeqIO.parse(path, "fasta"):
#             data.append({"id": record.id, "full_seq": str(record.seq)})
        
#     elif path.suffix.lower() == '.csv':
#         print(f"Detected CSV input: {path.name}")
#         df = pd.read_csv(path)
#         if 'id' not in df.columns or 'sequence' not in df.columns:
#             raise ValueError("CSV must contain 'id' and 'sequence' columns.")
#         for _, row in df.iterrows():
#             data.append({"id": str(row['id']), "full_seq": str(row['sequence'])})
#     else:
#         raise ValueError("Unsupported file format. Please provide a .fasta or .csv file.")
    
#     return data

# def main():
#     args = setup_cli()
    
#     # 1. Setup Directory Structure
#     base_dir = Path(args.outdir) / "input_tables"
#     chunk_dir = base_dir / "fasta_chunked"
    
#     base_dir.mkdir(parents=True, exist_ok=True)
#     chunk_dir.mkdir(parents=True, exist_ok=True)
    
#     # 2. Load and Process Data
#     raw_data = load_data(args.input)
#     if not raw_data:
#         print("No sequences found.")
#         return

#     processed_items = []
#     for item in raw_data:
#         # Clean the ID: Take first item before "__"
#         clean_id = item['id'].split('__')[0]
        
#         # Extract Monomer: Take first item before ":"
#         full_seq = item['full_seq']
#         monomer_seq = full_seq.split(':')[0]
        
#         processed_items.append({
#             "id": clean_id,
#             "original_id": item['id'],
#             "full_seq": full_seq,
#             "monomer_seq": monomer_seq,
#             "monomer_len": len(monomer_seq)
#         })

#     # 3. Sorting Logic (Descending by Monomer Length)
#     if args.sort_by_length:
#         print("Sorting sequences by binder (monomer) length...")
#         processed_items.sort(key=lambda x: x['monomer_len'], reverse=True)

#     monomer_records = []
#     full_csv_rows = []

#     for item in processed_items:
#         # Record for FASTA files (Cleaned ID, Monomer sequence)
#         monomer_records.append(SeqRecord(
#             Seq(item['monomer_seq']),
#             id=item['id'],
#             description=""
#         ))
        
#         # Record for CSV (Full sequence for MSA generation)
#         # Note: We keep the cleaned ID here too for consistency across the pipeline
#         full_csv_rows.append({"id": item['id'], "sequence": item['full_seq']})

#     # 4. Create colabsearch_input.csv
#     pd.DataFrame(full_csv_rows).to_csv(base_dir / "colabsearch_input.csv", index=False)
    
#     # 5. Create binder_monomer.fasta
#     monomer_fasta_path = base_dir / "binder_monomer.fasta"
#     SeqIO.write(monomer_records, monomer_fasta_path, "fasta")
    
#     # 6. Chunking Logic
#     total_seqs = len(monomer_records)
#     actual_chunks = min(args.num_chunks, total_seqs)
#     chunk_size = math.ceil(total_seqs / actual_chunks)
    
#     print(f"Splitting {total_seqs} sequences into {actual_chunks} chunks...")
    
#     for i in range(actual_chunks):
#         start_idx = i * chunk_size
#         end_idx = min(start_idx + chunk_size, total_seqs)
#         chunk_slice = monomer_records[start_idx:end_idx]
        
#         if not chunk_slice:
#             break
            
#         chunk_filename = chunk_dir / f"part_{i+1}.fasta"
#         SeqIO.write(chunk_slice, chunk_filename, "fasta")

#     print(f"Done! Outputs deposited in {base_dir}")

# if __name__ == "__main__":
#     main()