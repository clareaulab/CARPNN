import os
import shutil
import pandas as pd
from pathlib import Path
from math import ceil
import argparse

def split_directory_contents(source_dir, dest_dir, n_chunks, output_csv):
    '''
    Splits the contents of a directory into N chunks of subdirectories.
    Writes the mapping of original paths to new paths in a csv.
    '''
    source_dir = Path(source_dir)
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    all_items = sorted([item for item in source_dir.iterdir() if item.is_file()])
    chunk_size = ceil(len(all_items) / n_chunks)

    path_mappings = []

    for i in range(n_chunks):
        chunk_items = all_items[i*chunk_size:(i+1)*chunk_size]
        chunk_subdir = dest_dir / f"chunk_{i+1}"
        chunk_subdir.mkdir(parents=True, exist_ok=True)

        for item in chunk_items:
            new_path = chunk_subdir / item.name
            shutil.copy(str(item), new_path)
            path_mappings.append({
                "original_path": str(item.resolve()),
                "new_path": str(new_path.resolve())
            })

    df = pd.DataFrame(path_mappings)
    df.to_csv(output_csv, index=False)
    print(f"✅ Mapping saved to {output_csv}")

def main():
    parser = argparse.ArgumentParser(description="Split a directory into N chunks of subdirectories and write path mappings.")
    parser.add_argument("-src", type=str, required=True, help="Path to the source directory.")
    parser.add_argument("-dest", type=str, required=True, help="Destination base directory where chunk folders will be created.")
    parser.add_argument("-n", type=int, required=True, help="Number of chunks to split into.")
    parser.add_argument("-csv", type=str, default="path_mapping.csv", help="Output CSV path to save original -> new path mappings.")

    args = parser.parse_args()

    split_directory_contents(args.src, args.dest, args.n, args.csv)

if __name__ == "__main__":
    main()