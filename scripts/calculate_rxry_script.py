import argparse
import os
import re
import numpy as np
import pandas as pd
from scipy.stats import norm
# function to calculate ratio of X to autosomes
def calculate_Rt(idxstats, total_ref, total_map):
    rts = []
    for i in range(len(idxstats)):
        rts.append((idxstats[i][1]/total_map)/(idxstats[i][0]/total_ref))
    return rts

# read master file
def read_master_file(master_file):
    with open(master_file, 'r') as file:
        idxstats_files = file.readlines()
    return [line.strip() for line in idxstats_files if line.strip()]

def extract_sample_id(file_path: str) -> str:
    """
    Returns the base filename (minus extension) as the sample ID.
    Example: 'my_sample_001.idxstats' -> 'my_sample_001'
    """
    base_name = os.path.basename(file_path)      # e.g. 'my_sample_001.idxstats'
    root, ext = os.path.splitext(base_name)      # ('my_sample_001', '.idxstats')
    return root

# Read in metadata file
def read_metadata(metadata_path):
    """Read metadata file. """
    return pd.read_csv(metadata_path, sep='\t')

# Find sample ID column in metadata
def find_sample_id_column(metadata):
    possible_column_names = [
        'sample-id', 'sampleid', 'sample id',
        'id', 'featureid', 'feature id', 'feature-id', 'Run', 'SRA'
    ]
    for col in possible_column_names:
        if col in metadata.columns:
            return col
    raise ValueError("No valid sample ID column found in metadata file. Expected one of: " + ", ".join(possible_column_names))

def main():
    parser = argparse.ArgumentParser(description='Calculate Rx and Ry')
    parser.add_argument('--scaffolds', dest="scaffold_ids_file", required=True, type=str, help='File containing scaffold IDs of interest')
    parser.add_argument('--metadata', required=True, help="Path to the QIIME2 metadata file")
    parser.add_argument('--master_file', required=True, help="Path to the master file with idxstats paths")
    parser.add_argument('--homogametic_id', dest="x_id", required=True, type=str, help='Specify scaffold ID for homogametic chromosome (eg. In XY sex determination systems, homogametic chromosome is X,while in ZW sex determination systems, homogametic chromosome is Z )')
    parser.add_argument('--heterogametic_id', dest="y_id", required=True, type=str, help='Specify scaffold ID for heterogametic chromosome (eg. In XY sex determination systems, heterogametic chromosome is Y,while in ZW sex determination systems, heterogametic chromosome is W)')
    parser.add_argument('--system', dest="system", required=True, type=str, choices=['XY', 'ZW'], help='Specify the sex determination system (XY or ZW)')
    parser.add_argument('--threshold', dest="threshold", required=True, type=float, help='Specify the threshold for the ratio')
    parser.add_argument('--output', dest="output_file", required=True, type=str, help='Output file to save the results')
    args = parser.parse_args()

    metadata = read_metadata(args.metadata)
    sample_id_col = find_sample_id_column(metadata)
    idxstats_files = read_master_file(args.master_file)

    results = []

    with open(args.scaffold_ids_file, 'r') as file:
        scaffold_ids = file.read().splitlines()

    for idxstats_file in idxstats_files:
        sample_id = extract_sample_id(os.path.basename(idxstats_file))
        idxstats = pd.read_table(idxstats_file, header=None, index_col=0)
        idxstats = idxstats.loc[scaffold_ids]

        sex_chromosomes = [args.x_id] if args.system == 'XY' else [args.x_id]
        autosomes = set(scaffold_ids) - {args.x_id, args.y_id} if args.system == 'XY' else [args.x_id, args.y_id]
        
        missing_chromosomes = [chrom for chrom in sex_chromosomes if chrom not in idxstats.index]
        missing_autosomes = [chrom for chrom in autosomes if chrom not in idxstats.index]

        if missing_chromosomes:
            results.append({
                'SCiMS sample ID': sample_id,
                'Status': f"Error: Essential chromosomes {missing_chromosomes} are missing."
            })
            continue
        if missing_autosomes:
            results.append({
                'SCiMS sample ID': sample_id,
                'Status': f"Error: autosomal chromosomes {missing_autosomes} are missing."
            })
            continue

        c1 = np.array(idxstats.iloc[:, 0], dtype=float)
        c2 = np.array(idxstats.iloc[:, 1], dtype=float)
        total_ref = np.sum(c1)
        total_map = np.sum(c2)

        Rt_values = calculate_Rt(idxstats.values[:, :2], total_ref, total_map)

        if args.system == 'XY':
            x_id = args.x_id
            y_id = args.y_id

            if x_id not in idxstats.index:
                raise ValueError(f"{x_id} not found in the filtered idxstats")

            x_index = idxstats.index.get_loc(x_id)
            y_index = idxstats.index.get_loc(y_id) if y_id in idxstats.index else None

            copy = [Rt_values[x] for x in range(len(Rt_values)) if x != x_index and (y_index is None or x != y_index)]
            tot = Rt_values[x_index] / np.array(copy)
            Rx = np.mean(tot)
            z_value = np.round(norm.ppf(1 - ( 1 - args.threshold)/2), 3)
            
            conf_interval = (np.std(tot) / np.sqrt(len(copy))) * z_value
            CI1 = Rx - conf_interval
            CI2 = Rx + conf_interval

            x_count = idxstats.loc[x_id].iloc[1]
            y_count = idxstats.loc[y_id].iloc[1]
            tot_y = x_count + y_count
            
            Ry = (1.0 * y_count) / tot_y
            conf_interval = z_value * np.sqrt((Ry * (1 - Ry)) / tot_y)
            CI1_y = Ry - conf_interval
            CI2_y = Ry + conf_interval

        else:
            z_id = args.x_id
            w_id = args.y_id

            if z_id not in idxstats.index:
                raise ValueError(f"{z_id} not found in the filtered idxstats")

            z_index = idxstats.index.get_loc(z_id)
            w_index = idxstats.index.get_loc(w_id) if w_id in idxstats.index else None

            copy = [Rt_values[x] for x in range(len(Rt_values)) if x != z_index and (w_index is None or x != w_index)]
            tot = Rt_values[z_index] / np.array(copy)
            Rz = np.mean(tot)
            z_value = np.round(norm.ppf(1 - ( 1 - args.threshold)/2), 3)

            conf_interval = (np.std(tot) / np.sqrt(len(copy))) * z_value
            CI1 = Rz - conf_interval
            CI2 = Rz + conf_interval

            z_count = idxstats.loc[z_id].iloc[1]
            w_count = idxstats.loc[w_id].iloc[1]
            tot = z_count + w_count
            Rw = (1.0 * w_count) / tot

            conf_interval = z_value * np.sqrt((Rw * (1 - Rw)) / tot)
            CI1_y = Rw - conf_interval
            CI2_y = Rw + conf_interval

        results.append({
            'SCiMS sample ID': sample_id,
            'Rx': np.round(Rx, 4) if args.system == 'XY' else np.round(Rz, 4),
            'Rx 95% CI': (np.round(CI1, 3), np.round(CI2, 3)),
            'Ry': np.round(Ry, 4) if args.system == 'XY' else np.round(Rw, 4),
            'Ry 95% CI': (np.round(CI1_y, 3), np.round(CI2_y, 3)),
        })

    results_df = pd.DataFrame(results)
    
    merged_df = pd.merge(metadata, results_df, left_on=sample_id_col, right_on='SCiMS sample ID', how='left')
    merged_df.to_csv(args.output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()