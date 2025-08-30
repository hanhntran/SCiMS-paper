import argparse
import os
import re
import numpy as np
import pandas as pd
from scipy.stats import norm
from warnings import warn

# read idxstats files in a directory
def read_idxstats_files(idxstat_dir):
    idxstats_files = [os.path.join(idxstat_dir, f) for f in os.listdir(idxstat_dir) if f.endswith('.idxstats')]
    return idxstats_files

# Function to calculate Rt values
def calculate_Rt(idxstats, total_ref, total_map):
    rts = []
    for i in range(len(idxstats)):
        rts.append((idxstats[i][1] / total_map) / (idxstats[i][0] / total_ref))
    return np.array(rts)

# Function to calculate Rx and Ry values
def calculate_Rx_Ry(idxstats_files, scaffold_ids, x_id, y_id, system):
    rts = []
    with open(scaffold_ids, 'r') as f:
        scaffold_ids = f.read().splitlines()
    for idxstats_file in idxstats_files:
        sample_id, ext = os.path.splitext(os.path.basename(idxstats_file))
        idxstats = pd.read_table(idxstats_file, header=None, index_col=0)
        idxstats = idxstats.loc[scaffold_ids]

        sex_chromosomes = [x_id] if system == 'XY' else [x_id]
        autosomes = set(scaffold_ids) - {x_id, y_id} if system == 'XY' else [x_id, y_id]
        
        missing_chromosomes = [chrom for chrom in sex_chromosomes if chrom not in idxstats.index]
        missing_autosomes = [chrom for chrom in autosomes if chrom not in idxstats.index]

        if missing_chromosomes:
            rts.append({
                'Sample ID': sample_id,
                'Status': f"Error: Essential chromosomes {missing_chromosomes} are missing."
            })
            continue
        if missing_autosomes:
            rts.append({
                'Sample ID': sample_id,
                'Status': f"Error: autosomal chromosomes {missing_autosomes} are missing."
            })
            continue

        c1 = np.array(idxstats.iloc[:, 0], dtype=float)
        c2 = np.array(idxstats.iloc[:, 1], dtype=float)
        total_ref = np.sum(c1)
        total_map = np.sum(c2)
        Rt_values = calculate_Rt(idxstats.values[:, :2], total_ref, total_map)

        if system == 'XY':
            x_id = x_id
            y_id = y_id

            if x_id not in idxstats.index:
                raise ValueError(f"{x_id} not found in the filtered idxstats")
            x_index = idxstats.index.get_loc(x_id)
            y_index = idxstats.index.get_loc(y_id) if y_id in idxstats.index else None
            copy = [Rt_values[x] for x in range(len(Rt_values)) if x != x_index and (y_index is None or x != y_index)]
            tot = Rt_values[x_index] / np.array(copy)
            Rx = np.mean(tot)
            conf_interval = (np.std(tot) / np.sqrt(len(copy))) * 1.96
            CI1 = Rx - conf_interval
            CI2 = Rx + conf_interval

            x_count = idxstats.loc[x_id].iloc[1]
            y_count = idxstats.loc[y_id].iloc[1]
            tot_y = x_count + y_count
            Ry = (1.0 * y_count) / tot_y
            conf_interval = 1.96 * np.sqrt((Ry * (1 - Ry)) / tot_y)
            CI1_y = Ry - conf_interval
            CI2_y = Ry + conf_interval
            rts.append({
                'Sample': sample_id,
                'Rx': np.round(Rx, 3),
                'Rx 95% CI': f"{np.round(CI1, 3)}, {np.round(CI2, 3)}",
                'Ry': np.round(Ry, 3),
                'Ry 95% CI': f"{np.round(CI1_y, 3)}, {np.round(CI2_y, 3)}"
            })
        elif system == 'ZW':
            z_id = x_id
            w_id = y_id
            z_index = idxstats.index.get_loc(z_id)
            w_index = idxstats.index.get_loc(w_id)
            copy = [Rt_values[x] for x in range(len(Rt_values)) if x != z_index and x != w_index]
            Rz = np.mean(Rt_values[z_index] / copy)
            CI1 = Rz - 1.96 * np.std(Rt_values[z_index] / copy) / np.sqrt(len(copy))
            CI2 = Rz + 1.96 * np.std(Rt_values[z_index] / copy) / np.sqrt(len(copy))

            z_count = idxstats.loc[z_id].iloc[1]
            w_count = idxstats.loc[w_id].iloc[1]
            tot_w = z_count + w_count
            Rw = (1.0 * w_count) / tot_w
            conf_interval = 1.96 * np.sqrt((Rw * (1 - Rw)) / tot_w)
            CI1_w = Rw - conf_interval
            CI2_w = Rw + conf_interval
            rts.append({
                'Sample': sample_id,
                'Rz': np.round(Rz, 3),
                'Rz 95% CI': f"{np.round(CI1, 3)}, {np.round(CI2, 3)}",
                'Rw': np.round(Rw, 3),
                'Rw 95% CI': f"{np.round(CI1_w, 3)}, {np.round(CI2_w, 3)}"
            })
    return rts

def main():
    parser = argparse.ArgumentParser(description='Calculate Rx and Ry')
    parser.add_argument('--scaffolds', dest="scaffold_ids_file", required=True, type=str, help='File containing scaffold IDs of interest')
    parser.add_argument('--idxstats_dir', required=True, help="Path to idxstats files")
    parser.add_argument('--homogametic_id', dest="x_id", required=True, type=str, help='Specify scaffold ID for homogametic chromosome (eg. In XY sex determination systems, homogametic chromosome is X,while in ZW sex determination systems, homogametic chromosome is Z )')
    parser.add_argument('--heterogametic_id', dest="y_id", required=True, type=str, help='Specify scaffold ID for heterogametic chromosome (eg. In XY sex determination systems, heterogametic chromosome is Y,while in ZW sex determination systems, heterogametic chromosome is W)')
    parser.add_argument('--system', dest="system", required=True, type=str, choices=['XY', 'ZW'], help='Specify the sex determination system (XY or ZW)')
    parser.add_argument('--threshold', dest="threshold", required=True, type=float, help='Specify the threshold for the ratio')
    parser.add_argument('--output', dest="output_file", required=True, type=str, help='Output file to save the results')
    args = parser.parse_args()

    if args.system == 'XY':
        idxstats_files = read_idxstats_files(args.idxstats_dir)
        rts = calculate_Rx_Ry(idxstats_files, args.scaffold_ids_file, args.x_id, args.y_id, args.system)
        results_df = pd.DataFrame(rts)
        results_df.to_csv(args.output_file, index=False, sep="\t")
        print(f"Results saved to {args.output_file}")
    elif args.system == 'ZW':
        z_id = args.x_id
        w_id = args.y_id
        idxstats_files = read_idxstats_files(args.idxstats_dir)
        rts = calculate_Rx_Ry(idxstats_files, args.scaffold_ids_file, z_id, w_id, args.system)
        results_df = pd.DataFrame(rts)
        results_df.to_csv(args.output_file, index=False, sep="\t")
        print(f"Results saved to {args.output_file}")

if __name__ == "__main__":
    main()
    print("Done")