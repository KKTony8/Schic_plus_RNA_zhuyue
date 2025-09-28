# -*- coding: utf-8 -*-
import pandas as pd
import glob
import os
import numpy as np

# Input folder containing allValidPairs files
input_dir = "/data5/GPT/Wuky/Higashi/Higashi/demo_data/WKYthptxtnew/WKYthptxtnew"

# Output file path
output_file = "/data5/GPT/Wuky/Higashi/Higashi/demo_data/datathp.txt"

# Bin size (e.g., 500kb)
bin_size = 500000

# Get all allValidPairs files
all_files = sorted(glob.glob(os.path.join(input_dir, "*.allValidPairs")))

data_list = []

for cell_id, f in enumerate(all_files):
    print(f"Processing {f} -> cell_id {cell_id}")
    try:
        # Load columns: chr1, pos1, chr2, pos2
        df = pd.read_csv(f, sep="\t", header=None, usecols=[1, 2, 4, 5])
        df.columns = ["chrom1", "pos1", "chrom2", "pos2"]

        # Map positions to bins
        df["pos1"] = (df["pos1"] // bin_size) * bin_size
        df["pos2"] = (df["pos2"] // bin_size) * bin_size

        # Aggregate contact counts per bin pair
        df_grouped = df.groupby(
            ["chrom1", "pos1", "chrom2", "pos2"]
        ).size().reset_index(name="count")

        # Normalized count = count / sqrt(total_count_in_cell)
        total_count = df_grouped["count"].sum()
        df_grouped["normalized_count"] = df_grouped["count"] / np.sqrt(total_count)

        # Add cell_id
        df_grouped["cell_id"] = cell_id

        # Reorder columns
        df_grouped = df_grouped[
            ["cell_id", "chrom1", "chrom2", "pos1", "pos2", "count", "normalized_count"]
        ]
        data_list.append(df_grouped)

    except Exception as e:
        print(f"⚠️ Failed to read {f}, skipped. Error: {e}")

# Merge all cells
combined_df = pd.concat(data_list, ignore_index=True)

# Save in Higashi input format
combined_df.to_csv(output_file, sep="\t", index=False)
print(f"✅ Done! File saved to {output_file}")

