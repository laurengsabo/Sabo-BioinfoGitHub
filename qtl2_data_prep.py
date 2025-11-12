#!/usr/bin/env python3
"""
qtl2_data_prep_final_v2.py

Converts map.csv and geno_clean.csv into an R/qtl2-compatible CSV file.

Adjustments:
1. Phenotype (0/1) is the absolute first data column (Column A).
2. Genotype codes (AA, AB, BB) are converted to A, H, B, or -.
3. Optional argument to subsample markers per chromosome.
"""

import sys
import csv
import pandas as pd
import numpy as np
import argparse # Use argparse for cleaner handling of optional parameters

def convert_genotype(genotype):
    """Converts two-letter genotype to R/qtl2 format (A, H, B, -)."""
    if pd.isna(genotype) or not isinstance(genotype, str):
        return '-'
        
    genotype = genotype.strip().upper()
    if genotype == 'AA':
        return 'A'
    elif genotype in ('AB', 'BA'):
        return 'H'
    elif genotype == 'BB':
        return 'B'
    else:
        return '-'

def map_chromosome_to_number(chr_full_name):
    """Maps NC_036780.1 to 1, NC_036781.1 to 2, etc."""
    try:
        # Extract the number part: 'NC_036780.1' -> 36780
        number_part = int(chr_full_name.split('_')[1].split('.')[0])
        # Base number for chr 1 is 36780
        return number_part - 36780 + 1
    except:
        return np.nan

def subsample_markers(combined_df, max_markers_per_chr=None):
    """
    Subsamples markers to a maximum amount per chromosome,
    selecting markers that are equally spaced on the cM map.
    """
    if max_markers_per_chr is None:
        return combined_df
    
    # Ensure chr_num is ready for grouping
    combined_df['chr_num'] = combined_df['chr_num'].astype(int)
    
    subsampled_markers = []
    
    # Iterate over each chromosome group
    for chr_num, group in combined_df.groupby('chr_num'):
        num_markers = len(group)
        
        if num_markers <= max_markers_per_chr:
            # If the group is small enough, take all markers
            subsampled_markers.append(group)
        else:
            # Otherwise, select equally spaced markers
            # Create indices for the desired number of equally spaced markers
            # np.linspace provides the indices, np.round converts them to integers
            indices = np.round(np.linspace(0, num_markers - 1, max_markers_per_chr)).astype(int)
            
            # Select the markers at these indices
            subsampled_group = group.iloc[indices]
            subsampled_markers.append(subsampled_group)

    # Recombine the subsampled groups
    if not subsampled_markers:
        return pd.DataFrame() # Return empty if no data
        
    return pd.concat(subsampled_markers).sort_values(by=["chr_num", "cM"])


def prepare_qtl2_csv(map_file, geno_file, output_file="qtl2_inputfull.csv", max_markers_per_chr=None):
    
    print(f"Reading map data from: {map_file}")
    print(f"Reading genotype data from: {geno_file}")
    if max_markers_per_chr:
         print(f"Sampling mode: Limiting to {max_markers_per_chr} markers per chromosome.")
    
    # 1. Read and Process map.csv (with dtype fix and skiprows)
    try:
        # Skip header row (row 1) using skiprows=1
        map_df = pd.read_csv(
            map_file,
            header=None,
            names=['marker_full', 'chr_full', 'pos', 'cM'],
            skiprows=1, # Fix for the 'cM' string conversion error
            dtype={'pos': 'object', 'cM': 'float64'},
            low_memory=False
        )
        map_df['pos'] = pd.to_numeric(map_df['pos'], errors='coerce').astype('Int64')
    except Exception as e:
        print(f"Error reading map file: {e}")
        sys.exit(1)

    # Extract marker ID and map chromosome number
    map_df["real_marker"] = map_df["marker_full"].apply(lambda x: str(x).split("_")[-1])
    map_df["chr_num"] = map_df["chr_full"].apply(map_chromosome_to_number)
    map_df = map_df.set_index('marker_full')

    # 2. Read and Process geno_clean.csv
    try:
        geno_df = pd.read_csv(geno_file, index_col=0)
        geno_df.index.name = 'marker_full'
        geno_df.columns = geno_df.columns.str.strip()
    except Exception as e:
        print(f"Error reading genotype file: {e}")
        sys.exit(1)

    # 3. Merge, Sort, and Subsample Data
    
    common_markers = map_df.index.intersection(geno_df.index)
    map_subset = map_df.loc[common_markers, ["real_marker", "chr_num", "cM"]].copy()
    geno_subset = geno_df.loc[common_markers].copy()
    
    combined_df = pd.concat([map_subset, geno_subset], axis=1)
    
    # Apply subsampling logic
    combined_df = subsample_markers(combined_df, max_markers_per_chr)
    
    # Ordered lists of marker metadata
    marker_order = combined_df.index.tolist()
    marker_names = combined_df["real_marker"].tolist()
    chr_nums = combined_df["chr_num"].astype(int).astype(str).tolist()
    cm_values = combined_df["cM"].astype(str).tolist()

    # Ordered list of sample names
    sample_cols = geno_subset.columns.tolist()

    # 4. Determine Phenotype (0 = female, 1 = male) 
    phenotype_map = {}
    for s in sample_cols:
        s_lower = s.lower().strip()
        if s_lower.endswith("f"):
            phenotype_map[s] = 0
        elif s_lower.endswith("m"):
            phenotype_map[s] = 1
        else:
            phenotype_map[s] = '' # Missing phenotype

    # 5. Genotype Conversion and Transposition 
    
    # Select the genotype columns in the sorted/subsampled marker order
    geno_matrix_raw = combined_df.loc[marker_order, sample_cols]
    
    # Transpose: Samples are rows, markers are columns
    geno_matrix_t = geno_matrix_raw.T 
    
    # Apply the genotype conversion
    geno_matrix_final = geno_matrix_t.applymap(convert_genotype)

    # 6. Final Output Assembly (R/qtl2 Format)
    
    final_rows = []
    
    # Row 1: Header for all columns (Phenotype Name + Marker Names)
    # The first column *must* be the phenotype name for the R/qtl2 data file
    final_rows.append(['sex_pheno'] + marker_names)
    
    # Row 2: Chromosome numbers (preceded by 'chr')
    final_rows.append(chr_nums)
    
    # Row 3: cM distances (preceded by 'pos')
    final_rows.append(cm_values)

    # Rows 4+: Individual samples (Phenotype, then Genotypes)
    for sample in sample_cols:
        # 1. Phenotype (0 or 1) - THIS IS COLUMN A (first column of data rows)
        # 2. Genotypes (A, H, B, -) - Remaining columns
        
        genotypes = geno_matrix_final.loc[sample].tolist()
        
        # Final row structure: [Phenotype, Marker 1 Genotype, Marker 2 Genotype, ...]
        final_row = [phenotype_map[sample]] + genotypes
        final_rows.append(final_row)

    # Write to CSV without index and without header
    output_df = pd.DataFrame(final_rows)
    output_df.to_csv(output_file, index=False, header=False)

    print(f"\nR/qtl2 formatted CSV written to: {output_file}")
    print(f"Total Markers in Output: {len(marker_names)} | Total Samples: {len(sample_cols)}")


if __name__ == "__main__":
    # Setup argument parser to handle optional max_markers parameter
    parser = argparse.ArgumentParser(
        description="Prepare R/qtl2 CSV from map and genotype data.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("map_file", help="Path to the map CSV file.")
    parser.add_argument("geno_file", help="Path to the genotype CSV file.")
    parser.add_argument("output_file", nargs='?', default="qtl2_input.csv", 
                        help="Optional: Output file path (default: qtl2_input.csv).")
    parser.add_argument("--max-markers", type=int, default=None,
                        help="Optional: Maximum number of markers to select per chromosome for subsampling.")
    
    args = parser.parse_args()

    prepare_qtl2_csv(
        map_file=args.map_file, 
        geno_file=args.geno_file, 
        output_file=args.output_file,
        max_markers_per_chr=args.max_markers
    )