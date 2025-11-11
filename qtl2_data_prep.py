#!/usr/bin/env python3
"""
qtl2_data_prep_final.py

Converts map.csv and geno_clean.csv into an R/qtl2-compatible CSV file,
following the listeria.csv format.

Adjustments made:
1. Phenotype (0/1) is the first data column (Column A).
2. Genotype codes (AA, AB, BB) are converted to A, H, B, or - (nan).
3. Dtype handling added for map.csv to suppress warnings.
"""

import sys
import csv
import pandas as pd
import numpy as np # Needed for NaN check

def convert_genotype(genotype):
    """Converts two-letter genotype to R/qtl2 format (A, H, B, -)."""
    if pd.isna(genotype) or not isinstance(genotype, str):
        return '-' # Handle NaN or non-string
        
    genotype = genotype.strip().upper()
    if genotype == 'AA':
        return 'A'
    elif genotype in ('AB', 'BA'):
        return 'H'
    elif genotype == 'BB':
        return 'B'
    else:
        return '-' # Catch any other unexpected value

def map_chromosome_to_number(chr_full_name):
    """Maps NC_036780.1 to 1, NC_036781.1 to 2, etc."""
    try:
        # Extract the number part: 'NC_036780.1' -> 36780
        number_part = int(chr_full_name.split('_')[1].split('.')[0])
        # Base number for chr 1 is 36780
        return number_part - 36780 + 1
    except:
        return np.nan # Use NaN for malformed names

def prepare_qtl2_csv(map_file, geno_file, output_file="qtl2_input.csv"):
    
    print(f"Reading map data from: {map_file}")
    print(f"Reading genotype data from: {geno_file}")
    
    # --- 1. Read and Process map.csv (with dtype fix) ---
    try:
        # Explicitly set dtypes for 'pos' and 'cM' to avoid DtypeWarning and ensure numeric types
        map_df = pd.read_csv(
            map_file,
            header=None, # Assuming no header row based on original example
            names=['marker_full', 'chr_full', 'pos', 'cM'],
            skiprows=1,
            dtype={'pos': 'object', 'cM': 'float64'}, # Read pos as object for now, then convert
            low_memory=False
        )
        # Convert pos after reading (to handle potential scientific notation if present)
        map_df['pos'] = pd.to_numeric(map_df['pos'], errors='coerce').astype('Int64')
        
    except FileNotFoundError:
        print(f"Error: Map file not found at {map_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading map file: {e}")
        sys.exit(1)

    # Extract marker ID (number after last underscore)
    map_df["real_marker"] = map_df["marker_full"].apply(lambda x: str(x).split("_")[-1])
    # Map full chromosome name to simple number (1-22)
    map_df["chr_num"] = map_df["chr_full"].apply(map_chromosome_to_number)
    
    # Set index to the full marker name for merging
    map_df = map_df.set_index('marker_full')

    # --- 2. Read and Process geno_clean.csv ---
    try:
        # Read genotype data, using the first column as marker index
        geno_df = pd.read_csv(geno_file, index_col=0)
        geno_df.index.name = 'marker_full' # Rename index for clean merge
        geno_df.columns = geno_df.columns.str.strip() # Clean up sample names
        
    except FileNotFoundError:
        print(f"Error: Genotype file not found at {geno_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading genotype file: {e}")
        sys.exit(1)

    # --- 3. Merge and Sort Data ---
    
    # Use the full marker names (index) to ensure alignment
    # Select markers that are present in both files and align them
    common_markers = map_df.index.intersection(geno_df.index)
    
    # Filter and combine relevant data
    map_subset = map_df.loc[common_markers, ["real_marker", "chr_num", "cM"]].copy()
    geno_subset = geno_df.loc[common_markers].copy()
    
    # Combine (this is implicitly an inner join on the index)
    combined_df = pd.concat([map_subset, geno_subset], axis=1)
    
    # Sort by chromosome and distance
    combined_df = combined_df.sort_values(by=["chr_num", "cM"])
    
    # Ordered lists of marker metadata
    marker_order = combined_df.index.tolist()
    marker_names = combined_df["real_marker"].tolist()
    chr_nums = combined_df["chr_num"].astype(int).astype(str).tolist() # Convert to string for output
    cm_values = combined_df["cM"].astype(str).tolist() # Convert to string for output

    # Ordered list of sample names
    sample_cols = geno_subset.columns.tolist()

    # --- 4. Determine Phenotype (0 = female, 1 = male) ---
    phenotype_map = {}
    for s in sample_cols:
        s_lower = s.lower().strip()
        if s_lower.endswith("f"):
            phenotype_map[s] = 0
        elif s_lower.endswith("m"):
            phenotype_map[s] = 1
        else:
            phenotype_map[s] = '' # Missing phenotype

    # --- 5. Genotype Conversion and Transposition ---
    
    # Select the genotype columns in the sorted marker order
    geno_matrix_raw = combined_df.loc[marker_order, sample_cols]
    
    # Transpose: Samples are now rows, markers are columns
    geno_matrix_t = geno_matrix_raw.T 
    
    # Apply the genotype conversion function to every cell in the transposed matrix
    geno_matrix_final = geno_matrix_t.applymap(convert_genotype)

    # --- 6. Final Output Assembly (R/qtl2 Format) ---
    
    # Initialize the list of final rows
    final_rows = []
    
    # Row 1: Marker names (preceded by 'id' or empty based on listeria.csv)
    # The listeria.csv file uses 'id' as the column header for the sample/phenotype column
    final_rows.append(['id'] + marker_names)
    
    # Row 2: Chromosome numbers (preceded by 'chr')
    final_rows.append(['chr'] + chr_nums)
    
    # Row 3: cM distances (preceded by 'pos')
    final_rows.append(['pos'] + cm_values)

    # Rows 4+: Individual samples
    for sample in sample_cols:
        # 1. Sample ID (e.g., MCYHF1_010_m) - This is the index of the row in the output file
        # 2. Phenotype (0 or 1) - This is the first data column (Column A as per your request)
        # 3. Genotypes (A, H, B, -) - Remaining columns
        
        # NOTE: The listeria.csv format uses the first data column as the trait/phenotype.
        # It is common practice to drop the sample ID after the header rows.
        # We will follow your instruction: Column A = Phenotype, Column B+ = Marker Data
        
        # Get the row of genotypes for the current sample
        genotypes = geno_matrix_final.loc[sample].tolist()
        
        # Append the final row: [Sample_ID (for compatibility), Phenotype, Genotypes...]
        # R/qtl2 usually expects a key column in the first data row. We use the Sample ID here
        # and then overwrite the header in the first cell of the data rows with the phenotype.
        final_row = [phenotype_map[sample]] + genotypes
        
        # For the R/qtl2 file, the first column must be a unique identifier for the row,
        # which is the Sample ID, followed by the traits/phenotypes.
        # The first column *header* is 'id', and the first *data* column is the first trait (phenotype).
        # We need to include the sample ID to have a full set of columns.
        
        final_row_with_id = [sample] + final_row
        final_rows.append(final_row_with_id)

    # Convert to DataFrame for easy CSV writing without headers
    output_df = pd.DataFrame(final_rows)
    
    # Write to CSV without index and without header
    output_df.to_csv(output_file, index=False, header=False)

    print(f"\n✅ R/qtl2 formatted CSV written to: {output_file}")
    print(f"Markers: {len(marker_names)} | Samples: {len(sample_cols)}")
    print("Structure: [Sample ID | Phenotype | Marker 1 Genotype | Marker 2 Genotype | ...]")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        # Use default names for ease of testing
        map_file = "map.csv"
        geno_file = "geno_clean.csv"
        output_file = "qtl2_input.csv"
        print("Using default filenames: map.csv, geno_clean.csv, qtl2_input.csv.")
        print("Usage: python qtl2_data_prep_final.py <map_file> <geno_file> [output_file]")
    else:
        map_file = sys.argv[1]
        geno_file = sys.argv[2]
        output_file = sys.argv[3] if len(sys.argv) > 3 else "qtl2_input.csv"

    prepare_qtl2_csv(map_file, geno_file, output_file)