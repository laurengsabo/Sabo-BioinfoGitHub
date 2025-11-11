import csv

# Purpose: Input geno.csv and map.csv files for the creation of an annotated genotype matrix for r/qtl2 processing.

# Inputs
geno_file = input("Enter genotype CSV filename (e.g., geno_clean.csv): ")
map_file = input("Enter map CSV filename (e.g., map.csv): ")
output_file = input("Enter output filename (e.g., combined.csv): ")

# Load map data
map_data = {}
with open(map_file, newline='') as map_csv:
    reader = csv.DictReader(map_csv)
    for row in reader:
        marker = row.get('marker', '').strip()
        if marker:
            chrom = row.get('chr', '').replace('NC_03678', '').replace('.1', '').strip()
            map_data[marker] = {
                'chr': chrom,
                'cM': row.get('cM', '').strip()
            }

# Read geno.csv
with open(geno_file, newline='') as geno_csv:
    reader = csv.reader(geno_csv)
    header = next(reader)

    # Extract sample names (ignore first column)
    samples = header[1:]

    # Determine phenotype for each sample
    phenotype_map = {}
    for sample in samples:
        if sample.endswith('f'):
            phenotype_map[sample] = '1'
        elif sample.endswith('m'):
            phenotype_map[sample] = '0'
        else:
            phenotype_map[sample] = '-'

    # Prepare output file
    with open(output_file, 'w', newline='') as out_csv:
        writer = csv.writer(out_csv)
        writer.writerow(['phenotype', 'marker', 'chr', 'cM'] + samples)

        # Process genotype data
        for row in reader:
            marker_full = row[0].strip()
            geno_values = row[1:]

            # Extract marker base (before underscore)
            marker_parts = marker_full.split('_')
            base_marker = '_'.join(marker_parts[:2])  # e.g. "NC_036780.1_99202"
            chrom = base_marker.split('_')[0]
            chrom_num = chrom.replace('NC_03678', '').replace('.1', '')

            cM = map_data.get(base_marker, {}).get('cM', '-')

            # Convert genotypes
            converted_genotypes = []
            for g in geno_values:
                g = g.strip().upper() if g else ''
                if g == 'AA':
                    converted_genotypes.append('A')
                elif g in ('AB', 'BA'):
                    converted_genotypes.append('H')
                elif g == 'BB':
                    converted_genotypes.append('B')
                else:
                    converted_genotypes.append('-')

            # Phenotype: average of all samples (but you want one per row)
            # Let's just take the first sample's phenotype for simplicity:
            first_pheno = phenotype_map[samples[0]]

            # Write row
            writer.writerow([first_pheno, base_marker, chrom_num, cM] + converted_genotypes)

print(f"\n File '{output_file}' created successfully.")
