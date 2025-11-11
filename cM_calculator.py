import csv
import sys

if len(sys.argv) >= 3:
    map_file = sys.argv[1]
    cM_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) >= 4 else "map_with_cM.csv"
else:
    map_file = input("Enter the name of the map.csv file: ").strip()
    cM_file = input("Enter the name of the cM_lengths.csv file: ").strip()
    output_file = input("Enter a name for the output file (default: map_with_cM.csv): ").strip() or "map_with_cM.csv"

def normalize_headers(reader):
    reader.fieldnames = [h.strip().lower() for h in reader.fieldnames]
    return reader

cM_dict = {}
with open(cM_file, newline='', encoding='utf-8-sig') as f:
    reader = csv.DictReader(f)
    reader = normalize_headers(reader)
    for row in reader:
        chrom = row.get("chrom") or row.get("chr")
        chrom = chrom.strip()
        length = float(row["length"])
        cM = float(row["cm"])
        cM_dict[chrom] = {"length": length, "cM": cM}

# --- Process map.csv and compute cM values ---
with open(map_file, newline='', encoding='utf-8-sig') as infile, open(output_file, "w", newline='', encoding='utf-8') as outfile:
    reader = csv.DictReader(infile)
    reader = normalize_headers(reader)

    # Validate columns
    if not all(k in reader.fieldnames for k in ["marker", "pos"]) or not any(k in reader.fieldnames for k in ["chr", "chrom"]):
        raise KeyError("map.csv must contain columns 'marker', 'chr' (or 'chrom'), and 'pos'")

    fieldnames = ["marker", "chr", "pos", "cM"]
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()

    for row in reader:
        marker = row["marker"].strip()
        chrom = (row.get("chr") or row.get("chrom")).strip()
        pos = float(row["pos"])

        if chrom in cM_dict:
            length = cM_dict[chrom]["length"]
            max_cM = cM_dict[chrom]["cM"]
            cM_value = (pos / length) * max_cM
        else:
            cM_value = ""

        writer.writerow({
            "marker": marker,
            "chr": chrom,
            "pos": int(pos),
            "cM": round(cM_value, 6) if cM_value != "" else ""
        })

print(f"\nNew file '{output_file}' created successfully.")
