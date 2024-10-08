import os
import subprocess
import argparse

def run_blastx(protein_faa, gene_fna, output_file):
    """
    Run BLASTX between a protein and nucleotide file.
    
    Args:
        protein_faa (str): Path to the input protein file in .faa format.
        gene_fna (str): Path to the input gene file in .fna format.
        output_file (str): Path to the output file where BLAST results will be stored.
    """
    # BLASTX command: translates nucleotide sequences and searches against a protein db
    blastx_cmd = [
        "blastx",
        "-query", gene_fna,
        "-subject", protein_faa,
        "-out", output_file,
        "-outfmt", "6",  # Tabular output
        "-evalue", "1e-5",  # Set e-value threshold
        "-max_target_seqs", "1",  # Report only the best match
    ]
    
    # Run the BLASTX command
    try:
        subprocess.run(blastx_cmd, check=True)
        print(f"BLASTX completed successfully. Results saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running BLASTX: {e}")

def parse_blast_output(output_file):
    """
    Parse BLASTX output to extract matching genes.
    
    Args:
        output_file (str): Path to the BLAST output file.
    
    Returns:
        List of gene IDs.
    """
    genes = set()
    with open(output_file, 'r') as file:
        for line in file:
            columns = line.strip().split("\t")
            if len(columns) > 1:
                gene_id = columns[0]  # Assuming the first column is the gene ID from the nucleotide file
                genes.add(gene_id)
    
    return list(genes)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run BLASTX to find genes matching a protein sequence.")
    parser.add_argument("protein_faa", help="Input protein .faa file")
    parser.add_argument("gene_fna", help="Input gene .fna file")
    parser.add_argument("output_file", help="Output file for BLAST results")
    
    args = parser.parse_args()
    
    # Run BLASTX
    run_blastx(args.protein_faa, args.gene_fna, args.output_file)
    
    # Parse BLAST output and print matching genes
    genes = parse_blast_output(args.output_file)
    print(f"Genes present in {args.gene_fna}:")
    for gene in genes:
        print(gene)
    
    print("done!")

