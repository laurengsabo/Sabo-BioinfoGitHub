#!/usr/bin/env python3
"""
Description: Script to split FASTA contigs and process them using multiprocessing.
@author: Lauren Sabo | @version: 2.0
Notes for: Nikesh Kumar
"""

import os
import sys
import subprocess as sp
    # Using "from multiprocessing import Process" you can type "Process" into the splitJobs Function below, 
        # if you were to just type "import multiprocessing", you would have to type "multiprocessing.Process(...)" in the function instead.
from multiprocessing import Process
from pyfaidx import Fasta

def splice(fasta_loc, size):
    """
    1. Filter out 'NC' sequences (Chromosomes) to keep only contigs.
    2. Chunk the list into sub-lists of a specific 'run size'.
    """
    fasta = Fasta(fastaLoc)
    # 1. Filtering
    contig_keys = [c for c in fasta.keys() if not c.startswith('NC')]

    # 2 & 3. Chunking (List Comprehension)
    return [contig_keys[i : i + int(size)] for i in range(0, len(contig_keys), int(size))]
    
def make_file(contig, output_dir):
    """For each contig, creates an empty file named after the contig within the chosen dir."""
    # Using os.path.join for cross-platform safety
    dest_path = os.path.join(output_dir, contig)
    sp.call(["touch", dest_path])

### SPLIT JOBS FUNCTION
    # 1. The function starts off by using the SPLICE Function and saving the lists of lists into a variable called "scopes".
    #       Remember: Each of the lists' lengths within "scopes" are equal to the run size you inputted
    # 2. A directory is then made with your desired name in your desired location.
    # 3. Now... (here comes the multiprocessing magic)
    #       For each of the lists within the "scopes" list, we are going to tell the computer to run the contigs within 
    #           each list, simultaneously. For example, if each list within "scopes" contains a list of length 5, then 5 
    #           contigs will run at once, and then the next group will run together, and so on.
    #       A. To do this, we must create a new list (AKA "jobs") with our scopes's lists + each of the list's contigs and their 
    #           commands. For example, if our scope is currently [A,B,C], then the altered list (AKA "jobs") 
    #           will be [do(A), do(B), do(C)]. We have to make a new list of "scopes" and not alter the current one. It's 
    #           simpler.
    #       B. Once we have successfully copied over the "jobs" list with all of the inner lists' contigs + the contigs' commands, 
    #           now we run it. Since it is a for-loop, we're going to run each of the scopes sequentially, and the n-number
    #           of items within each scope will run together.
    #       
def split_jobs(size, fasta_loc, destination_root, folder_name):
    """
    Orchestrates the multiprocessing logic.
    Groups jobs into 'batches' based on size to avoid overwhelming the CPU.
    """
    # 1. Get the chunked list of contigs
    scopes = splice(fasta_loc, size)

    # 2. Setup output directory using the 'os' module best practices
    output_path = os.path.join(destination_root, folder_name)
    os.makedirs(output_path, exist_ok=True)

    # 3. Multiprocessing Magic
    for scope in scopes:
        jobs = []
        # Start a batch of processes
        for contig in scope:
            p = Process(target=make_file, args=(contig, output_path))
            jobs.append(p)
            p.start()

        # Wait for this specific batch to finish before starting the next
        # This prevents the computer from trying to run 10,000 things at once
        for j in jobs:
            j.join()
            

def main():
    # Check if the user provided all 4 required arguments
    if len(sys.argv) != 5:
        print("Usage: python script.py <run_size> <fasta_path> <output_root> <project_name>")
        print("Example: python script.py 22 ./data.fna ./results my_experiment")
        sys.exit(1) # Exit with an error code

    # Assigning arguments to descriptive variables
    run_size       = sys.argv[1]
    fasta_path     = sys.argv[2]
    output_root    = sys.argv[3]
    project_name   = sys.argv[4]

    print(f"Starting job: {project_name} with batch size {run_size}...")
    
    split_jobs(run_size, fasta_path, output_root, project_name)
    
    print("Process complete.")

if __name__ == "__main__":
    main()
