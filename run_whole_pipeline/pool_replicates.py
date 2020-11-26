#!/usr/bin/env python
__description__ = \
"""
Pool the counts of peptide sequences between two experimental replicates.
"""
__author__ = "Michael J. Harms"
__date__ = "2020-11-25"
__usage__ = "./pool_replicates.py count_file_1 count_file_2 output_file"

import sys

def read_counts_file(count_file,count_dict=None):
    """
    Read a hops_counts output file and return as a dictionary keying peptide
    sequence to counts.

    count_file: string pointing to file
    count_dict: existing dictionary keying peptide to counts. If None, 
                start a new dictionary.
    """

    # Deal with intput count_dict
    if count_dict is None:
        count_dict = {}

    # Go through lines in count_file
    with open(count_file,'r') as f:
        for line in f:

            # skip comments
            if line.startswith("#"):
                continue

            # Grab counts
            col = line.strip().split()
            seq = col[0]
            counts = int(col[1])

            # Update dictionary
            try:
                count_dict[seq] += counts
            except KeyError:
                count_dict[seq] = counts


    return count_dict

def write_counts_file(output_file,count_dict):
    """
    Write counts to a hops_counts style file.
    
    output_file: string pointing to output file
    count_dict: dictionary to write out
    """
   
    f = open(output_file,'w') 
    for c in count_dict:
        f.write(f"{c} {count_dict[c]}\n")
    f.close()
        
             

def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]

    try:
        rep_1_file = argv[0]
        rep_2_file = argv[1]
        output_file = argv[2]
    except IndexError:
        err = f"incorrect arguments. usage:\n\n{__usage__}\n\n"
        raise ValueError(err)

    count_dict = read_counts_file(rep_1_file)
    count_dict = read_counts_file(rep_2_file,count_dict) 

    write_counts_file(output_file,count_dict)    

if __name__ == "__main__":
    main()
