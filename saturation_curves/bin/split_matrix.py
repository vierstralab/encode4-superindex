#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd

def main(samples_order, binary, masterlist, prefix_samples):
   
    index_samples = pd.read_table(samples_order, header=None, names=['ag_id'])
    subset = pd.read_table(prefix_samples, sep=",", header=None).T
    subset.columns = ['ag_id']

    # Find which columns the samples of interest are in
    matching_aggregations = index_samples[index_samples['ag_id'].isin(subset['ag_id'])]

    # Get the col numbers 
    col_nums = matching_aggregations.index

    aggregations = col_nums.to_list()

    # Load numpy matrix
    binary_np = np.load(binary)

    #Subset binary matrix columns
    prefix_binary = binary_np[:, aggregations].astype(int)

    #Subset binary matrix rows
    row_sums = np.sum(prefix_binary, axis=1)
    prefix_binary_subset = prefix_binary[row_sums >= 1]
    
    #Subset DHS Masterlist
    dhs_masterlist = pd.read_table(masterlist)
    dhs_masterlist_subset = dhs_masterlist[row_sums >= 1]

    #Save npy and masterlist subsets
    np.save('prefix_binary_matrix.npy', prefix_binary_subset)
    dhs_masterlist_subset.to_csv('dhs_masterlist_subset.bed', index=False, sep='\t', header=True) 

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file1', type=str, help="Path to the first file")
    parser.add_argument('file2', type=str, help="Path to the second file")
    parser.add_argument('file3', type=str, help="Path to the third file")
    parser.add_argument('file4', type=str, help="Path to the fourth file")
    
    args = parser.parse_args()
    
    main(args.file1, args.file2, args.file3, args.file4)
