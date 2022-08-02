#!/usr/bin/env python3
'''
2018-03-09
Ver.0.2 ---- Russell Xie
This version has updated the plotting function, including:
1) Reorder the genes from chr1 to chrY in the correct order.
2) Plot neiboring chromosomes in different color.
3) Plot a vertical line to indicate the target position.

Ver.0.3 ---- Russell Xie
Bug Fix: 
1) move several steps out from the paralelle functions to increase speed. 
2) fix the bug for inapprorpriate reading of chrX and chrY

Ver.0.4 ---- Russell Xie
Add option for identifying genes that are up-regulated.

Ver.0.5 ---- Russell Xie
Improve the multithreading performance
~30s per region using 48 cores

Ver.0.6 ---- Russell Xie
Apply test for all 1023 combinations of 10 sgRNAs in every enhancer region. 
Return a matrix containing all p values

Ver.0.7 ---- Russell Xie
Normalize matrix by CPM

Ver.0.8 ---- Russell Xie
Normalize the data to get rid of the batch effect

Ver.0.9 ---- Russell Xie
Update to match the annotation for hg38, only plot the genes from chr1-22,X

Ver.0.9.1 ---- Russell Xie
Add options whether to use cpm or meta_cell method for the normalization

Ver.0.9.2 ---- Russell Xie
Only consider the genes expressed in more than 1% of the cells

Ver.0.9.3 ---- Russell Xie
Output both the cdf and sf pval of the hypergeometric test

Note: requires the 'plot_annotation.txt' file to get the correct
gene position information.
'''
import os
import sys
import re
import collections
import argparse
import tables
import itertools
import matplotlib
import time
matplotlib.use('agg')

#change matplotlib font to type II for illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse

from multiprocessing import Pool
from collections import defaultdict
from scipy import sparse, io
from scipy.sparse import csr_matrix
from multiprocessing import Pool
from random import shuffle

np.random.seed(0)


#fine the indecies of cells which contain the input sgRNA
def find_sgrna_cells(sgRNA_df, transcriptome_df, sgRNA_dict):
    cell_bc = sgRNA_df.loc[:, (sgRNA_df.loc[sgRNA_dict].sum() > 0)].T.index
    cell_index = []
    for i in cell_bc:
        current_idx = np.argwhere(transcriptome_df.columns == i)
        if  current_idx.size > 0:
            cell_index.append(current_idx.item())
    return [x for x in set(cell_index)]

def find_non_zero_cells(sgrna_df, tanscriptome_df):
    cell_bc = sgrna_df.T.loc[sgrna_df.sum() != 0].index
    cell_index = []
    for i in cell_bc:
        current_idx = np.argwhere(transcriptome_df.columns == i)
        if  current_idx.size > 0:
            cell_index.append(current_idx.item())
    return [x for x in set(cell_index)]

def hypergeo_test(non_zero_array, sgrna_idx, i):
    #find indecies of cells in which expression of given gene is
    #equal or less than the median of this gene in the whole population
    median_cell_idx  = np.argwhere(non_zero_array <= np.median(non_zero_array))

    #find the same cells subset in the cells with a given sgRNA
    overlap_cell_idx = np.intersect1d(median_cell_idx, sgrna_idx)
    
    #calculate the median fold change
    other_idx = np.setxor1d(sgrna_idx, range(len(non_zero_array)))
    
    fc = (np.mean(non_zero_array[sgrna_idx]) + 0.01) / (np.mean(non_zero_array[other_idx]) + 0.01)
    
    #perform hypergeometric test, get the upper tail
    k = len(overlap_cell_idx)
    M = len(non_zero_array)
    n = len(median_cell_idx)
    N = len(sgrna_idx)
    try:
        pval_up = stats.hypergeom.logcdf(k, M, n, N).item()
    except:
        pval_up = float('nan')
        
    try:
        pval_down = stats.hypergeom.logsf(k, M, n, N).item()
    except:
        pval_down = float('nan')
    
    return pval_down, pval_up, fc

def perform_DE(sgrna_idx, input_array, idx, num_processes, pval_list_down, pval_list_up, fc_list):
    nonzero_pval_list_up = []
    nonzero_pval_list_down = []
    nonzero_fc_list = []
    with Pool(processes=num_processes) as p:
        for pval_down, pval_up, fc in p.starmap(hypergeo_test, zip(
                input_array,
                itertools.repeat(sgrna_idx),
                idx)
        ):
            nonzero_pval_list_down.append(pval_down)
            nonzero_pval_list_up.append(pval_up)
            nonzero_fc_list.append(fc)
    for i in idx:
        pval_list_up[i] = nonzero_pval_list_up.pop(0)
        pval_list_down[i] = nonzero_pval_list_down.pop(0)
        fc_list[i] = nonzero_fc_list.pop(0)
    return len(sgrna_idx), pval_list_up, pval_list_down, fc_list


def FDR(x):
    """
    Assumes a list or numpy array x which contains p-values for multiple tests
    Copied from p.adjust function from R  
    """
    o = [i[0] for i in sorted(enumerate(x), key=lambda v:v[1],reverse=True)]
    ro = [i[0] for i in sorted(enumerate(o), key=lambda v:v[1])]
    q = sum([1.0/i for i in range(1,len(x)+1)])
    l = [q*len(x)/i*x[j] for i,j in zip(reversed(range(1,len(x)+1)),o)]
    l = [l[k] if l[k] < 1.0 else 1.0 for k in ro]
    return np.asarray(l)

def create_combo(sgrna_list):
    combo_list = []
    for i in range(1, 11):
        my_array = list(itertools.combinations(enumerate(sgrna_list), i))
        for j in range(len(my_array)):
            idx = list(zip(*my_array[j]))[0]
            combo_list.append(idx)
    #return a list of tuple, which idicate the indexing of each combination
    return combo_list


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--hdf5', dest='input_hdf5', required=True,
                        type=str,
                        help='specify the hdf5 file (output from 10X pipeline).')
    parser.add_argument('-s', '--sgrna', dest='input_sgrna', required=True,
                        type=str,
                        help='specify the sgrna summary file.')
    parser.add_argument('-l', '--sgrna_list', dest='sgrna_list', required=True,
                        type=str,
                        help='specify the sgRNAs need to be tested.')
    parser.add_argument('-o', '--output_dir', dest='output_dir', required=False,
                        default = '.',
                        help='specify an output directory, default is current dir.')
    parser.add_argument('-t', '--threads', dest = 'threads', required=False,
                        type=int,
                        default=1,
                        help='set number of barcode comparison threads. \
                        The default is 1')
    parser.add_argument('-n', '--norm', dest = 'norm_method', required=False,
                        type=str,
                        default='cpm',
                        help='choose normalization methods: CPM only or \
                        normalize to the Meta-cell.')

    args = parser.parse_args()

    num_processing = args.threads
    
    #check the normalization method
    norm = args.norm_method
    if (norm != 'cpm') and (norm != 'metacell'):
        print("Incorrect normalization method. Has to be either 'cpm' or 'metacell'.", file = sys.stderr, flush=True)
        sys.exit(0)

    #read the 10X hdf5 file
    print("Reading HDF5 File...", file = sys.stderr, flush=True)
    sub_df_file = args.input_hdf5
    sub_df = pd.read_pickle(sub_df_file)
    
    output_dir = args.output_dir

    #read the plotting annotation
    v2_FILE = '/project/GCRB/Hon_lab/shared/former_members/s160875/03.analysis/Mosaic-seq/CROP-DE-analysis_10X-66K_no_downsampling-CPM.hg38/combine_10sgRNAs-volcano/generate_annotations/plot_annotation.txt'
    annot_df = pd.read_csv(v2_FILE,
                       header=None,
                       sep='\t',
                       names=['idx', 'gene_names', 'chromosome', 'pos', 'strand', 'color_idx', 'chr_idx'])
    
    #filter the genes expressed in more than 10% of the cells
    nonzero_idx = np.where(np.sum(sub_df > 0, axis = 1) > 1)[0]
    nonzero_idx = nonzero_idx[nonzero_idx != (len(sub_df.columns) - 1)]
    
    print("Finished nonzero_idx generation", file = sys.stderr, flush=True)
    
    idx = list(set(nonzero_idx) & set(annot_df.idx))
    del nonzero_idx
    del annot_df
        
    #normalize the matrix.
    [g, c] = sub_df.shape
    cpm_matrix = np.zeros((g, c))

    uniq_id = set()

    for x in sub_df.columns:
        uniq_id.add(x[-1])
        
    if (norm == 'cpm'):
        cpm_matrix = np.array(sub_df / sub_df.sum(axis = 0) * 1000000)
    elif (norm == 'metacell'):
        for lib in sorted(uniq_id):
            print("Normalizing Batch " + lib, file = sys.stderr, flush=True)
            index = [i for i,e in enumerate(sub_df.columns) if e.split('-')[1] == lib]
            one_cell_cpm = np.array((sub_df.iloc[:,index].sum(axis = 1) + 1)\
                            / np.sum(sub_df.iloc[:,index].values) * 1e6).flatten()
            for i in index:
                cell_cpm = sub_df.iloc[:,i] / np.sum(sub_df.iloc[:,i].values) * 1e6
                cpm_matrix[:,i] = cell_cpm / one_cell_cpm

    del uniq_id          
    print("Finished normalization", file = sys.stderr, flush=True)
    
    #create input ndarray
    input_array = cpm_matrix[idx]
    del cpm_matrix

    #load the sgRNA file
    print("Loading sgRNA df.", file = sys.stderr, flush=True)
    sgrna_df = args.input_sgrna
    sgrna_df_adj = pd.read_pickle(sgrna_df)

    #perform hypergeometric test for every single gene in the dataframe
    print("Loading sgRNA dictionary.", file = sys.stderr, flush=True)
    sgrna_dict  = {}
    sgrnas_file = args.sgrna_list
    with open(sgrnas_file) as f:
        for line in f:
            region_id, sgrna_string = line.strip().split("\t")
            sgrnas = sgrna_string.split(";")
            sgrna_dict.update({region_id : sgrnas})


    for k in sgrna_dict:
        print("Region in processing: " + k[0:], file = sys.stderr, flush=True)

        #idx index of cells containing the given sgRNA
        sgrna_idx = find_sgrna_cells(sgrna_df_adj, sub_df, sgrna_dict[k])

        #force the up-tail p-vals of all zero expressed genes to be zero. (actual p-val is 1)
        pval_list_down = np.zeros(len(sub_df.index))
        pval_list_up = np.zeros(len(sub_df.index))

        fc_list = np.ones(len(sub_df.index))
       
        print("DE analysis: " + k[0:], file = sys.stderr, flush=True)
        #perform the differential gene analysis by using Virtual FACS
        num_sgrna_cell, pval_list_up, pval_list_down, fc_list = perform_DE(
            sgrna_idx,
            input_array,
            idx,
            num_processing,
            pval_list_down,
            pval_list_up,
            fc_list
        )

        #save all the output
        io.savemat(
            '%s/%s-up_log-pval'%(output_dir, k[0:]),
            {'matrix':pval_list_up}
        )

        io.savemat(
            '%s/%s-down_log-pval'%(output_dir, k[0:]),
            {'matrix':pval_list_down}
        )
        io.savemat('%s/%s-%s-foldchange'%(output_dir, k[0:], num_sgrna_cell),
                   {'matrix':fc_list})

if __name__ == '__main__':
    main()
