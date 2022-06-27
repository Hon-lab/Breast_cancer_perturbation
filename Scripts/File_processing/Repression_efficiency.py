#!/usr/bin/env python3
#This script is used for filtering the CellPlex multiplex, no CellPlex label cells and sgRNA outlier cells. 
#Yihan Wang, Jan 27, 2021

import os
import re
import sys
import collections
import argparse
import tables
import itertools
import matplotlib
import numba

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse

from multiprocessing import Pool
from collections import defaultdict
from scipy import sparse, io
from scipy.sparse import csr_matrix

def find_sgrna_cells(sgRNA_df, transcriptome_df, sgRNA_dict):
    cell_bc = sgRNA_df.loc[:, (sgRNA_df.loc[sgRNA_dict].sum() > 0)].T.index
    cell_index = []
    for i in cell_bc:
        current_idx = np.argwhere(transcriptome_df.columns == i)
        if  current_idx.size > 0:
            cell_index.append(current_idx.item())
    return [x for x in set(cell_index)]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--transcriptome_df', dest='transcriptome_df', required=True,
                        type=str,
                        help='specify the filtered sub_df pickle file')

    parser.add_argument('-s', '--sgrna', dest='input_sgrna', required=True,
                        type=str,
                        help='specify the filtered sgrna pkl file.')
    
    parser.add_argument('-d', '-dict', dest='dict', required=True,
                        type=str, 
                        help='specify the sgRNA dictionary txt file.')

    parser.add_argument('-r', '--region', dest='region', required=True,
                        type=str,
                        help='specify the quired regions and their target genes to calculate repression efficiency.')

    parser.add_argument('-o', '--output_file', dest='output_file', required=True,
                        help='specify output file.')
        

    args = parser.parse_args()
    TRANSCRIPTOME_DF = args.transcriptome_df
    SGRNA = args.input_sgrna
    SGRNA_DICT = args.dict
    TARGET_FILE = args.region
    OUTPUT_FILE = args.output_file

    print('Loading data.', file=sys.stderr, flush=True)
    sub_df = pd.read_pickle(TRANSCRIPTOME_DF)
    sgrna_df_bool = pd.read_pickle(SGRNA) > 0
    
    #load sgRNA dictionary
    sgrna_dict  = {}
    with open(SGRNA_DICT) as f:
        for line in f:
            region_id, sgrna_string = line.strip().split('\t')
            sgrnas = sgrna_string.split(';')
            sgrna_dict.update({region_id : [i.upper() for i in sgrnas]})

    #load regions and their target genes
    region_dict = {}
    with open(TARGET_FILE) as ft:
        for line in ft:
            region, gene = line.strip().split('\t')
            region_dict.update({region: gene})
    
    #load annotation df 
    v2_FILE = '/project/GCRB/Hon_lab/shared/former_members/s160875/03.analysis/Mosaic-seq/CROP-DE-analysis_10X-66K_no_downsampling-CPM.hg38/combine_10sgRNAs-volcano/generate_annotations/plot_annotation.txt'
    annot_df = pd.read_csv(v2_FILE,
                   header=None,
                   sep='\t',
                   names=['idx', 'gene_names', 'chromosome', 'pos', 'strand', 'color_idx', 'chr_idx'])

    print('Finished loading data.', file=sys.stderr, flush=True)
    
    #normalization
    cpm_matrix = np.array(sub_df / sub_df.sum(axis = 0) * 1000000)
    print('Finished normalization.', file=sys.stderr, flush=True)

    #calculate repression efficiency
    fc_list = []
    for i in region_dict.keys():
        target_gene = region_dict[i]
        print('Calculating gene: ' + str(target_gene) + ' repression in region: ' + str(i), file=sys.stderr, flush=True)
        target_gene_idx = annot_df[annot_df['gene_names'] == target_gene].idx.values[0]

        sgrna_idx = find_sgrna_cells(sgrna_df_bool, sub_df, sgrna_dict[i])
        non_zero_array = cpm_matrix[target_gene_idx]
        other_idx = np.setxor1d(sgrna_idx, range(len(non_zero_array)))
        sgrna_cells_expression = cpm_matrix[target_gene_idx][sgrna_idx]
        other_cells_expression = cpm_matrix[target_gene_idx][other_idx]
        fc = (np.mean(sgrna_cells_expression) + 0.01) / (np.mean(other_cells_expression) + 0.01)
        fc_list.append(fc)
    
    #write to output file
    output_file=open(OUTPUT_FILE, 'w')
    for i in region_dict.keys():
        target_gene = region_dict[i]
        fc = fc_list[list(region_dict.keys()).index(i)]
        output_file.write(str(i) + ': ' + str(target_gene) + ': ' + str(fc) + '\n')
    
    output_file.close()


if __name__ == '__main__':
    main()
