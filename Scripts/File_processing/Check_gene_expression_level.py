#!/usr/bin/env python3
#This script is used for checking the gene expression level(cpm)
#Yihan Wang, Feb 17, 2022

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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--transcriptome_df', dest='transcriptome_df', required=True,
                        type=str,
                        help='specify the filtered sub_df pickle file')

    parser.add_argument('-g', '--gene', dest='gene', required=True,
                        type=str,
                        help='specify the quired genes.')

    parser.add_argument('-o', '--output_file', dest='output_file', required=True,
                        help='specify output file.')
        

    args = parser.parse_args()
    TRANSCRIPTOME_DF = args.transcriptome_df
    GENE_FILE = args.gene
    OUTPUT_FILE = args.output_file

    print('Loading data.', file=sys.stderr, flush=True)
    sub_df = pd.read_pickle(TRANSCRIPTOME_DF)

    #load quired genes
    GENE_LIST = []
    with open(GENE_FILE) as ft:
        for line in ft:
            gene = line.strip()
            GENE_LIST.append(gene)
    
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

    #calculate average cpm and write to output file
    output_file=open(OUTPUT_FILE, 'w')
    
    for gene in GENE_LIST:
        gene_idx = annot_df[annot_df['gene_names'] == gene]['idx'].values[0]
        ave_cpm = np.mean(cpm_matrix[gene_idx])
        perc_cell = np.sum(cpm_matrix[gene_idx] > 0) /len(cpm_matrix[gene_idx])
        output_file.write(str(gene) + ': ' + str(ave_cpm) + ': ' + str(perc_cell) + '\n')
    
    output_file.close()


if __name__ == '__main__':
    main()
