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

def process_feature_df(mtx, barcodes, feature):
    gene_name = []
    for i in feature:
        gene_name.append(i[1])
    transcriptome_df = pd.DataFrame.sparse.from_spmatrix(data = mtx.tocsr(),
                                                         columns = barcodes, index = gene_name)
    return transcriptome_df

def filter_umi (df, copy=False):
    df = df.copy() if copy else df
    feature_cutoff = [[turn_point(i, df)] for i in list(df.index)]
    
    for i in range(0, len(feature_cutoff)):
        ZERO_HTO = df.iloc[i, :].loc[df.iloc[i, :] <= feature_cutoff[i][0]].index
        df.at[df.index[i], ZERO_HTO] = 0
    return df, feature_cutoff

def turn_point(sgRNA_name, df):
    sgRNA_count  = df.T.filter(items=[sgRNA_name]).sum(axis=1).sort_values(ascending=False)
    sgRNA_cumsum = sgRNA_count.cumsum()

    #get the total cell number of this sgRNA
    cell_num = np.argwhere(sgRNA_count > 0).size

    #calculate the turning point by using the max derivative
    turning_point = sgRNA_cumsum.loc[((sgRNA_cumsum.diff()) / sgRNA_count.sum() > (1/cell_num))].shape
    
    return(sgRNA_count.iloc[turning_point])

def get_sgrna_per_cell(df, return_mean=True, return_median=True):
    CorrSgrnaPerCell = np.sum(df > 0, 0)
    sgrna_mean = np.mean(CorrSgrnaPerCell)
    sgrna_median = np.median(CorrSgrnaPerCell)
    print("Average sgRNA Number per Cell is: " + str(sgrna_mean), file=sys.stderr, flush=True)
    
    if return_mean & return_median:
        return CorrSgrnaPerCell, sgrna_mean, sgrna_median
    elif return_mean:
        return CorrSgrnaPerCell, sgrna_mean
    elif return_median:
        return CorrSgrnaPerCell, sgrna_median
    else:
        return CorrSgrnaPerCell

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--feature_bc', dest='feature_bc', required=True,
                        type=str,
                        help='specify the output folder from cellranger pipeline (filtered_feature_bc_matrix folder)')

    parser.add_argument('-s', '--sgrna', dest='input_sgrna', required=True,
                        type=str,
                        help='specify the sgrna pkl file. If multiple libs, combine together')

    parser.add_argument('-o', '--output_directory', dest='output_directory', required=True,
                        help='specify output directory.')
        

    args = parser.parse_args()
    WORK_DIR = args.feature_bc
    SGRNA = args.input_sgrna
    OUTPUT_DIR = args.output_directory

    print('Loading transcriptome data.', file=sys.stderr, flush=True)
    
    #calculate HTO and filter
    Count_list = []
    MTX = io.mmread(WORK_DIR + 'matrix.mtx.gz')
    barcodes = np.loadtxt(WORK_DIR + 'barcodes.tsv.gz', dtype=str)
    features = np.loadtxt(WORK_DIR + 'features.tsv.gz', dtype=str)
    HTO_df = process_feature_df(MTX, barcodes, features).iloc[-12:,:]
    HTO_df_dense = HTO_df.sparse.to_dense()
    HTO_df_adj, _ = filter_umi(HTO_df_dense, copy=True)
    HTO_count = np.sum(HTO_df_adj.iloc[0:12, :] > 0, axis=0).values
    Count_list.append(HTO_count)

    HTO_df_adj_bool = HTO_df_adj > 0
    singlet_HTO_df = HTO_df_adj_bool.T[(HTO_df_adj_bool.sum(axis=0).values == 1)]

    singlet_rate = collections.Counter(Count_list[0])[1] / len(barcodes)
    print('singlet rate: ' + str(singlet_rate), file=sys.stderr, flush=True)

    sgRNA_df = pd.read_pickle(SGRNA)
    print('The size of sgRNA df before filter: ' + str(sgRNA_df.shape), file=sys.stderr, flush=True)
    sgRNA_df_bool = sgRNA_df > 0

    #split the cells into different antibodies
    CellPlex = ['CMO301', 'CMO302', 'CMO303', 'CMO304', 'CMO305', 'CMO306', 'CMO307', 'CMO308', 'CMO309', 'CMO310', 'CMO311', 'CMO312']

    singlet_cell_ID = []
    for i in CellPlex:
        cell_ID = set(list(singlet_HTO_df.index[singlet_HTO_df[i] == True])).intersection(set(sgRNA_df_bool.columns))
        singlet_cell_ID.append(cell_ID)

    sgRNA_num = []
    for i in np.arange(len(singlet_cell_ID)):
        sgRNA_num.append(list(np.sum(sgRNA_df_bool[singlet_cell_ID[i]], axis=0).values))

    for i in np.arange(len(sgRNA_num)):
        print('CellPlex ' + str(i) + ' sgRNA median: ' + str(np.median(sgRNA_num[i])), file=sys.stderr, flush=True)
    
    del sgRNA_df_bool
    trans_df = process_feature_df(MTX, barcodes, features)
    sub_df = trans_df[HTO_df_adj.columns[(HTO_df_adj > 0).sum(axis=0).values == 1].values].iloc[:-12,:]
    print('Cells before filter HTO: ' + str(len(trans_df.columns)), file=sys.stderr, flush=True)
    del HTO_df
    del HTO_count
    del HTO_df_adj
    del HTO_df_adj_bool
    del HTO_df_dense



    new_columns = set(sgRNA_df.columns.values).intersection(set(sub_df.columns.values))
    new_sub_df = sub_df[new_columns]
    del sub_df
    del trans_df
    print('Cells after filter HTO: ' + str(len(new_columns)), file=sys.stderr, flush=True)
    new_sgRNA_df = sgRNA_df[new_columns]
    del sgRNA_df
    #filter sgRNA outliers 
    CorrSgrnaPerCell, sgrna_mean , sgrna_median = get_sgrna_per_cell(new_sgRNA_df, return_mean=True, return_median=True)
    print('sgRNA mean before filter sgRNA outlier: ' + str(sgrna_mean), file=sys.stderr, flush=True)
    print('sgRNA median before filter sgRNA outlier: ' + str(sgrna_median), file=sys.stderr, flush=True)
    Q1 = CorrSgrnaPerCell.quantile(0.25)
    Q3 = CorrSgrnaPerCell.quantile(0.75)
    IQR = Q3 - Q1

    boston_df_out = new_sgRNA_df.loc[:,((CorrSgrnaPerCell < (Q1 - 1.5 * IQR)) |(CorrSgrnaPerCell > (Q3 + 1.5 * IQR))).values == False]
    del new_sgRNA_df
    CorrSgrnaPerCell_out, sgrna_mean_out , sgrna_median_out = get_sgrna_per_cell(boston_df_out, return_mean=True, return_median=True)
    print('sgRNA mean after filter sgRNA outlier: ' + str(sgrna_mean_out), file=sys.stderr, flush=True)
    print('sgRNA median after filter sgRNA outlier: ' + str(sgrna_median_out), file=sys.stderr, flush=True)
    print('Cells after filter sgRNA outliers: ' + str(len(boston_df_out.columns)), file=sys.stderr, flush=True)


    #write output file
    stats_output = open(OUTPUT_DIR + 'stats.txt', 'w')
    stats_output.write('singlet rate: ' + str(singlet_rate))
    stats_output.write('\n')
    stats_output.write('Cells after filter HTO: ' + str(len(new_columns)))
    stats_output.write('\n')
    stats_output.write('sgRNA mean before filter sgRNA outlier: ' + str(sgrna_mean))
    stats_output.write('\n')
    stats_output.write('sgRNA median before filter sgRNA outlier: ' + str(sgrna_median))
    stats_output.write('\n')
    stats_output.write('sgRNA mean after filter sgRNA outlier: ' + str(sgrna_mean_out))
    stats_output.write('\n')
    stats_output.write('sgRNA median after filter sgRNA outlier: ' + str(sgrna_median_out))
    stats_output.write('\n')
    stats_output.write('Cells after filter duplicates and sgRNA outliers: ' + str(len(boston_df_out.columns)))
    stats_output.close()


    new_sub_df.loc[:,boston_df_out.columns].to_pickle(OUTPUT_DIR + 'Singlet_sub_df.pkl')
    boston_df_out.to_pickle(OUTPUT_DIR + 'Singlet_sgRNA_df.pkl')
if __name__ == '__main__':
    main()
