#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 11:42:55 2021

@author: Dilon Shapiro
"""

#!/usr/bin/env python
"""Main script for EA-ML pipeline."""
import datetime
import shutil
import time
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from joblib import delayed, Parallel
from pkg_resources import resource_filename
from scipy import stats
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
from vcf import parse_gene_eaml
from visualize import mcc_hist, mcc_scatter, manhattan_plot
from weka_wrapper import eval_gene1


def compute_gene_dmatrix(gene,reference, data_fn, targets_fn, write_data, min_af, max_af, af_field, EA_parser, expdir):
    """
    Computes the full design matrix from an input VCF

    Args:
        gene (str): HGSC gene symbol

    Returns:
        DataFrame: EA design matrix for gene-of-interest
    """
    targets = pd.read_csv(targets_fn, header=None, dtype={0: str, 1: int}).set_index(0).squeeze().sort_index()
    if reference=='hg19':
        ref = pd.read_csv('./refs/hg19-refGene.protein-coding.txt', delimiter='\t', header=0, index_col='name2')
    elif reference=='hg38':
        ref = pd.read_csv('./refs/hg38-refGene.protein-coding.txt', delimiter='\t', header=0, index_col='name2')
    gene_reference = ref.loc[gene]
    dmatrix = parse_gene_eaml(data_fn, gene, gene_reference, list(targets.index), min_af,max_af, af_field, EA_parser)
    # if write_data==1:
    #     written = False
    #     while written is False:
    #         try:
    #             dmatrix.to_hdf(expdir+'dmatrices.h5', key=gene, complevel=5, complib='zlib', format='fixed')
    #             written = True
    #         except HDF5ExtError:
    #             pass
    return dmatrix



def eval_gene(gene,reference,data_fn,targets_fn,expdir,write_data, min_af, max_af, af_field, EA_parser,seed=111,cv=10,weka_path='./weka-3-8-5/'):
    """
    Parses input data for a given gene and evaluates it using Weka

    Args:
        gene (str): HGSC gene symbol

    Returns:
        dict(float): Mapping of classifier to MCC from cross validation
    """
    class_params = {
        'PART': '-M 5 -C 0.25 -Q 1',
        'JRip': '-F 3 -N 2.0 -O 2 -S 1 -P',
        'RandomForest': '-I 10 -K 0 -S 1',
        'J48': '-C 0.25 -M 5',
        'NaiveBayes': '',
        'Logistic': '-R 1.0E-8 -M -1',
        'IBk': '-K 3 -W 0 -A \".LinearNNSearch -A \\\".EuclideanDistance -R first-last\\\"\"',
        'AdaBoostM1': '-P 100 -S 1 -I 10 -W .DecisionStump',
        'MultilayerPerceptron': '-L 0.3 -M 0.2 -N 500 -V 0 -S 0 -E 20 -H a'
    }
    targets = pd.read_csv(targets_fn, header=None, dtype={0: str, 1: int}).set_index(0).squeeze().sort_index()
    gene_dmatrix = compute_gene_dmatrix(gene,reference, data_fn,targets_fn, write_data, min_af, max_af, af_field, EA_parser,expdir)
    mcc_results = eval_gene1(gene, gene_dmatrix, targets, class_params, 111, 10, Path(expdir), weka_path='./weka-3-8-5/')
    (Path(expdir) / f'tmp/{gene}.arff').unlink()  # clear intermediate ARFF file after gene scoring completes
    return gene, mcc_results




def report_results(raw_results,expdir):
    """Summarize and rank gene scores"""
    mcc_df_dict = defaultdict(list)
    for gene, mcc_results in raw_results:
        mcc_df_dict['gene'].append(gene)
        for clf, mcc in mcc_results.items():
            mcc_df_dict[clf].append(mcc)
    mcc_df = pd.DataFrame(mcc_df_dict).set_index('gene')
    clfs = mcc_df.columns
    mcc_df['mean'] = mcc_df.mean(axis=1)
    mcc_df['std'] = mcc_df[clfs].std(axis=1)
    mcc_df.sort_values('mean', ascending=False, inplace=True)
    mcc_df.to_csv(expdir + 'classifier-MCC-summary.csv')
    full_results = mcc_df
    mcc_df = full_results[['mean', 'std']]
    nonzero = mcc_df.loc[mcc_df[f'mean'] != 0].copy()
    nonzero.rename(columns={'mean': 'MCC'}, inplace=True)
    nonzero['logMCC'] = np.log(nonzero.MCC + 1 - np.min(nonzero.MCC))
    nonzero['zscore'] = (nonzero.logMCC - np.mean(nonzero.logMCC)) / np.std(nonzero.logMCC)
    nonzero['pvalue'] = stats.norm.sf(abs(nonzero.zscore)) * 2
    nonzero['qvalue'] = multipletests(nonzero.pvalue, method='fdr_bh')[1]
    stats_df = nonzero
    stats_df.to_csv(expdir + 'meanMCC-results.nonzero-stats.rankings')
    nonzero_results = stats_df
    return full_results,nonzero_results


def visualize(full_results, nonzero_results, expdir, reference):
    """
    Generate summary figures of EA-ML results, including a Manhattan plot of p-values and scatterplots and
    histograms of MCC scores
    """
    mcc_scatter(full_results, column='mean', dpi=150).savefig(expdir + f'meanMCC-scatter.png')
    mcc_hist(full_results, column='mean', dpi=150).savefig(expdir + f'meanMCC-hist.png')
    mcc_scatter(nonzero_results, column='MCC', dpi=150).savefig(expdir + 'meanMCC-scatter.nonzero.png')
    mcc_hist(nonzero_results, column='MCC', dpi=150).savefig(expdir + 'meanMCC-hist.nonzero.png')
    manhattan_plot(nonzero_results, reference, dpi=150).savefig(expdir + 'MCC-manhattan.svg')




