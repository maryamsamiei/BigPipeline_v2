"""
Author: Saeid Parvandeh November 2019
"""

import numpy as np
import pandas as pd
import math

from skrebate.turf import TuRF
from skrebate.relieff import ReliefF
import statsmodels.stats.multitest as fdr
import scipy


def EPIMUTESTR(X=None, top_features=200, nn=None, discrete_threshold=10, verbose=False, n_cores=1, estimator='relief', pct=0.5):
    X = X.loc[:, (X != 0).any(axis=0)]
    features, labels = X.drop('class', axis=1).values, X['class'].values
    features = np.nan_to_num(features)
    headers = list(X.drop("class", axis=1))
    if nn == None:
        nn = math.floor(0.154 * (X.shape[1] - 1))
    if (estimator=='TuRF'):
        # Total Unduplicated Reach and Frequency (TURF)
        fs = TuRF(core_algorithm="ReliefF", n_features_to_select=top_features, n_neighbors=nn, pct=pct, verbose=verbose, n_jobs=n_cores)
        fs.fit(features, labels, headers)
    elif (estimator == 'relief'):
        # ReliefF stand alone
        fs = ReliefF(n_features_to_select=top_features, n_neighbors=nn, discrete_threshold=discrete_threshold, verbose=verbose, n_jobs=n_cores)
        fs.fit(features, labels)

    scoreDict = dict(zip(X.drop('class', axis=1).columns, fs.feature_importances_))
    scoreDict_sorted = {i[1]: i[0] for i in sorted(zip(scoreDict.values(), scoreDict.keys()), reverse=True)}
    scores_list = list(scoreDict_sorted.values())
    pos_scores_list = [n for n in scores_list if n > 0]
    # calculate the P value and adjusted P value
    gene_scores = np.sqrt(pos_scores_list)
    gene_scores_mean = np.mean(gene_scores)
    gene_scores_sd = np.std(gene_scores)
    pvals = []
    for score in gene_scores:
        pvals.append(scipy.stats.norm(gene_scores_mean, gene_scores_sd).sf(score))
    # Benjamini/Hachberg FDR correction
    qvals = fdr.multipletests(np.asarray(pvals), method='fdr_bh', is_sorted=True)

    geneList = dict(zip(scoreDict_sorted.keys(), zip(pvals, qvals[1])))

    return geneList
