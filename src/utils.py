import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import scipy.stats as ss


def do_stats(data):
    df2 = data.copy()
    dists = list(df2['distance'])
    glist = list(df2['gene'])
    dists2 = list(-df2['distance'])
    glist2 = list(df2['gene'])
    disttot = dists + dists2
    glists = glist + glist2
    dffdr = pd.DataFrame({'distance': disttot, 'gene': glists})
    mean = np.mean(dffdr['distance'])
    std = np.std(dffdr['distance'])
    dffdr['z'] = (dffdr['distance'] - mean) / std
    dffdr['p'] = norm.sf(dffdr['z']) * 2
    dffdr['fdr'] = multipletests(dffdr['p'], method='fdr_bh')[1]
    return dffdr


def fdrlistgen(data, thr=0.01, full=False):
    df = data.copy()
    dffdr = do_stats(df)
    if not full:
        dffdr = dffdr.drop(dffdr.index[dffdr['fdr'] > thr])
        dffdrlist = list(dffdr['gene'])
        fdrnum = len(dffdrlist)
        print('Number of genes passing FDR < {} is {}'.format(thr, fdrnum))
    if full:
        dffdr = dffdr.head(df.shape[0])
    return dffdr
