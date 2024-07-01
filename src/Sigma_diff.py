#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 01:41:35 2024

@author: maryamsamieinasab
"""
import numpy as np
import pandas as pd
import random
from scipy import stats
from statsmodels.stats.multitest import multipletests

def sigma(design_matrix, gene_length, genes,samples):
    cases = samples[samples.iloc[:,0]==1].index.astype(str).tolist()
    conts = samples[samples.iloc[:,0]==0].index.astype(str).tolist()
    total_samples = samples.index.astype(str).tolist()
    
    design_matrix_case = design_matrix.loc[cases,genes]
    design_matrix_control = design_matrix.loc[conts,genes]
    
    SumEA_genes_case = np.sum(design_matrix_case,axis=0)
    SumEA_genes_control = np.sum(design_matrix_control,axis=0)
    
    expected_energy_case = np.sum(SumEA_genes_case)/np.sum(gene_length['gene_length'])
    sigma_matrix_case = pd.DataFrame(np.zeros((len(genes), 3)), index=genes, columns=['sumEA', 'gene_length', 'sigma'])
    sigma_matrix_case['sumEA']= SumEA_genes_case
    sigma_matrix_case['gene_length'] = gene_length
    sigma_matrix_case['sigma'] = expected_energy_case/(sigma_matrix_case['sumEA']/sigma_matrix_case['gene_length'])
    
    ## mu calculation
    expected_energy_control = np.sum(SumEA_genes_control)/np.sum(gene_length['gene_length'])
    sigma_matrix_control = pd.DataFrame(np.zeros((len(genes), 3)), index=genes, columns=['sumEA', 'gene_length', 'sigma'])
    sigma_matrix_control['sumEA']= SumEA_genes_control
    sigma_matrix_control['gene_length'] = gene_length
    sigma_matrix_control['sigma'] = expected_energy_control/(sigma_matrix_control['sumEA']/sigma_matrix_control['gene_length'])
    sigma_matrix = pd.DataFrame(np.zeros((len(genes), 2)), index=genes, columns=['sigma_case', 'sigma_control'])
    sigma_matrix['sigma_case'] = sigma_matrix_case['sigma'].copy()
    sigma_matrix['sigma_control'] = sigma_matrix_control['sigma'].copy()
    
    distance_matrix = pd.DataFrame(np.zeros((len(genes), 1)), index=genes, columns=['distance'])
    distance_matrix['distance'] = sigma_matrix_control['sigma']-sigma_matrix_case['sigma']
    
    for i in range(1000):
        cases1 = random.sample(total_samples, len(cases))
        controls1 = list(set(total_samples) - set(cases1))
        design_matrix_case = design_matrix.loc[cases1,genes]
        design_matrix_control = design_matrix.loc[controls1,genes]

        SumEA_genes_case = np.sum(design_matrix_case,axis=0)
        SumEA_genes_control = np.sum(design_matrix_control,axis=0)
    
        expected_energy_case = np.sum(SumEA_genes_case)/np.sum(gene_length['gene_length'])
        sigma_matrix_case['sumEA']= SumEA_genes_case
        sigma_matrix_case[str(i)] = expected_energy_case/(sigma_matrix_case['sumEA']/sigma_matrix_case['gene_length'])
    
        expected_energy_control = np.sum(SumEA_genes_control)/np.sum(gene_length['gene_length'])
        sigma_matrix_control['sumEA']= SumEA_genes_control
        sigma_matrix_control[str(i)] = expected_energy_control/(sigma_matrix_control['sumEA']/sigma_matrix_control['gene_length'])
    
        distance_matrix[str(i)] = sigma_matrix_control[str(i)]-sigma_matrix_case[str(i)]
    
    distance_matrix1 = distance_matrix.drop(columns='distance')
    mean = np.mean(distance_matrix1, axis=1)
    std = np.std(distance_matrix1, axis=1)
    distance_matrix['zscore'] = (distance_matrix['distance']-mean)/std 
    distance_matrix['pvalue'] = stats.norm.sf(abs(distance_matrix.zscore)) * 2
    distance_matrix['qvalue'] = multipletests(distance_matrix.pvalue, method='fdr_bh')[1]
    
    return sigma_matrix, distance_matrix       

