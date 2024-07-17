#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 15:59:00 2021

@author: maryamsamieinasab
"""

import re
from joblib import Parallel, delayed
from tqdm import tqdm
import numpy as np
import pandas as pd
from pysam import VariantFile
import argparse
from utils import *
from EPIMUTESTR import *
from Wavelet_Functions import *
from EAML_Functions import *
from Sigma_diff import *
from EAPathway_Functions import EA_Pathway_Wrapper
from EA_Pathways_VCFparser import vcf_parser
from vcf import *
import os


 
def parse_args():
    """
    Parses the Big pipeline arguments.
    """
    parser = argparse.ArgumentParser(description="Big Pipeline arguments")
    parser.add_argument('--VCF', nargs='?', default='./',help='Location of Cohort VCF')
    parser.add_argument('--ref', nargs='?', default='hg38', choices=('hg19', 'hg38'), help='genome reference file')
    parser.add_argument('--savepath', nargs='?', default='./',help='save path for output')
    parser.add_argument('--samples', nargs='?', default='./samples.txt',help='samples file path')
    parser.add_argument('--maxaf', type=float, default=0.01, help='maximum allele frequency cutoff')
    parser.add_argument('--minaf', type=float, default=0, help='minimum allele frequency cutoff')
    parser.add_argument('--transcript', default='canonical', choices=('all', 'max', 'mean', 'canonical'),help='how to parse EA scores from different transcripts')
    parser.add_argument('--Ann', default='VEP', choices=('VEP', 'ANNOVAR'),help='EA annotation method')
    parser.add_argument('--refPop', nargs='?', default='./refs/UKB_200K_WES.AC.AF.12102020.txt',help='text file containing reference population variants')
    parser.add_argument('--Groups', nargs='?', default='./refs/Reactomes2019_less100_NoSpecialChr.csv',help='biological groups of interest')
    parser.add_argument('--minAC',type=int , default=1, help='Min Allele Count Threshold for EA-Pathway Analysis')
    parser.add_argument('--maxAC',type=int , default=5, help='Max Allele Count Threshold for EA-Pathway Analysis')
    parser.add_argument('--cores', type=int, default=1, help='number of CPUs to use for multiprocessing')
    parser.add_argument('--pipeline', default='All', choices=('All', 'ML', 'Reactome', 'STRING','GoTerms', 'EAML', 'EPI', 'Wavelet', 'sigma'),help='which pipeline to be run')
    parser.add_argument('--JVMmemory', default='Xmx2g', help='memory argument for each Weka JVM')
    parser.add_argument('--chrX', type=int,default=1, help='1 if there is sex chromosome in the VCF, 0 if there is no sex chromosome in the VCF')
    parser.add_argument('--GeneLength', nargs='?', default='./refs/gene_length.csv',help='gene length file path')
    return parser.parse_args()


def main(args):
    if args.chrX==1:
        if args.Ann=='ANNOVAR':
            if args.ref=='hg19':
                ref = pd.read_csv('./refs/refGene-lite_hg19.May2013.txt', delimiter='\t', header=0, index_col='gene')
            elif args.ref=='hg38':
                ref = pd.read_csv('./refs/refGene-lite_hg38.June2017.txt', delimiter='\t', header=0, index_col='gene')
        elif args.Ann=='VEP':
            if args.ref=='hg19':
                ref = pd.read_csv('./refs/ENSEMBL-lite_GRCh37.v75.txt', delimiter='\t', header=0, index_col='gene')
            elif args.ref=='hg38':
                ref = pd.read_csv('./refs/ENSEMBL-lite_GRCh38.v94.txt', delimiter='\t', header=0, index_col='gene')
    elif (args.chrX==0) & (args.ref=='hg38') & (args.Ann=='VEP'):
        ref = pd.read_csv('./refs/ENSEMBL-lite_GRCh38.v94.noX.txt', delimiter='\t', header=0, index_col='gene')
    ref = ref[~ref.index.duplicated(keep='first')]    
    text_file = open(args.savepath+'arguments.txt', 'w')
    text_file.write('You are running '+str(args.pipeline)+' pipeline for:'+'\n'+str(args.minaf)+'<Allele Frequency<'+str(args.maxaf)+ '\n'+ 
                    'Genome Reference:'+str(args.ref))
    text_file.close()
    print('You are running', args.pipeline,'pipeline for:\n',args.minaf,'<Allele Frequency<',args.maxaf, '\n', 'Genome Reference:',args.ref)
    
    if args.pipeline=='All' or args.pipeline=='Reactome' or args.pipeline=='STRING' or args.pipeline=='GoTerms':
##code for parsing VCF into EA-Pathways input file
        sample_file = pd.read_csv(args.samples, index_col=0,header=None)
        cases = sample_file[sample_file.iloc[:,0]==1].index.tolist()
        conts = sample_file[sample_file.iloc[:,0]==0].index.tolist()
        os.makedirs(args.savepath+'Pathway_output', exist_ok = True)
        Pathway_output_path =  args.savepath+'Pathway_output/'
        vcf_parser(args.VCF, Pathway_output_path, args.refPop, args.minAC, args.maxAC, cases, conts, args.Ann)
        print('Prased VCF for EAPathway is completed')
   
        for ac in range(args.minAC, args.maxAC + 1):
            # ac+=1
            sample_input_df_Cases = pd.read_csv(args.savepath+'Pathway_output/'+'Input_files/' + 'Cases_PathwaysInput_AC' + str(ac) + '.csv', header=0)
            sample_input_df_Controls = pd.read_csv(args.savepath+'Pathway_output/'+'Input_files/' + 'Controls_PathwaysInput_AC' + str(ac) + '.csv', header=0)
            os.makedirs(Pathway_output_path + 'AC' + str(ac), exist_ok = True)
            output_dir = Pathway_output_path + 'AC' + str(ac)+'/'
            if args.pipeline=='Reactome' or args.pipeline=='All':
                Reactome_input_df = pd.read_csv('./refs/Reactome2023_Greater5Less100_03032023.csv', header=None)
                EA_Pathway_Wrapper(sample_input_df_Cases, Reactome_input_df, output_dir, 'Reactome','Cases', args.cores)
                EA_Pathway_Wrapper(sample_input_df_Controls, Reactome_input_df, output_dir, 'Reactome','Controls', args.cores)
            if args.pipeline=='STRING' or args.pipeline=='All':
                STRING_input_df = pd.read_csv('./refs/STRINGv11_Greater5Less100_02022021.csv', header=None) 
                EA_Pathway_Wrapper(sample_input_df_Cases, STRING_input_df, output_dir, 'STRING','Cases', args.cores)
                EA_Pathway_Wrapper(sample_input_df_Controls, STRING_input_df, output_dir, 'STRING','Controls', args.cores)
            if args.pipeline=='GoTerms' or args.pipeline=='All':
                Goterms_input_df = pd.read_csv('./refs/GOterms_Greater5Less100_07202023.csv', header=None) 
                EA_Pathway_Wrapper(sample_input_df_Cases, Goterms_input_df, output_dir, 'GoTerms','Cases', args.cores)
                EA_Pathway_Wrapper(sample_input_df_Controls, Goterms_input_df, output_dir, 'GoTerms','Controls', args.cores)
## Sigma Diff Analysis    
    if args.pipeline=='All' or args.pipeline=='sigma':   
        samples = pd.read_csv(args.samples, index_col=0,header=None)
        cases = samples[samples.iloc[:,0]==1].index.astype(str).tolist()
        conts = samples[samples.iloc[:,0]==0].index.astype(str).tolist()
        total_samples = samples.index.astype(str).tolist()
        print('\n Sigma Diff analysis started')
        os.makedirs(args.savepath+'Sigma_output', exist_ok = True)
        Sigma_output_path =  args.savepath+'Sigma_output/' 
      #### Generate sumEA matrix for Sigma Diff Analysis
        if args.Ann=='ANNOVAR':
            matrix = Parallel(n_jobs=args.cores)(delayed(parse_ANNOVAR_Sigma)(args.VCF, gene, ref.loc[gene], total_samples, min_af=args.minaf, max_af=args.maxaf,af_field='AF',EA_parser='canonical') for gene in tqdm(ref.index.unique()))
        if args.Ann=='VEP':
            matrix = Parallel(n_jobs=args.cores)(delayed(parse_VEP_Sigma)(args.VCF, gene, ref.loc[gene], total_samples, max_af=args.maxaf, min_af=args.minaf,af_field='AF') for gene in tqdm(ref.index.unique()))
        design_matrix = pd.concat(matrix, axis=1)
        ## reading gene length file
        gene_length = pd.read_csv(args.GeneLength, index_col=0)
        genes = set(design_matrix.columns.tolist()).intersection(set(gene_length.index.tolist()))
        gene_length = gene_length.loc[genes]
        
        sigma_matrix, distance_matrix  = sigma(design_matrix, gene_length, genes,samples)
        distance_matrix.to_csv(Sigma_output_path+'distance_matrix.tsv', sep='\t', header=True, index=True) 
        sigma_matrix.to_csv(Sigma_output_path+'sigma.tsv', sep='\t', header=True, index=True) 
        print('\n Sigma Diff analysis completed')    
    
    if args.pipeline=='All' or args.pipeline=='ML' or args.pipeline=='EAML':
## EA-ML Analysis 
        print('\n EAML analysis started')
        os.makedirs(args.savepath+'tmp', exist_ok = True)
        os.makedirs(args.savepath+'EAML_output', exist_ok = True)
        EAML_output_path =  args.savepath+'EAML_output/'
        gene_results = Parallel(n_jobs=args.cores)(delayed(eval_gene)(gene,reference=args.ref, data_fn=args.VCF,targets_fn=args.samples, expdir=args.savepath,
                                                                      min_af=args.minaf, max_af=args.maxaf, af_field='AF', EA_parser=args.transcript, EA_Ann=args.Ann, seed=111,cv=10,weka_path='./weka-3-8-5', memory=args.JVMmemory) for gene in tqdm(ref.index.unique()))
        raw_results = gene_results
        full_results,nonzero_results = report_results(raw_results,EAML_output_path)
        print('\n EAML analysis completed')
    
    if args.pipeline=='All' or args.pipeline=='ML' or args.pipeline=='Wavelet' or args.pipeline=='EPI':   
#### Generate pEA matrix for EA-Wavelet and EPIMUTESTR Analysis
        sample_file = pd.read_csv(args.samples, index_col=0,header=None)
        cases = sample_file[sample_file.iloc[:,0]==1].index.astype(str).tolist()
        conts = sample_file[sample_file.iloc[:,0]==0].index.astype(str).tolist()
        samples = cases + conts
        if args.Ann=='ANNOVAR':
            gene_dfs = Parallel(n_jobs=args.cores)(delayed(parse_ANNOVAR)(vcf_fn=args.VCF, gene=gene,gene_ref=ref.loc[gene],samples=samples, min_af=args.minaf, max_af=args.maxaf, af_field='AF', EA_parser=args.transcript) for gene in tqdm(ref.index.unique()))
            design_matrix = pd.concat(gene_dfs, axis=1) 
        elif args.Ann=='VEP':
            gene_dfs = Parallel(n_jobs=args.cores)(delayed(parse_VEP)(vcf_fn=args.VCF, gene=gene,gene_ref=ref.loc[gene],samples=samples, min_af=args.minaf, max_af=args.maxaf, af_field='AF', EA_parser=args.transcript) for gene in tqdm(ref.index.unique()))
            design_matrix = pd.concat(gene_dfs, axis=1) 
        design_matrix.to_csv(args.savepath+'input_matrix.csv', header=True, index=True)    
## EAWavelet Analysis   
        data, genelist, caselist, contlist, shift, case, cont = init(args.savepath+'input_matrix.csv', args.samples)
        if args.pipeline!='EPI':           
            print('\n EA-Wavelet Analysis Started') 
            case, cont, genelistnet = create_graphs(data, case, cont, caselist, contlist, genelist)
            TotChi, shift, case_node_map, cont_node_map = run_wavelet_decomp(case, cont, genelistnet, shift)
            distances = PCA_distance(TotChi, case_node_map, cont_node_map, shift)
            dffdr = fdrlistgen(data=distances, thr=1)
            os.makedirs(args.savepath+'EAWavelet_output', exist_ok = True)
            dffdr.to_csv(args.savepath+'EAWavelet_output/'+'wavelet_output.csv', header=True, index=False)
            print('\n EA Wavelet Analysis Complete')
        if args.pipeline!='Wavelet': 
## EPIMUTESTR Analysis 
            Data = data.reindex(sample_file.index.astype(str).tolist())
            Data['class'] = sample_file.iloc[:,0].values
            print('\n EPIMUTESTR Analysis Started')
            geneList = EPIMUTESTR(X=Data, top_features=200, n_cores=args.cores)
            print('\n EPIMUTESTR Analysis Complete')
            os.makedirs(args.savepath+'EPIMUTESTR_output', exist_ok = True)
            EPI_output_path =  args.savepath+'EPIMUTESTR_output/'
            fh = open(EPI_output_path + 'EPI_output.tsv', 'w')
            for key, value in geneList.items():
                fh.write(str(key) + '\t' + str(value[0]) + '\t' + str(value[1]) + '\n')
            fh.close()

if __name__ == "__main__":
    args = parse_args()
    main(args)


