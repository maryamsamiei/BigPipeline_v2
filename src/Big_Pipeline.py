#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 15:06:08 2021

@author: maryamsamieinasab

Main script to run BigPipeline

"""

import re
from joblib import Parallel, delayed
from tqdm import tqdm
import numpy as np
import pandas as pd
from pysam import VariantFile
import argparse
# from graphwave_py3 import *
from utils import *
from EPIMUTESTR import *
from Wavelet_Functions import *
from EAML_Functions import *
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
    parser.add_argument('--refPop', nargs='?', default='./refs/UKB_200K_WES.AC.AF.12102020.txt',help='text file containing reference population variants')
    parser.add_argument('--Groups', nargs='?', default='./refs/Reactomes2019_less100_NoSpecialChr.csv',help='biological groups of interest')
    parser.add_argument('--minAC',type=int , default=1, help='Min Allele Count Threshold for EA-Pathway Analysis')
    parser.add_argument('--maxAC',type=int , default=5, help='Max Allele Count Threshold for EA-Pathway Analysis')
    parser.add_argument('--cores', type=int, default=1, help='number of CPUs to use for multiprocessing')
    parser.add_argument('--writedata',type=int, default=0, help='keep design matrix after analysis')
    parser.add_argument('--pipeline', default='All', choices=('All', 'ML', 'Pathways', 'EAML', 'EPI', 'Wavelet'),help='which pipeline to be run')

    return parser.parse_args()


def main(args):
    
    if args.ref=='hg19':
        ref = pd.read_csv('./refs/hg19-refGene.protein-coding.txt', delimiter='\t', header=0, index_col='name2')
    elif args.ref=='hg38':
        ref = pd.read_csv('./refs/hg38-refGene.protein-coding.txt', delimiter='\t', header=0, index_col='name2')
    text_file = open(args.savepath+'arguments.txt', 'w')
    text_file.write('You are running '+str(args.pipeline)+' pipeline for:'+'\n'+str(args.minaf)+'<Allele Frequency<'+str(args.maxaf)+ '\n'+str(args.transcript)+' transcript'+ '\n'+ 
                    'Genome Reference:'+str(args.ref))
    text_file.close()
    print('You are running', args.pipeline,'pipeline for:\n',args.minaf,'<Allele Frequency<',args.maxaf, '\n', args.transcript,'transcript \n', 'Genome Reference:',args.ref)
    
    if args.pipeline=='All' or args.pipeline=='Pathways':
##code for parsing VCF into EA-Pathways input file
        sample_file = pd.read_csv(args.samples, index_col=0,header=None)
        cases = sample_file[sample_file.iloc[:,0]==1].index.tolist()
        conts = sample_file[sample_file.iloc[:,0]==0].index.tolist()
        os.makedirs(args.savepath+'Pathway_output', exist_ok = True)
        Pathway_output_path =  args.savepath+'Pathway_output/'
        vcf_parser(args.VCF, Pathway_output_path, args.refPop, args.minAC, args.maxAC, cases, conts)
        print('Prased VCF for EAPathway is completed')
    
        Reactome_input_df = pd.read_csv('./refs/Reactome2020_Greater5Less100_02022021.csv', header=None)
        STRING_input_df = pd.read_csv('./refs/STRINGv11_Greater5Less100_02022021.csv', header=None)    
        for ac in range(args.minAC, args.maxAC + 1):
            # ac+=1
            sample_input_df_Cases = pd.read_csv(args.savepath+'Pathway_output/'+'Input_files/' + 'Cases_PathwaysInput_AC' + str(ac) + '.csv', header=0)
            sample_input_df_Controls = pd.read_csv(args.savepath+'Pathway_output/'+'Input_files/' + 'Controls_PathwaysInput_AC' + str(ac) + '.csv', header=0)
            os.makedirs(Pathway_output_path + 'AC' + str(ac), exist_ok = True)
            output_dir = Pathway_output_path + 'AC' + str(ac)+'/'
            EA_Pathway_Wrapper(sample_input_df_Cases, Reactome_input_df, output_dir, 'Reactome','Cases', args.cores)
            EA_Pathway_Wrapper(sample_input_df_Controls, Reactome_input_df, output_dir, 'Reactome','Controls', args.cores)
            EA_Pathway_Wrapper(sample_input_df_Cases, STRING_input_df, output_dir, 'STRING','Cases', args.cores)
            EA_Pathway_Wrapper(sample_input_df_Controls, STRING_input_df, output_dir, 'STRING','Controls', args.cores) 
    
    if args.pipeline=='All' or args.pipeline=='ML' or args.pipeline=='EAML':
## EA-ML Analysis 
        print('\ nEAML analysis started')
        os.makedirs(args.savepath+'tmp', exist_ok = True)
        os.makedirs(args.savepath+'EAML_output', exist_ok = True)
        EAML_output_path =  args.savepath+'EAML_output/'
        gene_results = Parallel(n_jobs=args.cores)(delayed(eval_gene)(gene,reference=args.ref, data_fn=args.VCF, 
                                                                      targets_fn=args.samples, expdir=args.savepath, write_data=args.writedata, min_af=args.minaf, max_af=args.maxaf, af_field='AF', EA_parser=args.transcript,seed=111,cv=10,weka_path='./weka-3-8-5') for gene in tqdm(ref.index.unique()))
        raw_results = gene_results
        full_results,nonzero_results = report_results(raw_results,EAML_output_path)
        print('\n EAML analysis completed')
    
    if args.pipeline=='All' or args.pipeline=='ML' or args.pipeline=='Wavelet' or args.pipeline=='EPI':   
#### Generate pEA matrix for EA-Wavelet and EPIMUTESTR Analysis
        sample_file = pd.read_csv(args.samples, index_col=0,header=None)
        cases = sample_file[sample_file.iloc[:,0]==1].index.astype(str).tolist()
        conts = sample_file[sample_file.iloc[:,0]==0].index.astype(str).tolist()
        samples = cases + conts
        gene_dfs = Parallel(n_jobs=args.cores)(delayed(parse_gene)(vcf_fn=args.VCF, gene=gene,gene_reference=ref.loc[gene],samples=samples, min_af=args.minaf, max_af=args.maxaf, af_field='AF', EA_parser=args.transcript) for gene in tqdm(ref.index.unique()))
        design_matrix = pd.concat(gene_dfs, axis=1) 
        design_matrix.to_csv(args.savepath+'input_matrix.csv', header=True, index=True)    
## EAWavelet Analysis   
        data, genelist, caselist, contlist, shift, case, cont = init(args.savepath+'input_matrix.csv', args.samples)
        if args.pipeline!='EPI':           
            print('\n EA-Wavelet Analysis Started') 
            case, cont, genelistnet = create_graphs(data, case, cont, caselist, contlist, genelist)
            TotChi, shift, case_node_map, cont_node_map = run_wavelet_decomp(case, cont, genelistnet, shift)
            distances = PCA_distance(TotChi, case_node_map, cont_node_map, shift)
            dffdr = fdrlistgen(data=distances, thr=0.1)
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
