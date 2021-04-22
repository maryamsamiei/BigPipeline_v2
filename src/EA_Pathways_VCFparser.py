#!/usr/bin/env python3

"""
@author: Jenn Asmussen
Pilot script to extract variant information from VCF for running EA-Pathways
"""
import pandas as pd 
import numpy as np 
from pysam import VariantFile
import csv
import time
import  os


def convert_UKB_AF(x):
    try:
        item_new = float(x)
    except ValueError:
        item_new = 1.1
    return item_new

def getRefPopVariants(refPopVariantInfoFile):
    col_names = ['chr', 'pos', 'ref', 'alt', 'ref_AC', 'ref_AF']
    col_type = {'chr': str, 'pos': str, 'ref': str, 'alt': str, 'ref_AC': int, 'ref_AF': str}
    refVariant_df = pd.read_csv(refPopVariantInfoFile, sep='\t', names=col_names, dtype=col_type)

    refVariant_df['identifier'] = refVariant_df['chr'] + '-' + refVariant_df['pos'] + '-' + \
                                         refVariant_df['ref'] + '-' + refVariant_df['alt']
    #refVariant_df['ref_AF_updated'] = refVariant_df['ref_AF'].apply(convert_UKB_AF)
    refVariant_df.drop_duplicates('identifier', keep=False, inplace=True)
    return refVariant_df

def getRefVariantsUnderThreshold(refVariant_df,refPopVariant_lowerThreshold,refPopVariantThreshold):
    refVariantThreshold_df = refVariant_df.loc[(refVariant_df['ref_AC'] <= int(refPopVariantThreshold)) & (refVariant_df['ref_AC'] != 0) &(refVariant_df['ref_AC'] >= int(refPopVariant_lowerThreshold))]
    refVariantThreshold_dict = dict(zip(refVariantThreshold_df.identifier, refVariantThreshold_df.ref_AC))
    return refVariantThreshold_dict

def selectTranscriptSubEA(transcript, sub, EA, gene):
    if type(gene) == tuple:
        final_gene = gene[0]
    else:
        final_gene = gene
    final_EA_sub_transcript = [EA[0], sub[0], transcript[0], final_gene]
    return final_EA_sub_transcript

def variant_class(final_EA):
    if final_EA == 'STOP':
        variant_class = 'stopgain SNV'
    elif final_EA == 'fs-indel':
        variant_class = 'fs-indel'
    elif final_EA == 'indel':
        variant_class = 'indel'
    elif final_EA == 'no_STOP':
        variant_class = 'stop loss'
    elif final_EA == 'silent':
        variant_class = 'synonymous SNV'
    elif final_EA == 'START_loss':
        variant_class = 'start loss'
    else:
        final_EA_float = float(final_EA)
        variant_class = 'nonsynonymous SNV'
    return variant_class

def createFinalVariantMatrix(parsedVCFVariantsMatrix, refPopVariantThreshold_dict):
    snvs_to_filter = ['.', 'UNKNOWN']
    EA_variants_to_drop = ['no_trace', '.', 'no_gene', 'no_action']
    EA_Reactome_dict = {'silent': '', 'STOP': '', 'fs-indel': '', 'indel': '', 'no_STOP': '', 'START_loss': ''}

    parsedVCFVariantsMatrixCleaned = parsedVCFVariantsMatrix.copy()
    parsedVCFVariantsMatrixCleaned = parsedVCFVariantsMatrixCleaned.loc[~parsedVCFVariantsMatrixCleaned.gene.isin(snvs_to_filter)]
    parsedVCFVariantsMatrixCleaned['final_EA_sub_transcript'] = parsedVCFVariantsMatrixCleaned.apply(lambda x: selectTranscriptSubEA(x['NM'], x['sub'], x['EA'], x['gene']), axis=1)
    parse_cols = parsedVCFVariantsMatrixCleaned['final_EA_sub_transcript'].apply(pd.Series)
    parse_cols = parse_cols.rename(columns=lambda x: 'value_' + str(x))
    parsedVCFVariantsMatrixCleaned_final = pd.concat([parsedVCFVariantsMatrixCleaned[:], parse_cols[:]], axis=1)
    parsedVCFVariantsMatrixCleaned_final.rename(columns={'value_0': 'Final_EA', 'value_1': 'Final_Sub', 'value_2': 'Final_Transcript', 'value_3':'Final_Gene'}, inplace=True)
    parsedVCFVariantsMatrixCleaned_final = parsedVCFVariantsMatrixCleaned_final.loc[~parsedVCFVariantsMatrixCleaned_final['Final_EA'].isin(EA_variants_to_drop)]
    parsedVCFVariantsMatrixCleaned_final['Variant_classification'] = parsedVCFVariantsMatrixCleaned_final['Final_EA'].apply(variant_class)
    parsedVCFVariantsMatrixCleaned_final.rename(columns={'Final_Gene': 'gene_ID', 'Final_Sub': 'AAchange', 'Final_EA': 'Action'}, inplace=True)
    parsedVCFVariantsMatrixCleaned_final.replace({'Action': EA_Reactome_dict}, inplace=True)
    parsedVCFVariantsMatrixCleaned_final['identifier'] = parsedVCFVariantsMatrixCleaned_final['chr'] + '-' + \
                                                         parsedVCFVariantsMatrixCleaned_final['pos'] + '-' + \
                                                         parsedVCFVariantsMatrixCleaned_final['ref'] + '-' + \
                                                         parsedVCFVariantsMatrixCleaned_final['alt']
    parsedVCFVariantsMatrixCleaned_final['refPop_AC'] = parsedVCFVariantsMatrixCleaned_final['identifier'].map(refPopVariantThreshold_dict)

    return parsedVCFVariantsMatrixCleaned_final

def selectCanonicalNMIDtranscript(variant_matrix):

    def transcriptSet(v):
        return list(set(v))

    def selectCanonNMID(v):
        new_values = []
        for item in v:
            new_v = item[3:len(item)]
            new_values.append(new_v)
        min_value = min([int(x) for x in new_values])
        min_value = str(min_value)
        canonical = [x for x in new_values if min_value in x]
        canonical = 'NM_' + canonical[0]
        return canonical

    gene_transcriptID_dict = {k: g["Final_Transcript"].tolist() for k, g in variant_matrix.groupby("gene_ID")}
    gene_transcriptID_dict = {k: transcriptSet(v) for k, v in gene_transcriptID_dict.items()}
    gene_transcriptID_dict = {k: selectCanonNMID(v) for k, v in gene_transcriptID_dict.items()}

    finalNMIDs = list(gene_transcriptID_dict.values())
    variant_matrix_final = variant_matrix.loc[variant_matrix['Final_Transcript'].isin(finalNMIDs)]

    return variant_matrix_final


def createACoutputFiles(final_variant_can_df,outputPath, refPopVariant_lowerThreshold,refPopVariantThreshold, cases, controls):
    os.makedirs(outputPath+'Input_files', exist_ok = True)
    InputPath=outputPath+'Input_files/'

    for ac in range(int(refPopVariant_lowerThreshold), int(refPopVariantThreshold) + 1):
        if ac == 0:
            pass
        else:
            final_variant_can_df_AC = final_variant_can_df.copy()
            final_variant_can_df_AC = final_variant_can_df_AC.loc[
                final_variant_can_df_AC['refPop_AC'] <= int(ac)]

            final_case_df = final_variant_can_df_AC.loc[final_variant_can_df_AC['sample'].isin(cases)]
            final_control_df = final_variant_can_df_AC.loc[final_variant_can_df_AC['sample'].isin(controls)]

            final_case_df[['gene_ID', 'Variant_classification', 'AAchange', 'Action', 'sample', 'refPop_AC']].to_csv(
                InputPath  + 'Cases_PathwaysInput_AC' + str(ac) + '.csv',
                index=False)
            final_control_df[['gene_ID', 'Variant_classification', 'AAchange', 'Action', 'sample', 'refPop_AC']].to_csv(
                InputPath  + 'Controls_PathwaysInput_AC' + str(ac) + '.csv',
                index=False)


def vcf_parser(vcfFile, outputPath, refPopVariantFile, refPopVariant_lowerThreshold,refPopVariantThreshold, cases, controls):
    allPatients = cases+controls
    refPopVariant_df = getRefPopVariants(refPopVariantFile)
    print('Number of variants in ref population:', refPopVariant_df.shape[0])

    refPopVariantThreshold_dict = getRefVariantsUnderThreshold(refPopVariant_df, refPopVariant_lowerThreshold,refPopVariantThreshold)
    print(str(refPopVariant_lowerThreshold),'<=Number of variants with AC <=',str(refPopVariantThreshold),':',len(refPopVariantThreshold_dict))
    vcf = VariantFile(vcfFile)
    start = time.time()
    rows = []

    for var in vcf:
        vcfVariantID = str(var.chrom) + '-' +  str(var.pos) + '-' + str(var.ref) + '-' + str(var.alts[0])
        if vcfVariantID in refPopVariantThreshold_dict:
            for sample in allPatients:
                gt = var.samples[str(sample)]['GT']
                if 1 in gt and '.' not in gt:
                    rows.append([var.chrom, var.pos, var.ref, var.alts[0], sample, var.info['gene'],
                                var.info['NM'], var.info['sub'], var.info['EA'], var.samples[str(sample)]['GT']])
                else:
                    pass

    cols = ['chr','pos','ref','alt','sample','gene','NM','sub','EA','GT']
    col_type = {'chr': str, 'pos': str, 'ref': str, 'alt': str}
    parsedVCFforReactomes_df = pd.DataFrame(rows, columns = cols)
    parsedVCFforReactomes_df = parsedVCFforReactomes_df.astype(col_type)

    parsedVCFforReactomes_final_df = createFinalVariantMatrix(parsedVCFforReactomes_df,refPopVariantThreshold_dict)

    parsedVCFforReactomes_canfinal_df = selectCanonicalNMIDtranscript(parsedVCFforReactomes_final_df)

    createACoutputFiles(parsedVCFforReactomes_canfinal_df, outputPath, refPopVariant_lowerThreshold, refPopVariantThreshold, cases, controls)

    print('Time to parse and prep Reactome variants:', time.time() - start)
    
