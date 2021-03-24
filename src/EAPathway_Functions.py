"""
@author: Jenn Asmussen
Pilot script to extract variant information from VCF for running EA-Pathways
"""
import numpy as np
import pandas as pd
#from statistics import stdev
from statsmodels.stats.multitest import fdrcorrection
from scipy import stats as statistics
from scipy.optimize import curve_fit
import time
import random
import os
import multiprocessing as mp
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, AutoMinorLocator
import statistics


### EA_Pathways_Step1_Stats
def sample_and_groups_basic_stats(sample_input, groups_input, output_text_location, variant_types):

    start = time.time()
    text_file = open(output_text_location, 'w')

    # cancer stats
    relevant_variant_lst = variant_types

    # variant stats
    total_mutations = sample_input.shape[0]
    sample_variant_summary = sample_input.groupby('Variant_classification').count()

    total_relevant_mts = 0
    sample_summary_dict = {}
    for item in relevant_variant_lst:
        item_filter = sample_input['Variant_classification'] == item
        num_item = len(list(sample_input.loc[item_filter]['gene_ID']))
        sample_summary_dict[item] = num_item
        total_relevant_mts = total_relevant_mts + num_item

    # stats on biological groupings
    all_groups_genes = []
    for row in range(groups_input.shape[0]):
        group_genes = list(groups_input.loc[row, 2:groups_input.shape[1]].dropna())
        all_groups_genes.extend(group_genes)

    all_groups_genes_unique = set(all_groups_genes)

    # Output summary of stats to text file
    # Write cancer stats to text file
    text_file.write("Summary of Mutations in Input Samples:" + '\n' + '\n')
    #text_file.write("Total number of unique input samples: " + str(unique_patients) + '\n')
    text_file.write("Total number of mutations: " + str(total_mutations) + '\n')
    #text_file.write("Avg. number mutations per sample: " + str(avg_mt_per_patient) + '\n')
    #text_file.write("Stdev. of mutations per sample: " + str(std_mt_per_patient) + '\n')
    #text_file.write("Range of mutations per sample: " + str(min(lst_mt_per_patient))
                    #+ ' - ' + str(max(lst_mt_per_patient)) + '\n')
    text_file.write("Considered (relevant) SNVs: " + '\n')
    for item in relevant_variant_lst:
        text_file.write(item + '\n')
    text_file.write("Total number of relevant SNVs in input samples: "
                    + str(total_relevant_mts) + '\n')
    text_file.write('\n')
    text_file.write("Summary of relevant SNV data for input samples: " + '\n')
    for key, values in sample_summary_dict.items():
        text_file.write(key + ': ' + str(sample_summary_dict[key]) + '\n')

    #write stats on biological groupings
    text_file.write('\n' + '\n')
    text_file.write("Summary of Input Samples Mutations in Biological Groupings:" + '\n' + '\n')
    text_file.write('Total number of biological groups: ' + str(groups_input.shape[0]) + '\n')
    text_file.write("Number of unique genes in biological groups: " + str(len(all_groups_genes_unique)) + '\n')

    cancer_variant_gene_lst = []
    for variant in relevant_variant_lst:
        variant_filter = sample_input['Variant_classification'] == variant
        genes = set(list(sample_input[variant_filter]['gene_ID']))
        text_file.write('Total number of input samples ' + variant + 's in biological groups: '+
              str(len(genes.intersection(all_groups_genes_unique))) + '\n')
        cancer_variant_gene_lst.extend(genes.intersection(all_groups_genes_unique))

    #print("Time to generate basic stats: ", str(time.time() - start))
    text_file.write('Total input samples (relevant) SNVs in biological groups: ' + str(len(cancer_variant_gene_lst)) + '\n')
    text_file.write('Percentage input samples (relevant) SNVs in biological groups: ' + str(len(cancer_variant_gene_lst)/total_relevant_mts) + '\n')
    text_file.write('\n')
    text_file.write("Time to generate basic stats: " + str(time.time() - start))
    text_file.close()

    return all_groups_genes_unique, all_groups_genes

### EA_Pathways_Step2_Gene_x_EA_matrix
### # sample input genes used to build matrix: removed nonsyn SNV without EA scores, kept synon SNV
def generate_sample_gene_by_EAmatrix(user_defined_variants, sample_input, summary_txt_location):

    start = time.time()

    relevant_variant_filter = sample_input['Variant_classification'].isin(user_defined_variants)
    sample_input_filtered = sample_input.copy()
    sample_input_filtered = sample_input_filtered[relevant_variant_filter]

    sample_input_filtered.loc[sample_input_filtered['Variant_classification'] == 'stopgain SNV','Action'] = 100
    sample_input_filtered.loc[sample_input_filtered['Variant_classification'] == 'stop loss','Action'] = 100
    sample_input_filtered.loc[sample_input_filtered['Variant_classification'] == 'start loss','Action'] = 'synon'
    sample_input_filtered.loc[sample_input_filtered['Variant_classification'] == 'fs-indel', 'Action'] = 100
    sample_input_filtered.loc[sample_input_filtered['Variant_classification'] == 'synonymous SNV','Action'] = 'synon'
    sample_input_filtered.loc[sample_input_filtered['Variant_classification'] == 'indel', 'Action'] = 'synon'
    sample_input_filtered['Action'].fillna('no_EA', inplace = True)
    nonsyn_error_genes = list(sample_input_filtered[sample_input_filtered['Action'] == 'no_EA']['gene_ID'])

    EA_100 = [x for x in list(sample_input_filtered['Action']) if x == 100]
    syn_0 = [x for x in list(sample_input_filtered['Action']) if x == 'synon']
    nonsyn_missing_EA = [x for x in list(sample_input_filtered['Action']) if x == 'no_EA']
    nonsyn_with_EA = [x for x in list(sample_input_filtered['Action']) if type(x) != str and x < 100]

    # need to drop from sample_input_filtered df genes with 'no_EA'; these are "errored" genes in AK code
    nonsyn_missing_EA_filter = sample_input_filtered['Action'] != 'no_EA'
    sample_input_filtered_final = sample_input_filtered.copy()
    sample_input_filtered_final = sample_input_filtered_final[nonsyn_missing_EA_filter]

    # creating nested list of genes with corresponding EA scores to buid gene x EA matrix (after filtering 'no_EA' genes)
    # collecting list of all EA scores in cancer type (after filtering 'no_EA' genes)
    # building gene x EA matrix

    sample_input_genes = list(set(list(sample_input_filtered_final['gene_ID'])))

    sample_genes_with_EA_lsts = []
    sample_genes_all_EA_scores = []

    for gene in sample_input_genes:
        gene_filter = sample_input_filtered_final['gene_ID'] == gene
        gene_df = sample_input_filtered_final.copy()
        gene_df = gene_df[gene_filter]
        gene_action = list(gene_df['Action'])
        sample_genes_all_EA_scores.extend(gene_action)
        gene_action.insert(0, gene)
        sample_genes_with_EA_lsts.append(gene_action)

    length = max(map(len, sample_genes_with_EA_lsts))
    sample_input_gene_action_matrix = np.array([x+[None]*(length-len(x)) for x in sample_genes_with_EA_lsts])

    # text_file summary of gene x EA matrix generation process
    text_file = open(summary_txt_location, 'w')

    text_file.write("Summary of sample_input genes x EA matrix creation process:" + '\n' + '\n')
    text_file.write('Shape of sample_input matrix: ' + str(sample_input.shape) + '\n')
    text_file.write("Considered (relevant) SNVs: " + '\n')
    for item in user_defined_variants:
        text_file.write(item + '\n')
    text_file.write('\n')
    text_file.write('Number of stopgain SNV, fs-indels, and stop loss SNV annotated with EA = 100: '+ str(len(EA_100)) + '\n')
    text_file.write('Number of synonymous SNV, indels, and start loss SNV annotated with EA = synon: '+ str(len(syn_0)) + '\n')
    text_file.write('Number of nonsynonymous SNV w/o EA annotated with EA = no_EA: '+ str(len(nonsyn_missing_EA))+ '\n')
    text_file.write('Number of nonsynonymous SNV with EA scores: '+ str(len(nonsyn_with_EA)) + '\n')
    text_file.write('Shape of input sample matrix (relevant SNVs only): '+ str(sample_input_filtered.shape)+ '\n')
    text_file.write('Shape of input sample matrix (relevant SNVs only + drop no_EA SNVs): '+ str(sample_input_filtered_final.shape)+ '\n')
    text_file.write('\n')
    text_file.write('Total number of unique genes with SNVs in input samples after filtering: ' + str(len(sample_input_genes))+ '\n')
    text_file.write('Shape of gene x EA matrix: ' + str(sample_input_gene_action_matrix.shape) + '\n')
    text_file.write('\n')
    text_file.write("Time to generate gene x EA matrix: " + str(time.time() - start) + '\n')

    text_file.close()
    
    return sample_input_gene_action_matrix, sample_genes_all_EA_scores, sample_input_genes, list(set(nonsyn_error_genes))


### EA_Pathways_Step3_KS_Test_SingleGenes

from scipy import stats as statistics
def KS_test_individual_sample_genes(sample_genes_all_EA_scores, sample_input_gene_action_matrix, txt_summary_location, csv_summary_location):
    
    # v3: KS test on each mutated gene in cancer type; synon SNVs included in gene x EA matrix
    start = time.time()
    sample_genes_all_float_EA_scores = [x for x in sample_genes_all_EA_scores if x != 'synon']
    sample_genes_all_float_EA_scores = [float(x) for x in sample_genes_all_float_EA_scores]

    sig_single_gene_df = pd.DataFrame()
    sig_single_gene_df['gene'] = sample_input_gene_action_matrix[:,0]

    sig_single_gene_pvalues = []

    for row in range(sample_input_gene_action_matrix.shape[0]):
        gene_ea_scores = sample_input_gene_action_matrix[row, 1:sample_input_gene_action_matrix.shape[1]]
        gene_ea_scores = [x for x in gene_ea_scores if x != None]
        gene_ea_scores = [x for x in gene_ea_scores if x != 'synon']
        gene_ea_scores = [float(x) for x in gene_ea_scores]
        if len(gene_ea_scores) == 0:
            sig_single_gene_pvalues.append('No EA scores')
        else: 
            sig_single_gene_pvalues.append(statistics.mstats.ks_twosamp(gene_ea_scores, sample_genes_all_float_EA_scores, alternative='less')[1])

    sig_single_gene_df['p_value'] = sig_single_gene_pvalues

    sig_single_gene_pvalues_dropNoEA = [x for x in sig_single_gene_pvalues if x != 'No EA scores']

    q_value_df = sig_single_gene_df.copy()
    index_noIntEA = q_value_df[q_value_df['p_value']=='No EA scores'].index
    q_value_df.drop(index_noIntEA, inplace = True)
    q_value_df_pvalue_lst = list(q_value_df['p_value'])
    sig_single_gene_qvalues = fdrcorrection(q_value_df_pvalue_lst, alpha=0.05)[1]
    q_value_df['q_value'] = sig_single_gene_qvalues

    ks_merged_df = sig_single_gene_df.merge(q_value_df, how = 'outer', left_on='gene', right_on='gene')
    ks_merged_df.sort_values(by = 'p_value_y', inplace = True)
    ks_merged_df.drop(columns = 'p_value_y', inplace = True)
    ks_merged_df.rename(columns = {'p_value_x':'p_value'}, inplace = True)

    # generate sig_single_genes_lst from KS analysis performed on each gene
    filter_sig_genes = ks_merged_df['q_value'] < 0.1
    sig_single_genes_lst = list(ks_merged_df[filter_sig_genes]['gene'])

    ks_merged_df.to_csv(csv_summary_location, index = False)

    # text file summarizing KS test of each gene mutated in cancer type
    text_file = open(txt_summary_location, 'w')
    text_file.write("Summary of single gene KS tests for input samples:" + '\n' + '\n')
    text_file.write('Total number of EA annotations in input samples: ' +  str(len(sample_genes_all_EA_scores)) + '\n')
    text_file.write('Total number of integer EA scores in input samples: ' + str(len(sample_genes_all_float_EA_scores))+ '\n')
    text_file.write('Total number of unique mutated genes assessed by KS test: '+ str(len(sig_single_gene_pvalues_dropNoEA))+ '\n')
    text_file.write('Shape of gene x EA matrix: ' + str(sample_input_gene_action_matrix.shape) + '\n')
    text_file.write('Shape of sig_single_gene_matrix: '+ str(ks_merged_df.shape) + '\n')
    text_file.write('\n')
    text_file.write("Time to perform KS test for each unique mutated gene in input samples: " + str(time.time() - start) + '\n')
    text_file.close()
    
    return sample_genes_all_float_EA_scores, sig_single_genes_lst

### /EA_Pathways_Step4_PrepSamples4LOO_Analysis
### #v3 - Nested lists of groups, genes, and EA scores (includes synon SNV) 
## build two sets of nested lists, one with sig genes and one without sig genes
def PrepSamples4LOO_Analysis(sample_input_genes, all_groups_genes_unique, nonsyn_errored_genes, groups_input, sample_input_gene_action_matrix, sig_single_genes_lst, txt_summary_location):

    #first part of function - prep groups with EA scores for LOO analysis
    start=time.time()
    overlapping_sample_and_group_genes = set(sample_input_genes).intersection(all_groups_genes_unique)
    groups_input_id_lst = list(groups_input[0])

    groups_with_genes_and_EA_scores_lst = []
    groups_with_noSigGenes_and_EA_scores_lst = []
    genes_in_cohort_and_groups_with_mutation = []
    nonsyn_error_genes_in_cohort_and_groups = []
    #total_variants_in_group = []

    for row_outer in range(groups_input.shape[0]):
        #omit from analysis any groups with only one gene
        if len(list(groups_input.iloc[row_outer, 2:groups_input.shape[1]].dropna())) == 1:
            pass
        else:
            group_lst_all = []
            group_lst_noSigGenes = []
            group_lst_error_genes = []
            #group_lst_variants = []

            group_name = groups_input.iloc[row_outer,0]
            group_lst_all.append(group_name)
            group_lst_noSigGenes.append(group_name)
            #group_lst_variants.append(group_name)
            #group_lst_error_genes.append(group_name)

            group_genes = list(groups_input.iloc[row_outer, 2:groups_input.shape[1]].dropna())

            for gene in group_genes:
                gene_with_EA_lst = []
                gene_with_EA_lst.append(gene)

                for row_inner in range(sample_input_gene_action_matrix.shape[0]):
                    if gene == sample_input_gene_action_matrix[row_inner,0]:
                        gene_ea_scores = sample_input_gene_action_matrix[row_inner, 1:sample_input_gene_action_matrix.shape[1]]
                        gene_ea_scores = [x for x in gene_ea_scores if x != None]
                        gene_with_EA_lst.extend(gene_ea_scores)

                        if len(gene_ea_scores) != 0:
                            genes_in_cohort_and_groups_with_mutation.append(gene)
                        else:
                            pass

                group_lst_all.append(gene_with_EA_lst)

                if gene in sig_single_genes_lst:
                    pass
                else:
                    group_lst_noSigGenes.append(gene_with_EA_lst)

                if gene in nonsyn_errored_genes:
                    group_lst_error_genes.append(gene)
                else:
                    pass

            groups_with_genes_and_EA_scores_lst.append(group_lst_all)
            groups_with_noSigGenes_and_EA_scores_lst.append(group_lst_noSigGenes)
            nonsyn_error_genes_in_cohort_and_groups.append(group_lst_error_genes)

    #second part of function - generate summary matrix to collect analysis results
    #v3 - generate summary matrix prior to LOO analysis, synon SNV included in count

    #create lists to generate summary matrix
    group_names_lst = []
    group_length_lst = []
    group_genes_with_EAscores_lst = [] #v3 count includes synon SNV
    group_number_genes_with_EAscores_lst = [] #v3 count includes synon SNV
    group_sig_single_genes_lst = []
    group_total_variants_in_cohort_lst = []

    for group in groups_with_genes_and_EA_scores_lst:
        group_name = group[0]
        group_names_lst.append(group_name)

        group_length = len(group) - 1
        group_length_lst.append(group_length)

        group_genes_with_EAscores = []
        group_sig_genes = []
        group_variants = []

        for gene in group[1:]:
            if len(gene) == 1:
                pass
            else:
                group_genes_with_EAscores.append(gene[0])
                gene_variants = gene[1:len(gene)]
                gene_variants = [x for x in gene_variants if x != 'synon']
                group_variants.append(len(gene_variants))

            if gene[0] in sig_single_genes_lst:
                group_sig_genes.append(gene[0])
            else:
                pass

        group_genes_with_EAscores_lst.append(group_genes_with_EAscores)
        group_number_genes_with_EAscores_lst.append(len(group_genes_with_EAscores))
        group_sig_single_genes_lst.append(group_sig_genes)
        group_total_variants_in_cohort_lst.append(sum(group_variants))

    summary_matrix = pd.DataFrame()
    summary_matrix['group_name'] = group_names_lst
    summary_matrix['number_group_genes'] = group_length_lst

    group_sig_single_genes_lst_lengths = [len(x) for x in group_sig_single_genes_lst]
    functional_group_size_lst = [a - b for a, b in zip(group_length_lst, group_sig_single_genes_lst_lengths)]
    #summary_matrix['functional_group_size'] = functional_group_size_lst

    group_error_genes_lst_lengths = [len(x) for x in nonsyn_error_genes_in_cohort_and_groups]
    functional_group_size_lst2 = [a - b for a, b in zip(functional_group_size_lst, group_error_genes_lst_lengths)]
    summary_matrix['functional_group_size'] = functional_group_size_lst2
    summary_matrix['group_errored_genes'] = nonsyn_error_genes_in_cohort_and_groups


    summary_matrix['number_group_genes_with_EAscores'] = group_number_genes_with_EAscores_lst
    summary_matrix['group_genes_with_EAscores'] = group_genes_with_EAscores_lst
    summary_matrix['group_sig_genes'] = group_sig_single_genes_lst
    summary_matrix['total_group_variants'] = group_total_variants_in_cohort_lst

    text_file = open(txt_summary_location, 'w')
    text_file.write('Summary of Prepping Input Samples into Biological Groups for LOO Analysis' + '\n' + '\n')
    text_file.write('Number of biological groups: ' + str(len(groups_input_id_lst)) + '\n')
    text_file.write('Number of unique genes in biological groups: ' + str(len(all_groups_genes_unique)) + '\n')
    text_file.write('Number of unique genes with at least one mutation in input samples: ' + str(len(sample_input_genes)) + '\n')
    text_file.write('Overlap of unique sample input genes and unique biological groups genes: ' +
          str(len(overlapping_sample_and_group_genes)) + '\n')
    text_file.write('Number of genes in input samples and groups with >= 1 EA annotation: ' +
          str(len(set(genes_in_cohort_and_groups_with_mutation))) + '\n' + '\n')
    text_file.write('Time to prep sample for LOO analysis and generate summary matrix: ' + str(time.time() - start) + '\n')
    text_file.close()

    return summary_matrix, groups_with_noSigGenes_and_EA_scores_lst


### EA_Pathways_Step5_Sample_LOO_Analysis
### v3 - version corrects for synon SNV
### function will be distributed to different cores with list item
### currently written to take one item from nested list and compare to background EA scores
### output group name, original group KS test p-value, core genes, core genes KS test p-value
def input_samples_loo_multiprocessing(arg):
    simulation = arg[0]
    background = arg[1]
    name, ini_pvalue, core_gen, cor_pvalue = group_LOO_core_gene_analysis(simulation, background)
    return name, ini_pvalue, core_gen, cor_pvalue

def pool_loo_analysis_input_samples_fx(all_simulations, EA_background, cores):
    args = tuple(zip(all_simulations, [EA_background] * len(all_simulations)))
    pool = mp.Pool(processes=cores)
    output = pool.map(input_samples_loo_multiprocessing, args)
    pool.close()
    pool.join()
    return output


def group_LOO_core_gene_analysis(biological_group_noSigGenes_EAscores_item, background_EAscores):
    
    background_float_EAscores = [float(x) for x in background_EAscores]
    
    group_name = biological_group_noSigGenes_EAscores_item[0]
    # uncomment the following line if you want to see progress during LOO analysis
    # print(group_name)
    group_length = len(biological_group_noSigGenes_EAscores_item)
    
    # collect group EA scores
    group_all_EA_scores = []
    for gene in biological_group_noSigGenes_EAscores_item[1:group_length]:
        if len(gene) == 1:
            pass
        else:
            group_all_EA_scores.append(gene[1:len(gene)])
    
    group_all_EA_scores = [item for sublist in group_all_EA_scores for item in sublist]
    group_all_EA_scores = [x for x in group_all_EA_scores if x != 'synon']
    group_all_EA_scores = [float(x) for x in group_all_EA_scores]
    
    # perform KS test on group EA scores to generate original p-value
    if len(group_all_EA_scores) == 0:
        initial_group_pvalue = 1
    else:
        initial_group_pvalue = statistics.mstats.ks_twosamp(group_all_EA_scores, background_float_EAscores, alternative='less')[1]
    
    # perform LOO-KS analysis if initial_group_pvalue < 1, else pass
    group_core_genes = []
    group_core_EA_scores = []
    
    if initial_group_pvalue < 1:
        
        for test_gene in biological_group_noSigGenes_EAscores_item[1:group_length]:
            if len(test_gene) > 1:
                test_group = biological_group_noSigGenes_EAscores_item[1:group_length].copy()
                test_group.remove(test_gene)
                
                loo_group_EA_scores = []
                for gene_item in test_group:
                    if len(gene_item) == 1:
                        pass
                    else:
                        loo_group_EA_scores.append(gene_item[1:len(gene_item)])
                
                loo_group_EA_scores = [item for sublist in loo_group_EA_scores for item in sublist]
                loo_group_EA_scores = [x for x in loo_group_EA_scores if x != 'synon']
                loo_group_EA_scores = [float(x) for x in loo_group_EA_scores]
                
                # perform KS test on group subset to generate new p-value
                if len(loo_group_EA_scores) == 0:
                    group_core_genes.append(test_gene[0])
                    group_core_EA_scores.append(test_gene[1:len(test_gene)])
                else:
                    test_gene_loo_pvalue = statistics.mstats.ks_twosamp(loo_group_EA_scores, background_float_EAscores, alternative='less')[1]
                    
                    if test_gene_loo_pvalue > initial_group_pvalue:
                        group_core_genes.append(test_gene[0])
                        group_core_EA_scores.append(test_gene[1:len(test_gene)])
                    else:
                        pass   
            
            else:
                pass
            
    else:
        # generate placeholders for whatever is generated in previous if statement
        # append "no core genes" to list where core genes are being collected
        group_core_genes.append('No nonsyn SNV mutations in biological group')
    
    group_core_EA_scores = [item for sublist in group_core_EA_scores for item in sublist]
    group_core_EA_scores = [x for x in group_core_EA_scores if x != 'synon']
    
    # perform KS test on core gene EA scores to generate final p-value
    if len(group_core_EA_scores) == 0:
        core_group_pvalue = 1
    else:
        core_group_pvalue = statistics.mstats.ks_twosamp(group_core_EA_scores, background_float_EAscores, alternative='less')[1]

    return group_name, initial_group_pvalue, group_core_genes, core_group_pvalue   

### EA_Pathways_Step6_GenerateGroupSimulations
### function includes sig_single_genes when generating sim paths
### function uses all genes in groups (includes duplicates) to generate simulated pathways
### note - AK did not include sig_single_genes when building sim paths
### this code includes sig_single_genes when making simulations, but removes them for LOO analysis
### including sig_single_genes should make this more like real groups
### function generates nested list of simulated pathways with EA scores
### function input: num sims, gene x EA matrix, ALL group input genes, sig_single_genes
def generate_simulated_groups(total_simulations, simulation_path_size, gene_EA_matrix, ALL_group_input_genes, sig_single_genes_lst, errors):
    
    lst_all_genes = ALL_group_input_genes
    # the following line removes sig_single_genes from list of group input genes prior to building simulations
    lst_all_genes = [elem for elem in lst_all_genes if elem not in sig_single_genes_lst]
    # the following line removes errored genes from list of gorup input genes prior to building simulations
    lst_all_genes = [elem for elem in lst_all_genes if elem not in errors]

    lst_sim_pathways = []
    
    for sim in range(total_simulations):
        sim_genes = random.sample(lst_all_genes, k = simulation_path_size)
        if len(sim_genes) != len(set(sim_genes)):
            while len(sim_genes) != len(set(sim_genes)):
                new_sim_genes = []
                set_sim_genes = set(sim_genes)
                num_replacements = simulation_path_size - len(set_sim_genes)
                replacements = random.sample(lst_all_genes, k = num_replacements)
                new_sim_genes.extend(list(set_sim_genes))
                new_sim_genes.extend(replacements)
                sim_genes = new_sim_genes.copy()
            lst_sim_pathways.append(sim_genes)
        else:
            lst_sim_pathways.append(sim_genes)
    
    lst_final_sim_paths_NoSigGenes = []
    i = 1
    for sim_path in lst_sim_pathways:
        final_path = [str(simulation_path_size)+'_sim_'+str(i)]
        i += 1
        for gene in sim_path:
            if gene in sig_single_genes_lst:
                pass
            else:
                gene_name_with_EA_scores = []
                gene_name_with_EA_scores.append(gene)
                for row in range(gene_EA_matrix.shape[0]):
                    if gene_EA_matrix[row,0] == gene:
                        gene_ea_scores = gene_EA_matrix[row, 1:gene_EA_matrix.shape[1]]
                        gene_ea_scores = [x for x in gene_ea_scores if x != None]
                        gene_name_with_EA_scores.extend(gene_ea_scores)
                    else:
                        pass
                final_path.append(gene_name_with_EA_scores)
        lst_final_sim_paths_NoSigGenes.append(final_path)
    return lst_final_sim_paths_NoSigGenes

def build_sims_multiprocessing(arg):
    num_sims = arg[0]
    size_sims = arg[1]
    sample_gene_EA_matrix = arg[2]
    group_input_genes = arg[3]
    sig_genes = arg[4]
    error_genes = arg[5]
    simulations = generate_simulated_groups(num_sims, size_sims, sample_gene_EA_matrix, group_input_genes, sig_genes, error_genes)
    return simulations

def pool_fx_sims(number_simulations, simulation_size_lst, gene_EA_matrix, ALL_group_input_genes, sig_single_genes_lst, error_gene_lst, cores):
    args = tuple(zip(np.full(len(simulation_size_lst), number_simulations).tolist(),
                    simulation_size_lst, [gene_EA_matrix] * len(simulation_size_lst), 
                    [ALL_group_input_genes] * len(simulation_size_lst),
                     [sig_single_genes_lst] * len(simulation_size_lst),
                     [error_gene_lst] * len(simulation_size_lst)))
    pool = mp.Pool(processes=cores)
    output = pool.map(build_sims_multiprocessing, args)
    pool.close()
    pool.join()
    
    return output

### EA_Pathways_Step7_LOO_AnalysisSimGroups
def sims_loo_multiprocessing(arg):
    simulation = arg[0]
    background = arg[1]
    name, ini_pvalue, core_gen, cor_pvalue = group_LOO_core_gene_analysis(simulation, background)
    return name, ini_pvalue, core_gen, cor_pvalue

def pool_loo_analysis_sims_fx(all_simulations, EA_background, cores):
    args = tuple(zip(all_simulations, [EA_background] * len(all_simulations)))
    pool = mp.Pool(processes=cores)
    output = pool.map(sims_loo_multiprocessing, args)
    pool.close()
    pool.join()
    return output

def collect_sim_core_pvalues_and_percentiles(list_all_group_sizes, LOO_KS_sim_output):
    all_core_pvalues_lst = []
    all_core_percentiles_lst = []
    for sim_size in list_all_group_sizes:
        sim_name = "simulations_size_"+str(sim_size)
        sim_core_pvalues = []
        sim_core_pvalues.append(sim_name)
        for item in LOO_KS_sim_output:
            loo_simulation_output_name = item[0]
            result = loo_simulation_output_name.startswith(str(sim_size)+'_')
            if result == False:
                pass
            else:
                test_core_pvalue = item[3]
                sim_core_pvalues.append(test_core_pvalue)
        all_core_pvalues_lst.append(sim_core_pvalues)

        #prep list for conversion to percentile
        sim_core_pvalues_copy = sim_core_pvalues.copy()
        sim_core_pvalues_copy.remove(sim_name)
        sim_core_percentile = []
        sim_core_percentile.append(sim_name)
        for perc in range(100):
            sim_core_percentile.append(np.percentile(sim_core_pvalues_copy, perc))
        all_core_percentiles_lst.append(sim_core_percentile)
    
    return all_core_pvalues_lst, all_core_percentiles_lst

### EA_Pathways_Step8_CompareSamplesToSims

### define functions needed for threshold simulation
def func(x, a, b):
    return a * np.exp(-b * x)

def rsquared (x, y):
    slope, intercept, r_value, p_value, stderr = statistics.linregress(x, y)
    return r_value ** 2

# define functions to calculate threshold, foldbetter, qvalue filter columns, sig core genes
def qvalue_filter(row):
    if row['fdr_q_value_core_pathway'] < 0.05:
        return 1
    else:
        return 0

def threshold(row, a, b, threshold_dictionary):
    if int(row['functional_group_size']) <= 15 and int(row['functional_group_size']) >= 5: #modified for removing smaller pathway simulations
        thresh = threshold_dictionary[int(row['functional_group_size'])]
        return thresh
    else:
        thresh = func(float(row['functional_group_size']), a, b)
        return thresh

def foldbetter(row):
    fold_better = row['5th_percentile_threshold']/row['core_group_pvalue']
    return fold_better

def collect_core_genes_in_sig_groups(summary_matrix):
    significant_core_genes = []
    for row in range(summary_matrix.shape[0]):
        row_q_value_filter = summary_matrix.at[row, 'passed_q_value_filter']
        row_foldbetter = summary_matrix.at[row, 'fold_better']
        row_total_variants = summary_matrix.at[row, 'total_group_variants']
        
        if row_q_value_filter == 1 and row_foldbetter > 1 and row_total_variants >= 20:
            row_core_genes = [summary_matrix.at[row, 'core_genes']]
            significant_core_genes.extend(row_core_genes)
    
    significant_core_genes = [item for sublist in significant_core_genes for item in sublist]
    significant_core_genes_set = set(significant_core_genes)
    
    return significant_core_genes_set

# function comparing true data to simulated data
def compare_to_simulations(small_and_large_sims_lst, sim_percentile_matrix, summary_matrix):
    # collect 5th percentile p-values from simulations for threshold calculation
    xvalues = small_and_large_sims_lst
    yvalues = list(sim_percentile_matrix[5])
    xvalues = [float(x) for x in xvalues]
    yvalues = [float(x) for x in yvalues]

    # create threshold_dictionary to call exact threshold values for pathways <= 15
    threshold_dictionary = dict(zip(xvalues, yvalues))

    # calculate coefficients for func that was previously defined using x and y values
    x = np.array(xvalues)
    y = np.array(yvalues)
    popt, pcov = curve_fit(func, x, y, [0.2, 0.2])
    print('threshold=' + str(popt[0]) + '*EXP(-' + str(popt[1]) + '*OriginalGroupSize)')
    a = popt[0]
    b = popt[1]

    # calculate r-squared for ln transformed data
    # r-squared is generated by transforming data by natural log to make linear
    # but, actual data is fit to exponential curve and threshold is derived from expoential fitted data
    # should x undergo ln too? check without...
    ln_x = [np.log(i) for i in x]
    ln_y = [np.log(i) for i in y]
    print('r-squared for this equation is: ' + str(rsquared(x, y)))
    
    # had to generate extra steps for fdr correction because i set p-value = 1 for groups with no mutations
    # groups without mutations are removed prior to fdr correction
    # also, noticed that AK code is missing some of the communities...don't know what is happening
    # this is going to generate different q-values, unless i can determine location of discrepency

    # fdr correction of true sample data
    summary_matrix_fdr = summary_matrix.copy()
    summary_matrix_fdr = summary_matrix_fdr[~summary_matrix_fdr['core_genes'].apply(lambda x: 'No nonsyn SNV mutations in biological group' in x)]
    sample_core_group_pvalues = tuple(summary_matrix_fdr['core_group_pvalue'])
    fdr_core_group_pvalues = fdrcorrection(sample_core_group_pvalues, alpha = 0.05)[1]
    summary_matrix_fdr['fdr_q_value_core_pathway'] = fdr_core_group_pvalues
    summary_matrix = summary_matrix.merge(summary_matrix_fdr[['group_name', 'fdr_q_value_core_pathway']], on="group_name", how="outer")

    # compare true sample data to simulations and add to summary_matrix
    summary_matrix['passed_q_value_filter'] = summary_matrix.apply (lambda row: qvalue_filter(row), axis=1)
    summary_matrix['5th_percentile_threshold'] = summary_matrix.apply (lambda row: threshold(row, a, b, threshold_dictionary), axis=1)
    summary_matrix['fold_better'] = summary_matrix.apply (lambda row: foldbetter(row), axis=1)
    summary_matrix.sort_values(by = ['passed_q_value_filter', 'fold_better'], ascending = False, inplace = True)
    
    return summary_matrix

### EA_Pathways_Step9_MakeHistograms
def full_biological_group_historgrams(group_name, group_input_df, sig_genes, sample_gene_EA_matrix, output_directory):

    for row in range(group_input_df.shape[0]):
        if group_input_df.iloc[row,0] == group_name:
            full_group_genes = list(group_input_df.iloc[row, 2:group_input_df.shape[1]].dropna())
        else:
            pass

    full_group_genes = [x for x in full_group_genes if x not in sig_genes]

    full_group_genes_for_plot = []
    full_group_genes_EA_scores = []
    for gene in full_group_genes:
        for row in range(sample_gene_EA_matrix.shape[0]):
            if sample_gene_EA_matrix[row, 0] == gene:
                #print(gene)
                full_group_genes_for_plot.append(gene)
                gene_EA_scores = sample_gene_EA_matrix[row, 1:sample_gene_EA_matrix.shape[1]]
                #print(gene_EA_scores)
                gene_EA_scores = [x for x in gene_EA_scores if x != None]
                gene_EA_scores = ['0' if x == 'synon' else x for x in gene_EA_scores]
                gene_EA_scores = [float(x) for x in gene_EA_scores]
                full_group_genes_EA_scores.append(gene_EA_scores)

    n_bins = 10
    #plt.switch_backend('agg')
    fig, ax = plt.subplots()
    ax.hist(full_group_genes_EA_scores, n_bins, range = (0.0, 100.0), histtype = 'barstacked', stacked = True, label = full_group_genes_for_plot)
    ax.set_xlim(0,100)
    ax.set_ylim(0,30)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    plt.legend(loc = 'upper right')
    plt.title(str(group_name))
    plt.savefig(output_directory + str(group_name) + "_FullGroupHist_NoSigGenes.png", transparent = True)
    plt.close()

def core_biological_group_historgrams(group_name, significant_groups_df, sample_gene_EA_matrix, output_directory):
    for row in range(significant_groups_df.shape[0]):
        if significant_groups_df.iloc[row, 0] == group_name:
            group_name_info = significant_groups_df.iloc[row].copy()
            group_core_genes = group_name_info['core_genes']
            #print(group_core_genes)
            #print(len(group_core_genes))
            #print(type(group_core_genes))

    group_core_genes_EA_scores = []
    for gene in group_core_genes:
        for row in range(sample_gene_EA_matrix.shape[0]):
            if sample_gene_EA_matrix[row, 0] == gene:
                gene_EA_scores = sample_gene_EA_matrix[row, 1:sample_gene_EA_matrix.shape[1]]
                gene_EA_scores = [x for x in gene_EA_scores if x != None]
                gene_EA_scores = ['0' if x == 'synon' else x for x in gene_EA_scores]
                gene_EA_scores = [float(x) for x in gene_EA_scores]
                group_core_genes_EA_scores.append(gene_EA_scores)

    n_bins = 10
    #plt.switch_backend('agg')
    fig, ax = plt.subplots()
    ax.hist(group_core_genes_EA_scores, n_bins, range = (0.0, 100.0), histtype = 'barstacked', stacked = True, label = group_core_genes)
    ax.set_xlim(0,100)
    ax.set_ylim(0,30)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    plt.legend(loc = 'upper right')
    plt.title(str(group_name))
    plt.savefig(output_directory + str(group_name) + "_CoreGenesHist_NoSigGenes.png", transparent = True)
    plt.close()
###----------------------------------------------------------------------------------------------------------

def EA_Pathway_Wrapper(sample_input_df, groups_input_df, output_directory, method, case_control, number_cores):

    # os.makedirs(output_directory+'Hist', exist_ok = True)
    # histFolder = output_directory+'Hist/'
### output file locations
    os.makedirs(output_directory+ method + '_' + case_control, exist_ok = True)

    LogFile_summary_txt_location = output_directory+ method + '_' + case_control + '/Time4EachStep.txt'
    step1_summary_txt_location = output_directory+ method + '_' + case_control + '/Step1_BasicStats.txt'
    step2_summary_txt_location = output_directory+ method + '_' + case_control+ '/Step2_GeneEA_MatrixSummary.txt'
    step3_summary_txt_location = output_directory+ method + '_' + case_control +'/Step3_KS_TestSingleGeneSummary.txt'
    step3_summary_csv_location = output_directory+ method + '_' + case_control + '/Step3_KS_TestSingleGeneSummary.csv'
    step4_summary_txt_location = output_directory+ method + '_' + case_control+ '/Step4_SamplePrep4LOO_AnalysisSummary.txt'
    step7_RawSim_summary_csv_location = output_directory+ method + '_' + case_control + '/Step7_RawSimulation.csv'
    step7_SimCorePvalues_csv_location = output_directory+ method + '_' + case_control + '/Step7_SimCorePvalues.csv'
    step7_SimCorePvaluePerc_csv_location = output_directory+ method + '_' + case_control + '/Step7_SimCorePvaluePercentiles.csv'
    step8_FinalSummary_location = output_directory+ method + '_' + case_control+ '/Step8_FinalAnalysisSummary.csv'
    step8_SigCoreGenes_location = output_directory+ method + '_' + case_control +'/Step8_SigCoreGenesAbove5thPercentile.txt'
    step9_HistCoreGenes_location = output_directory+ method + '_' + case_control+ '/'

    # initiate start time for full analysis
    start_all = time.time()
    # variants considered in analysis
    relevant_variants = ['nonsynonymous SNV', 'stopgain SNV', 'synonymous SNV','stop loss','start loss','indel','fs-indel']

    # open log file to collect time for each step
    LogFile = open(LogFile_summary_txt_location, 'w')

    # Step1: Generate basic stats on samples and biological groups
    step1_start = time.time()
    print('Now Performing Step1: Generating basic statistics on input samples and biological groups')
    all_unique_group_genes, all_group_genes_with_duplicates = sample_and_groups_basic_stats(sample_input_df,
                                                                                        groups_input_df,
                                                                                        step1_summary_txt_location,
                                                                                        relevant_variants)
    LogFile.write('Time to perform Step1: ' + str(time.time() - step1_start) + '\n')

    # Step2: Generate sample genes x EA matrix
    step2_start = time.time()
    print('Now Performing Step2: Creating genes from input samples X EA scores matrix')
    sample_gene_EA_matrix, all_sample_genes_EAscores, all_unique_sample_genes, errored_genes = generate_sample_gene_by_EAmatrix(relevant_variants,
                                                                                                             sample_input_df,
                                                                                                             step2_summary_txt_location)
    LogFile.write('Time to perform Step2: ' + str(time.time() - step2_start) + '\n')

    # Step3: Perform KS test on each sample gene to identify individually significant sample genes
    step3_start = time.time()
    print('Now Performing Step3: Performing KS test on each gene in input samples to identify genes with significantly biased EA distributions')
    all_sample_genes_floatEAscores, sig_single_sample_genes_lst = KS_test_individual_sample_genes(all_sample_genes_EAscores,
                                                                                              sample_gene_EA_matrix,
                                                                                              step3_summary_txt_location,
                                                                                              step3_summary_csv_location)
    LogFile.write('Time to perform Step3: ' + str(time.time() - step3_start) + '\n')

    # Step4: Prep samples genes in biological groups with EA scores for LOO analysis and generate summary matrix to collect results
    step4_start = time.time()
    print('Now Performing Step4: Prepping input samples and biological groups for LOO analysis')
    final_summary_matrix, prepped_groups_noSigGenes_EAscores_lst = PrepSamples4LOO_Analysis(all_unique_sample_genes,
                                                                                        all_unique_group_genes,
                                                                                        errored_genes,
                                                                                        groups_input_df,
                                                                                        sample_gene_EA_matrix,
                                                                                        sig_single_sample_genes_lst,
                                                                                        step4_summary_txt_location)
    #final_summary_matrix.to_csv(step4_TestMatrix_location, index = False)
    LogFile.write('Time to perform Step4: ' + str(time.time() - step4_start) + '\n')

    # Step5: Perform LOO analysis on input samples and biological groups
    print('Now Performing Step5: Performing LOO analysis on input samples and biological groups')
    step5_start = time.time()
    loo_input_samples_output = pool_loo_analysis_sims_fx(prepped_groups_noSigGenes_EAscores_lst,
                                                  all_sample_genes_floatEAscores,
                                                  number_cores)
    loo_input_samples_df = pd.DataFrame(loo_input_samples_output,
                             columns=['group_name', 'original_group_pvalue', 'core_genes', 'core_group_pvalue'])
    final_summary_matrix = final_summary_matrix.merge(loo_input_samples_df, on = "group_name", how = "outer")
    #final_summary_matrix.to_csv(step5_TestMatrix_location, index = False)
    LogFile.write('Time to perform Step5: ' + str(time.time() - step5_start) + '\n')

    # Step6: Generate simulated biological groups with multiprocessing
    step6_start = time.time()
    print('Now Performing Step6: Generating simulated biological groups')
    lst_small_sims = np.arange(5, 16, 1).tolist() #remove simulations smaller than size 5
    lst_large_sims = np.arange(15, 51, 5).tolist()
    lst_large_sims.remove(15)
    small_and_large_sims = lst_small_sims + lst_large_sims

    # build simulations with parallelization
    all_small_sims = pool_fx_sims(1000, lst_small_sims, sample_gene_EA_matrix, all_group_genes_with_duplicates,
                              sig_single_sample_genes_lst, errored_genes, number_cores)

    start6_large_start = time.time()
    all_large_sims_100 = pool_fx_sims(1000, lst_large_sims, sample_gene_EA_matrix, all_group_genes_with_duplicates,
                                  sig_single_sample_genes_lst, errored_genes, number_cores)

    all_small_sims_flat = [item for sublist in all_small_sims for item in sublist]
    all_large_sims_100_flat = [item for sublist in all_large_sims_100 for item in sublist]
    all_sims_combined_flat = all_small_sims_flat + all_large_sims_100_flat
    LogFile.write('Time to perform Step6: ' + str(time.time() - step6_start) + '\n')

    # Step7: Perform LOO analysis on simulated groups
    step7_start = time.time()
    print('Now Performing Step7: LOO analysis on simulated biological groups')
    # LOO analysis with multiprocessing
    loo_simulation_output = pool_loo_analysis_sims_fx(all_sims_combined_flat, all_sample_genes_floatEAscores, number_cores)
    # save simulation LOO and KS analysis output to dataframe for QC purposes
    simulation_df = pd.DataFrame(loo_simulation_output,
                             columns=['sim_number', 'sim_initial_pvalue', 'sim_core_genes', 'sim_core_pvalue'])
    simulation_df.to_csv(step7_RawSim_summary_csv_location, index=False)
    # collecting core pvalues and core pvalue percentiles
    sim_loo_core_pvalues, sim_loo_core_pvalue_perc = collect_sim_core_pvalues_and_percentiles(small_and_large_sims,
                                                                                          loo_simulation_output)
    # save simulation core pvalues to df and export to csv
    sim_core_pvalue_df = pd.DataFrame(sim_loo_core_pvalues)
    sim_core_pvalue_df.to_csv(step7_SimCorePvalues_csv_location, index=False)
    # save simulation core pvalue percentiles to df and export to csv
    sim_percentile_df = pd.DataFrame(sim_loo_core_pvalue_perc)
    sim_percentile_df.to_csv(step7_SimCorePvaluePerc_csv_location, index=False)
    LogFile.write('Time to perform Step7: ' + str(time.time() - step7_start) + '\n')

    # Step8: Compare sample input/group data to simulated biological groups
    step8_start = time.time()
    print('Now Performing Step8: Comparison of sample input data to simulated biological groups')
    final_summary_matrix_updated = compare_to_simulations(small_and_large_sims, sim_percentile_df, final_summary_matrix)
    final_summary_matrix_updated.to_csv(step8_FinalSummary_location, index=False)
    summary_sig_core_genes = collect_core_genes_in_sig_groups(final_summary_matrix_updated)
    text_file = open(step8_SigCoreGenes_location, 'w')
    text_file.write('Core Genes in Pathways Above 5th Percentile: ' + '\n')
    text_file.write('Total genes: ' + str(len(summary_sig_core_genes)) + '\n')
    for gene in summary_sig_core_genes:
        text_file.write(gene + '\n')
    text_file.close()
    LogFile.write('Time to perform Step8: ' + str(time.time() - step8_start) + '\n')

    #Step9: Generate historgrams of significant biological groups
    # step9_start = time.time()
    # print('Now Performing Step9: Generating histograms of significant biological groups and core genes')
    # significant_groups_df = final_summary_matrix_updated.query("(passed_q_value_filter == 1) and (fold_better >= 1)")
    # print("Number significant pathways in final analysis: ", str(significant_groups_df.shape[0]))
    # significant_groups = list(significant_groups_df['group_name'])
    # for sig_group in significant_groups:
    #     full_biological_group_historgrams(sig_group, groups_input_df, sig_single_sample_genes_lst, sample_gene_EA_matrix, step9_HistCoreGenes_location)
    #     core_biological_group_historgrams(sig_group, significant_groups_df, sample_gene_EA_matrix, step9_HistCoreGenes_location)

    # LogFile.write('Time to perform Step9: ' + str(time.time() - step9_start) + '\n')

    # calculate time to complete full analysis
    print("Time to complete full analysis: ", time.time() - start_all)
    LogFile.write('Time to perform full analysis - ' + str(time.time() - start_all)  + '\n')
    LogFile.close()
