
# coding: utf-8

# In[ ]:



# Author: Vanessa Brisson (vbrisson@berkeley.edu)

# HabitatCorrectedNetwork.py is a python module for calculating correlation matrices for 
# network analysis that account for habitat filtering effects.  

import os
import numpy as np
from scipy import stats
from statsmodels.sandbox.stats import multicomp
import ConfigParser
import argparse

def read_OTU_table(filename):
    # Read in the OTU (ASV) abundance data file
    # Input OTU (ASV) table format:
    #   first row lists OTU (ASV) IDs
    #   subsequent rows start with the sample name followed by OTU counts for that sample
    print('Reading OTU table data')
    sample_names = []
    f = open(filename, 'r')
    first = True
    for line in f:
        columns = line.split()
        if first:
            master_OTU_list = columns
            OTU_table = np.zeros((0,len(master_OTU_list)))
            first = False
        else:
            sample_names.append(columns[0])
            abundances = map(float, columns[1:len(columns)])
            OTU_table = np.vstack((OTU_table, abundances))
    f.close()
    return(OTU_table, master_OTU_list, sample_names)

def read_sample_data(filename):
    # Read in the sample metadata
    # Input sample data format:
    #   first row gives column names for sample variables (must include sample_name)
    #   subsequent rows give data for each sample
    print('Reading sample data')
    subgroup_data = {}
    f = open(filename, 'r')
    first = True
    for line in f:
        columns = line.split()
        if first:
            subgroup_names = columns
            for sgn in subgroup_names:
                subgroup_data[sgn] = []
            first = False
        else:
            for idx in range(len(subgroup_names)):
                sgn = subgroup_names[idx]
                subgroup_data[sgn].append(columns[idx])
    f.close()
    return(subgroup_data)

def group_data(full_data_array, sample_data, split_group, sub_group):
    # Group the data for analysis
    #   Split the data into separate analyses by split_group
    #   Identify sample sub groups (habitats), based on sub_group
    print('Grouping data into separate networks by variable: ' + split_group)
    print('Performing habitat filtering correction by habitat variable: ' + sub_group)
    sample_names = sample_data['sample_name']
    if sub_group == 'none':
        sample_subgroup_info = np.full(len(sample_names), 'none').tolist()
    else:
        sample_subgroup_info = sample_data[sub_group]
    num_OTUs = full_data_array.shape[1]
    if split_group == 'none':
        split_data_arrays = {'all_data' : full_data_array}
        sample_subsets = {'all_data' : sample_subgroup_info}
    else:
        sample_split_group_info = sample_data[split_group]
        split_group_list = list(set(sample_split_group_info))
        split_data_arrays = dict([(name, np.zeros((0, num_OTUs))) for name in split_group_list])
        sample_subsets = dict([(name, []) for name in split_group_list])
        for idx in range(len(sample_names)):
            sg = sample_split_group_info[idx]
            split_data_arrays[sg] = np.vstack((split_data_arrays[sg], full_data_array[idx, :]))
            ss = sample_subgroup_info[idx]
            sample_subsets[sg].append(ss)
    return(split_data_arrays, sample_subsets)

def filter_OTU_prevalence(data_array, master_OTU_list, min_prev = 0.5):
    # Filter out low prevalence OTUs (ASVs) 
    #   min_prev is the minimum fraction of samples the OTU (ASV) must be present in to 
    #   be included in the analysis.
    print('Filtering out OTUs below minimum prevalence level of ' + str(min_prev))
    low_prevalence = []
    for OTU_idx in range(data_array.shape[1]):
        OTU_data = data_array[:,OTU_idx]
        if sum(OTU_data > 0)/float(len(OTU_data)) < min_prev:
            low_prevalence.append(OTU_idx)
    result = np.delete(data_array, low_prevalence, 1)
    OTU_list = (np.delete(np.array(master_OTU_list), low_prevalence)).tolist()
    return(result, OTU_list)

def clr_transform_data(data_array):
    # Centered Log Ratio Transform to compensate for compositionality effects
    print('Performing CLR transform')
    # clr(x) = ln(x/g(x)) = ln(x) - mean(ln(x))
    pseudo = min([0.5, np.min(0.5 * data_array[data_array > 0])])
    data_array_no_zeros = data_array + pseudo
    result = np.zeros(data_array_no_zeros.shape)
    gmeans = stats.gmean(data_array_no_zeros, 1)
    for sample_idx in range(data_array_no_zeros.shape[0]):
        for OTU_idx in range(data_array_no_zeros.shape[1]):
            result[sample_idx, OTU_idx] = np.log(data_array_no_zeros[sample_idx, OTU_idx] / gmeans[sample_idx])
    return(result)

def subtract_subgroup_means(data_array, sample_subgroups):
    # Perform habitat filterin correction by subtracting out within habitat (sub_group) mean 
    # abundance for each OTU (ASV) in each sample
    print('Performing habitat filtering correction')
    result = np.zeros(data_array.shape)
    subgroups = list(set(sample_subgroups))
    for OTU_idx in range(data_array.shape[1]):
        sgmeans = {}
        for sg in subgroups:
            sgmeans[sg] = np.mean(data_array[:,OTU_idx][np.array(sample_subgroups) == sg])
        for sample_idx in range(data_array.shape[0]):
            sg = sample_subgroups[sample_idx]
            result[sample_idx, OTU_idx] = data_array[sample_idx, OTU_idx] - sgmeans[sg]
    return(result)

def calculate_pearson_correlations(data_array):
    # Perform Pearson correlations
    print('Calculating Pearson correlatios')
    N_OTUs = data_array.shape[1]
    corr = np.zeros((N_OTUs, N_OTUs))
    pvals = np.zeros((N_OTUs, N_OTUs))
    for OTU_idx1 in range(N_OTUs):
        for OTU_idx2 in range(N_OTUs):
            if OTU_idx1 == OTU_idx2:
                r = 0
                p = 1
            else:
                r, p = stats.pearsonr(data_array[:,OTU_idx1], data_array[:, OTU_idx2])
            corr[OTU_idx1, OTU_idx2] = r
            pvals[OTU_idx1, OTU_idx2] = p
    return(corr, pvals)

def calculate_spearman_correlations(data_array):
    # Perform Spearman corelations
    print('Calculating Spearman correlations')
    N_OTUs = data_array.shape[1]
    corr, pvals = stats.spearmanr(data_array)
    for idx in range(N_OTUs):
        corr[idx,idx] = 0
        pvals[idx,idx] = 1
    return(corr, pvals)

def adjust_pvalues(pvals, method):
    # Adjust p-values for multiple comparisons 
    print('Adjusting p-values for multiple correlations')
    dim = pvals.shape[0]
    indices = np.triu_indices(dim, k=1)
    adjusted = multicomp.multipletests(pvals[indices], method = method)[1]
    new_pvals = np.ones((dim,dim))
    for i in range(len(adjusted)):
        idx1 = indices[0][i]
        idx2 = indices[1][i]
        new_pvals[idx1][idx2] = adjusted[i]
        new_pvals[idx2][idx1] = adjusted[i]
    return(new_pvals)

def threshold_significance(corr, pvals, min_corr, max_pval, direction, pval_correction):
    # Retain only correlations that meet selected ctiteria:
    #   p-value (adjusted if pval_correction = True) is below threshold set by max_pval
    #   correlation coefficient is above threshold set by min_cor
    #   correlation direction matches selected direction
    print('Selecting only correlations that meet the folowing criteria:')
    if pval_correction == True:
        print('   (FDR adjusted) pvalue <= ' + str(max_pval))
    else:
        print('   pvalue <= ' + str(max_pval))
    print('   correlation coeficient  >= ' + str(min_corr))
    if direction in ['pos', 'positive']:
        print('   positive correlations only')
    if direction in ['neg', 'negative']:
        print('   negative correlations only')
    significant = np.where(pvals <= max_pval, corr, 0)
    if direction == 'pos' or direction == 'positive':
        significant = np.where(significant >= min_corr, significant, 0) # POSITIVE CORRELATIONS ONLY
    elif direction == 'neg' or direction == 'negative':
        significant = np.where(significant <= 0 - min_corr, significant, 0) # NEGATIVE CORRELATIONS ONLY
    elif direction == 'both':
        significant = np.where(np.absolute(significant) >= min_corr, significant, 0)
    return(significant)

def valid_entry(x):
    # Check that entered value is a number etween 0 and 1
    try:
        float(x)
        if float(x) >= 0 and float(x) <= 1:
            return True
        else:
            return False
    except ValueError:
        return False

def custom_split(valid_groups):
    # Allow user to interactibely enter a split_group for dividing data into multiple networks
    split_group = 'none'
    split_YN = raw_input('Do you want to split the data into multiple networks (yes/no)? ')
    while not(split_YN in ['y','Y','yes','n','N', 'no']):
        split_YN = raw_input('Please answer "yes" or "no". ')
    if split_YN in ['y', 'Y', 'yes']:
        print('Enter the variable you wish to use to split the data. ')
        print('Valid choices are: ')
        print(valid_groups)
        split_group = raw_input()
        while not(split_group in valid_groups and split_group != 'none'):
            print('Please enter a valid variable name from this list: ')
            print(valid_groups)
            print('Or enter "none" if you do not wish to split the data.')
            split_group = raw_input()
    return(split_group)

def custom_HF_groups(valid_groups):
    # Allow user to interactiely enter a sub_group (habitat group) for HF correction
    sub_group = 'none'
    HF_YN = raw_input('Do you want to perform habitat filtering correction (yes/no)? ')
    while not(HF_YN in ['y','Y','yes','n','N', 'no']):
        HF_YN = raw_input('Please answer "yes" or "no". ')
    if HF_YN in ['y', 'Y', 'yes']:    
        print('Please enter the habitat filtering variable you wish to correct for. ')
        print('Valid choices are: ')
        print(valid_groups)
        sub_group = raw_input()
        while not(sub_group in valid_groups and sub_group != 'none'):
            print('Please enter a valid variable name from this list: ')
            print(valid_groups)
            sub_group = raw_input()
    return(sub_group)

def custom_param():  
    # Allow user to interactivly enter analysis parameters
    print()
    CLR = False
    CLR_YN = raw_input('Do you want to CLR transform the data to correct for compositional effects (yes/no)? ')
    while not(CLR_YN in ['y','Y','yes','n','N', 'no']):
        CLR_YN = raw_input('Please answer "yes" or "no". ')
    if CLR_YN in ['y', 'Y', 'yes']:
        CLR = True
    print()     
    pval_correction = False
    pval_YN = raw_input('Do you want to adjust the pvalues for multipe comparisons (yes/no)? ')
    while not(pval_YN in ['y','Y','yes','n','N', 'no']):
        pval_YN = raw_input('Please answer "yes" or "no". ')
    if pval_YN in ['y', 'Y', 'yes']:
        pval_correction = True
    print()    
    corr_method = raw_input('Select the desired correlation method (Spearman or Pearson): ')
    while not(corr_method in ['spearman', 'Spearman', 'pearson', 'Pearson']):
        corr_method = raw_input('Please enter either "Spearman" or "Pearson". ')
    print()    
    min_prev = raw_input('Please enter the minimum OTU prevalence (between 0 and 1). ')
    while not(valid_entry(min_prev)):
        min_prev = raw_input('Please enter a number between 0 and 1. ')
    min_prev = float(min_prev)
    print()
    min_corr = raw_input('Please enter the minimum correlation value (between 0 and 1). ')
    while not(valid_entry(min_corr)):
        min_corr = raw_input('Please enter a number between 0 and 1. ')
    min_corr = float(min_corr)
    print()
    max_pval = raw_input('Please enter the maximum p-value (between 0 and 1). ')
    while not(valid_entry(max_pval)):
        max_pval = raw_input('Please enter a number between 0 and 1. ')
    max_pval = float(max_pval)
    print()
    direction = 'pos'
    direction = raw_input('Which type(s) of correlation do you want to include (positive, negative, or both)? ')
    while not(direction in ['pos','positive','neg', 'negative', 'both']):
        direction = raw_input('Please enter "positive", "negative", or "both". ')
    print()    
    return(CLR, pval_correction, corr_method, min_prev, min_corr, max_pval, direction)

def check_parameters(CLR, pval_correction, corr_method, min_prev, min_corr, max_pval, direction):
    # Check that entered parameters are valid
    if not(type(CLR) == bool):
        print('ERROR: invalid entry.  "CLR" must be either "True" or "False"')
        return False
    if pval_correction == "True":
        pval_correction = True
    elif pval_correction == "False":
        pval_correction = False
    if not(type(pval_correction) == bool):
        print('ERROR: invalid entry.  "pval_correction" must be either "True" or "False"')
        return False
    if not(corr_method in ['spearman', 'Spearman', 'pearson', 'Pearson']):
        print('ERROR: invalid entry.  "corr_method" must be either "Spearman" or "Pearson"')
        return False
    if not(valid_entry(min_prev)):
        print('ERROR: invalid entry.  "min_prev" must be a number between 0 and 1')
        return False
    if not(valid_entry(min_corr)):
        print('ERROR: invalid entry.  "min_corr" must be a number between 0 and 1')
        return False
    if not(valid_entry(max_pval)):
        print('ERROR: invalid entry.  "max_pval" must be a number between 0 and 1')
        return False
    if not(direction in ['pos','positive','neg', 'negative', 'both']):
        print('ERROR: invalid entry.  "direction" must be "positive", "negative", or "both"')
        return False
    return True

def HF_analysis(OTU_filename, 
                sample_filename,
                setup = 'default', 
                setup_filename = '',
                split_group = 'none',
                sub_group = 'none',
                CLR = True,
                pval_correction = False,
                corr_method = 'spearman',
                min_prev = 0.66,
                min_corr = 0,
                max_pval = 0.05,
                direction = 'pos'):
    
    # Conduct Habitat Filtering Corrected correlation analysis using the entered parameters
    
    pval_method = 'fdr_bh'
    full_data_array, master_OTU_list, sample_names = read_OTU_table(OTU_filename)
    sample_data = read_sample_data(sample_filename)

    if setup == 'custom':
        split_group = custom_split(sample_data.keys())
        sub_group = custom_HF_groups(sample_data.keys())
        CLR, pval_correction, corr_method, min_prev, min_corr, max_pval, direction = custom_param()
    elif setup == 'config':
        if os.path.isfile(setup_filename):
            config = ConfigParser.RawConfigParser()
            try:
                config.read(setup_filename)
                if not(config.sections() == ['NetworkSplit', 'HabitatFiltering', 'Parameters']):
                    print('ERROR: ' + setup_filename + 
                          ' is not a valid config file for habitat filtering correction.')
                    return
                split_group = config.get('NetworkSplit', 'split_group')
                sub_group = config.get('HabitatFiltering', 'habitat_group')
                CLR = config.getboolean('Parameters', 'CLR')
                pval_correction = config.getboolean('Parameters', 'pval_correction')
                corr_method = config.get('Parameters', 'corr_method')
                min_prev = config.getfloat('Parameters', 'min_prev')
                min_corr = config.getfloat('Parameters', 'min_corr')
                max_pval = config.getfloat('Parameters', 'max_pval')
                direction = config.get('Parameters', 'direction')
            except ConfigParser.MissingSectionHeaderError:
                print('ERROR: ' + setup_filename + ' is not a valid config file.')
                return 
        else:
            print('ERROR: You must enter a valid file name for "setup_filename"')
            return
    
    if not(split_group == 'none' or split_group in (sample_data.keys())):
        print('ERROR: split_group must be "none" or correspond to a variable in the sample data.')
        print('Valid groups for these data:')
        print(sample_data.keys())
        return
    if not(sub_group == 'none' or sub_group in (sample_data.keys())):
        print('ERROR: sub_group must be "none" or correspond to a variable in the sample data.')
        print('Valid groups for these data:')
        print(sample_data.keys())
        return
    
    split_data_arrays, sample_subsets = group_data(full_data_array, sample_data, split_group, sub_group)
    split_data_arrays_filt = {}
    correlations = {}
    pvalues = {}
    OTU_subsets = {}
    
    for group in split_data_arrays.keys(): 
        print()
        print('Processing data:\t' + group)
        split_data_arrays_filt[group], OTU_subsets[group] = filter_OTU_prevalence(split_data_arrays[group], 
                                                                             master_OTU_list, 
                                                                             min_prev)
        my_data = split_data_arrays_filt[group].copy()
        if CLR == True: 
            my_data = clr_transform_data(my_data)
        my_data = subtract_subgroup_means(my_data, sample_subsets[group])
        if corr_method == 'spearman' or corr_method == 'Spearman':
            c, p = calculate_spearman_correlations(my_data)
        elif corr_method == 'pearson' or corr_method == 'Pearson':
            c, p = calculate_pearson_correlations(my_data)
        if pval_correction == True:
            p = adjust_pvalues(p, pval_method)
        c = threshold_significance(c, p, min_corr, max_pval, direction, pval_correction)
        correlations[group] = c
        pvalues[group] = p
    return (correlations, pvalues, OTU_subsets)

def save_output(correlations,
                pvalues,
                OTU_subsets,
                output_prefix):
    # Save correlation matrix, p-value matrix, and OTU (ASV) generated in the analysis
    for k in correlations.keys():
        np.savetxt(output_prefix + "_" + k + '_correlations.txt', correlations[k])
        np.savetxt(output_prefix + "_" + k + '_pvalues.txt', pvalues[k])
        f = open(output_prefix + "_" + k + '_ASVs.txt', 'w')
        for OTU in OTU_subsets[k]:
            f.write(OTU + '\n')
        f.close()
    return

def parse_input():
    #  Parse command line input
    OTU_filename = ''
    sample_filename = ''
    setup = 'default'
    setup_filename = ''
    split_group = 'none'
    sub_group = 'none'
    CLR = True
    pval_correction = False
    corr_method = 'spearman'
    min_prev = 0.66
    min_corr = 0
    max_pval = 0.05
    direction = 'pos'
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--ASV_filename', '-A', required=True)
    parser.add_argument('--sample_filename', '-S', required=True)
    parser.add_argument('--output_prefix', '-out', required=True)
    parser.add_argument('--setup', '-s')
    parser.add_argument('--setup_filename', '-sfn')
    parser.add_argument('--split_group', '-split')
    parser.add_argument('--habitat_group', '-hab')
    parser.add_argument('--CLR','-clr')
    parser.add_argument('--pval_adjustment','-padj')
    parser.add_argument('--corr_method', '-cm')
    parser.add_argument('--min_prev', '-prev', type=float)
    parser.add_argument('--min_corr', '-cor', type=float)
    parser.add_argument('--max_pval', '-pval', type=float)
    parser.add_argument('--direction', '-d')
    args = parser.parse_args()
    
    OTU_filename = args.ASV_filename
    sample_filename = args.sample_filename
    output_prefix = args.output_prefix
    if args.setup is not None:
        setup = args.setup
    if args.setup == 'config':
        if args.setup_filename is None:
            print('ERROR: You must enter a setup_filename to define parameters with a config file.')
            return
        else:
            setup_filename = args.setup_filename
    if args.split_group is not None:
        split_group = args.split_group
    if args.habitat_group is not None:
        sub_group = args.habitat_group
    if args.CLR is not None:
        CLR = args.CLR
        if CLR == "True":
            CLR = True
        elif CLR == "False":
            CLR = False
    if args.pval_adjustment is not None:
        pval_correction = args.pval_adjustment
        if pval_correction == "True":
            pval_correction = True
        elif pval_correction == "False":
            pval_correction = False
    if args.corr_method is not None:
        corr_method = args.corr_method
    if args.min_prev is not None:
        min_prev = args.min_prev
    if args.min_corr is not None:
        min_corr = args.min_corr
    if args.max_pval is not None:
        max_pval = args.max_pval
    if args.direction is not None:
        direction = args.direction
        
    return(OTU_filename,
           sample_filename,
           output_prefix,
           setup,
           setup_filename,
           split_group,
           sub_group,
           CLR, 
           pval_correction, 
           corr_method, 
           min_prev, 
           min_corr, 
           max_pval, 
           direction)

def main():
    (OTU_filename,
    sample_filename,
    output_prefix,
    setup,
    setup_filename,
    split_group,
    sub_group,
    CLR, 
    pval_correction, 
    corr_method, 
    min_prev, 
    min_corr, 
    max_pval, 
    direction) = parse_input()

    if check_parameters(CLR, pval_correction, corr_method, min_prev, min_corr, max_pval, direction):    
        (correlations, 
         pvalues, 
         OTU_subsets) = HF_analysis(OTU_filename, 
                                    sample_filename,
                                    setup, 
                                    setup_filename,
                                    split_group,
                                    sub_group,
                                    CLR,
                                    pval_correction,
                                    corr_method,
                                    min_prev,
                                    min_corr,
                                    max_pval,
                                    direction)
        save_output(correlations,
                    pvalues,
                    OTU_subsets,
                    output_prefix)

if __name__ == "__main__":
	main()

