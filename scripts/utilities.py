'''
utilities.py
Functions to clean and filter data
'''
import pandas as pd
import os
from pyfaidx import Fasta

def filter_low_exp(count_file, filter_by='cpm', filter_co=1, npass=3,
                   read_col='summed_est_counts', experiments=None, outname=None):
    '''
    Remove genes from dataset which are not expressed. Filter by either CPM or count.
    filter_by = specify either cpm or count
    filter_co = numerical cutoff to pass filter, i.e. 1 for 1 CPM
    npass = number of samples which must pass this cutoff for gene to be included.
    experiments = list of dicts to specify experiment cutoffs. Use this to test
    whether at least one experiment meets the pass cutoff.
    e.g. experiments = [{'RNAtype':'input', 'condition':'mock', 'co':2},
    {'RNAtype':'pd', 'condition':1, 'co':2}]]
    '''

    df = pd.read_csv(count_file)
    df['CPM'] = df[read_col]*1e6/df.groupby(['replicate', 'condition', 'RNAtype'])[read_col].transform('sum')
    #Find genes which pass in at least one experiment
    filtered_genes = set()

    if experiments is not None:
        for i in experiments:
            RNAtype, condition, npass = i['RNAtype'], i['condition'], i['npass']
            s_df = df.loc[(df['RNAtype'] == RNAtype) & (df['condition'] == condition), ['gene', 'RNAtype', 'condition', 'replicate', 'CPM', read_col]]

            if filter_by == 'cpm':
                s_df['pass_exp_filter'] = s_df['CPM'] >= filter_co
            else:
                s_df['pass_exp_filter'] = s_df[read_col] >= filter_co

            pass_counts = s_df.groupby('gene')['pass_exp_filter'].sum()
            #Get genes which pass in at least npass libraries
            passed_genes = pass_counts[pass_counts >= npass].copy()
            #To get genes which pass in at least one experiment, simply add index to the set
            filtered_genes.update(passed_genes.index)
    else:
        if filter_by == 'cpm':
            df['pass_exp_filter'] = df['CPM'] >= filter_co
        else:
            df['pass_exp_filter'] = df[read_col] >= filter_co

        pass_counts = df.groupby('gene')['pass_exp_filter'].sum()
        passed_genes = pass_counts[pass_counts >= npass].copy()
        filtered_genes.update(passed_genes.index)

    print(f'{len(filtered_genes)} passed the filter')
    #If outname is provided, then write the genes passing the filter to a file.
    if outname is not None:
        pd.Series(list(filtered_genes)).to_csv(f'{outname}_passed.csv', index=False, header=None)
    return filtered_genes

def load_dataset(count_file, passed_gene_file):
    '''
    Load the dataset with genes not passing the expression cutoffs removed.
    '''
    df = pd.read_csv(count_file, index_col='gene')
    passed_genes = pd.read_csv(passed_gene_file, header=None)[0]
    df = df[df.index.isin(passed_genes)].copy()
    return df

def get_txt_id(s):
    '''
    This function is needed to extract the txt id for the CDS file from Flybase.
    '''
    s2 = s.split('; ')
    txt_id = list(filter(lambda x: x.startswith('parent'), s2))[0].split('=')[-1].split(',')[1]
    return txt_id

def parse_fb_fasta(infile, extract_ids=False):
    '''
    Use pyfaidx Fasta() to load the Flybase sequences. Use the FBtr IDs.
    '''
    if extract_ids:
        txt_seqs = Fasta(infile, read_long_names=True, key_function=get_txt_id)
    else:
        txt_seqs = Fasta(infile)
    return txt_seqs

