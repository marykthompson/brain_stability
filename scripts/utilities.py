'''
utilities.py
Functions to clean and filter data
'''
import pandas as pd
from pyfaidx import Fasta
import numpy as np
import statsmodels.api as sm
from sklearn.preprocessing import PolynomialFeatures

def filter_low_exp(count_file, filter_co=1, npass=3,
                   filter_col='summed_est_counts', read_col='summed_est_counts', experiments=None, outname=None):
    '''
    Remove genes from dataset which are not expressed. Filter by either CPM or count.
    filter_col = specify the column name to filter on. Pass CPM to filter by calculated CPM value
    filter_co = numerical cutoff to pass filter, i.e. 1 for 1 CPM
    npass = number of samples which must pass this cutoff for gene to be included.
    experiments = list of dicts to specify experiment cutoffs. Use this to test
    whether at least one experiment meets the pass cutoff.
    e.g. experiments = [{'RNAtype':'input', 'condition':'mock', 'npass':2},
    {'RNAtype':'pd', 'condition':1, 'npass':2}]]
    '''

    df = pd.read_csv(count_file)
    df['CPM'] = df[read_col]*1e6/df.groupby(['replicate', 'condition', 'RNAtype'])[read_col].transform('sum')
    #Find genes which pass in at least one experiment
    filtered_genes = set()

    if experiments:
        for i in experiments:
            RNAtype, condition, npass = i['RNAtype'], i['condition'], i['npass']
            s_df = df.loc[(df['RNAtype'] == RNAtype) & (df['condition'] == condition), ['gene', 'RNAtype', 'condition', 'replicate', 'CPM', read_col]]
            s_df['pass_exp_filter'] = s_df[filter_col] >= filter_co
            # Get the number of libraries passed filter
            pass_counts = s_df.groupby('gene')['pass_exp_filter'].sum()
            # Get genes which pass in at least npass libraries
            passed_genes = pass_counts[pass_counts >= npass].copy()
            # To get genes which pass in at least one experiment, simply add index to the set
            filtered_genes.update(passed_genes.index)
    # if no experiments specified, perform filtering across all samples
    else:
        print('filtering all data')
        print('filter col', filter_col)
        print('npass', npass)
        df['pass_exp_filter'] = df[filter_col] >= filter_co
        pass_counts = df.groupby('gene')['pass_exp_filter'].sum()
        passed_genes = pass_counts[pass_counts >= npass].copy()
        filtered_genes.update(passed_genes.index)

    print(f'{len(filtered_genes)} passed the filter')
    #If outname is provided, then write the genes passing the filter to a file.
    if outname is not None:
        pd.Series(list(filtered_genes)).to_csv(f'{outname}_passed.csv', index=False, header=None)
    return filtered_genes

def calc_pseudocount_val(df, val_col=None, frac=0.1, groupby=None):
    '''
    Calculate a pseudocount value which is lower than all the real values
    '''
    # If a groupby list is passed, get mean of these groups
    if groupby:
        df = df.groupby(groupby).mean()
    a = np.unique(df[val_col])
    lowest_nonzero = a[a>0][0]
    return lowest_nonzero*frac

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
            
def feat_2_group(arr, nfeats):
    '''
    Convert the featurenames out to arrays specifying which group they belong to.
    nfeats is the number of input features before transformation.
    It is understood that 1 = 'intercept', not group 0,0,0
    However, I will assign it to group 0,0,0 for convenience, 
    and this will map Beta0.
    '''
    group_ids = []
    for i in arr:
        a = np.zeros(nfeats,)
        if i == '1':
            pass
        else:
            l = [int(j.lstrip('x')) for j in i.split(' ')]
            a[l] = 1
        group_ids.append(a)
    return group_ids

def combine_predictors(group_ids, predictors, beta_df=None):
    '''
    Given a list of arrays (group_ids), return either coefficient labels 
    (beta_df given) or a column of combined predictor names in correct order.
    '''
    out = pd.DataFrame(group_ids).astype(int)
    if beta_df is not None:
        out2 = pd.merge(out, beta_df, left_on=list(range(0, len(predictors))), right_on=list(range(0, len(predictors))))
        labels = out2['beta_label']
    else:
        labels = []
        pred_arr = np.array(predictors)
        for i in group_ids:
            l = np.compress(i, pred_arr)
            labels.append(':'.join(l)) 
    return labels

def run_OLS(df, target_column, predictors, interactions=False, beta_df=None):
    '''
    Run linear regression on the data given a target column and predictors.
    Use statsmodels.ols and pass in a dataframe with either the actual predictor names
    or a betas, which need to be specified via another mapping df.
    beta_df = df with integer labelled columns corresponding to the predictors (i.e. [0,1,2]),
    1/0 in each row and a column called 'beta_label'.
    '''
    X = df[predictors].values
    y = df[target_column].values
    if interactions:
        # l = combine_pred_names(predictors)
        poly = PolynomialFeatures(interaction_only=True, degree=len(predictors))
        # When using PolynomialFeatures, it already adds a constant.
        X_tr = poly.fit_transform(X)
        # Get group ids, this is in the from of arrays specifiying 1/0 for ech predictor
        group_ids = feat_2_group(poly.get_feature_names_out(), len(predictors))
        l = combine_predictors(group_ids, predictors, beta_df=beta_df)
        Xt = pd.DataFrame(X_tr, columns=l)
        mod = sm.OLS(y, Xt)
        res = mod.fit()
        res.summary()
    else:
        X1 = sm.add_constant(X.astype(int))
        if beta_df is not None:
            columns = beta_df['beta_label'][0:len(predictors)+1]
        else:
            columns = ['intercept'] + predictors
        X1 = pd.DataFrame(X1, columns=columns)
        mod = sm.OLS(y, X1)
        res = mod.fit()
        res.summary()
    return res

