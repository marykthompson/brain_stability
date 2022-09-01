import pandas as pd
from collections import defaultdict
import numpy as np

cols =  ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
rows = ['A','B','C','D','E','F','G','H']

def prep_df(df, start):
    '''
    Read in the plate layout
    '''
    small_df =  df.iloc[start + 1: start + 10, 0:13].copy()
    #replace header with 1-12
    cols = ['row', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
    small_df.columns = cols
    small_df.drop(small_df.index[0], inplace = True)
    small_df.set_index('row', inplace = True)
    return small_df

def plate_to_series(plate_df, name = ''):
    '''
    Given a plate with rows A-H and cols 1-12, return series with rowname: val
    I think that it adds the underscore to the name because the name is also the positional value
    i.e. '1' = 1
    '''
    d = {}
    cols = range(1, 13)
    for i in plate_df.itertuples():
        rowname = i.Index
        for col in cols:
            wellname = '{}{}'.format(rowname, col)
            d[wellname] = i[col]
    data_s = pd.Series(d)
    data_s.name = name
    return data_s

def read_plate_data(data_file):
    if data_file.endswith('.csv'):
        data_df = pd.read_csv(data_file)
    elif data_file.endswith('.xlsx'):
        #data_df = pd.read_excel(data_file)
        data_df = pd.read_excel(data_file)
    else:
        print('cant read file!')

    new_index = list(zip(rows, ['Cq']*len(rows)))
    data_df = data_df.set_index(['Unnamed: 0', 'Unnamed: 1']).reindex(new_index)
    data_df.index.set_names(['row', 'Cq'], inplace = True)
    data_df.index = data_df.index.droplevel('Cq')
    return data_df

def main(data_file, template_file, ctrl_primer, drop_samples = []):

    template_df = pd.read_excel(template_file)

    #get index where primer plate starts
    primer_start = template_df.index[template_df['qPCR template sheet'] == 'PRIMERS'].tolist()[0]
    sample_start = template_df.index[template_df['qPCR template sheet'] == 'SAMPLES'].tolist()[0]

    norm_start = template_df.index[template_df['qPCR template sheet'] == 'NORMALIZATION'].tolist()[0]
    frac_start = template_df.index[template_df['qPCR template sheet'] == 'FRACTIONS'].tolist()[0]
    graphs_start = template_df.index[template_df['qPCR template sheet'] == 'TOGRAPH'].tolist()[0]

    primer_df = prep_df(template_df, primer_start)
    sample_df = prep_df(template_df, sample_start)

    data_df = read_plate_data(data_file)

    #Combine Cq, sample, primer values into well df
    cq_s = plate_to_series(data_df, name = 'Cq')
    sample_s = plate_to_series(sample_df, name = 'sample')
    primer_s = plate_to_series(primer_df, name = 'primer')

    #Split into control and experimental primers
    well_df = pd.concat([sample_s, primer_s, cq_s], axis = 1)
    ctrl_df = well_df[well_df['primer'] == ctrl_primer].copy().groupby('sample').mean()
    ctrl_df.columns = ['ctrl_Cq']
    exp_df = well_df[well_df['primer'] != ctrl_primer].copy().groupby(['sample', 'primer']).mean()
    dCt_df = pd.merge(exp_df, ctrl_df, left_index = True, right_index = True)

    #adjust the samples by their fractional value, i.e. rel to ng total RNA used at start of experiment
    #because frac_df is from a mixed df, the dtype is a string, convert to float
    frac_df = template_df.iloc[frac_start + 1: graphs_start-2, 0:2]
    frac_df.columns = ['sample', 'fraction']
    frac_df.set_index('sample', inplace = True)
    frac_df['fraction'] = pd.to_numeric(frac_df['fraction'])

    dCt_df = pd.merge(dCt_df, frac_df, left_index = True, right_index = True)

    dCt_df['adj_Cq'] = dCt_df['Cq'] - (1/dCt_df['fraction']).apply(np.log2)
    dCt_df['adj_ctrl_Cq'] = dCt_df['ctrl_Cq'] - (1/dCt_df['fraction']).apply(np.log2)

    dCt_df['dCt'] = dCt_df['adj_Cq'] - dCt_df['adj_ctrl_Cq']
    #Create the ddCt df by normalizing to the control sample
    norm_df = template_df.iloc[norm_start + 1: frac_start-2, 0:2]
    norm_df.columns = ['numerator', 'denominator']
    norm_df2 = norm_df.set_index('numerator')
    norm_df2.index.name = 'sample'

    #add the denominator info to the dCt
    numerator_df = pd.merge(dCt_df, norm_df2, left_index = True, right_index = True)
    numerator_df.set_index('denominator', append =  True, inplace = True)
    denominator_df = dCt_df[dCt_df.index.isin(norm_df['denominator'].values, level='sample')].copy()
    denominator_df.index.set_names('denominator', level = 'sample', inplace = True)
    #rename sample to denominator, then can merge it in
    ddCt_df = pd.merge(numerator_df, denominator_df, left_index = True, right_index = True, suffixes = ('_num', '_denom'))
    ddCt_df['ddCt'] = ddCt_df['dCt_num'] - ddCt_df['dCt_denom']
    ddCt_df['fold_change'] = 2**-ddCt_df['ddCt']
    result_df = ddCt_df[['ddCt', 'fold_change']].copy()
    #drop samples, this is for example to drop plate controls that we don't want to plot
    result_df = result_df[~result_df.index.get_level_values('sample').isin(drop_samples)].copy()
    return result_df

if __name__ == '__main__':
    main(sys.argv[1:])
