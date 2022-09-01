##### Annotation utilities #####
#https://stackoverflow.com/questions/4172448/is-it-possible-to-break-a-long-line-to-multiple-lines-in-python
import pandas as pd
from collections import defaultdict

def parse_ann_file(ID_file):
    '''
    Take in a file of type fbgn_annotation_ID_fb_*_*.tsv from Flybase
    and write a conversion file of old_ID->current_ID.
    In the newer files, there is a column organism abbreviation.
    '''
    # df = pd.read_csv(ID_file, sep='\t', skiprows=5, header=None,
    # names=['gene_symbol', 'organism', 'primary_FBg', 'secondary_FBgs',
    #        'annotation_ID', 'secondary_annotation_IDs'])
    # df = df[df['organism'] == 'Dmel'].copy()
    #store old_id: {set of new IDs}
    df = pd.read_csv(ID_file, sep='\t', skiprows=4)
    cols = df.columns
    if 'organism_abbreviation' in df:
        df = df[df['organism_abbreviation'] == 'Dmel'].copy()

    con_dict = defaultdict(set)
    sym_dict = {}
    for index, i in df.iterrows():
        sym_dict[i['primary_FBgn#']] = i['##gene_symbol']
        if pd.isnull(i['secondary_FBgn#(s)']):
            con_dict[i['primary_FBgn#']] = set([i['primary_FBgn#']])
        else:
            sec_ids = i['secondary_FBgn#(s)'].split(',')
            for j in sec_ids:
                con_dict[j].add(i['primary_FBgn#'])
            #also add the current one mapping to the current one
            con_dict[i['primary_FBgn#']] = set([i['primary_FBgn#']])
    return con_dict, sym_dict

def update_ids(genes, to_anns, from_version='', id_type='FB'):
    '''
    Update an ID set from an old version of Flybase (from_anns) to a new version
    (to_anns). In the case of old ID split to mutliple new IDs, try to match on
    the gene symbol, which hopefully will not have also changed.
    from_anns and to_anns should be the fbgn_annotation_ID_fb_*_*.tsv file.
    First convert symbol to id if id_type == 'symbol'
    Also add option to not have from_anns, useful if it's not known.
    '''
    #ensure the genes are unique, otherwise they will appear as split annotations
    genes = set(genes)

    #get conversion from any old -> new and also current symbols
    con_dict, new_sym = parse_ann_file(to_anns)
    old_sym = None
    if from_version != '':
        #get the old symbols as well
        _, old_sym = parse_ann_file(from_version)

        if id_type == 'symbol':
            sym2id = {v:k for k,v in old_sym.items()}
            genes = [sym2id[i] if i in sym2id else i for i in genes]
    else:
        if id_type == 'symbol':
            raise Exception('symbol mapping cannot be used with unspecified from version')

    rows = []
    for i in genes:
        res = []
        for j in con_dict[i]:
            res.append(j)
        d = {}
        d[i] = [res]
        rows.append(pd.DataFrame.from_dict(d, orient='index'))
    df = pd.concat(rows)
    df.columns = ['new_ID']
    df = df.explode('new_ID')
    df['new_sym'] = df['new_ID'].map(new_sym)
    if old_sym is not None:
        df['old_sym'] = df.index.map(old_sym)
        #report ones not found, i.e. no old symbol.
        #Implies that genome version supplied by the user for the old IDs is incorrect
        n_notfound = len(df[pd.isnull(df['old_sym'])])
        print('number of genes not found = %s' % n_notfound)

    n_split = len(df[df.index.duplicated(keep=False)])
    n_merge = len(df[df['new_ID'].duplicated(keep=False)].dropna())

    #return split and merge df for checking
    # split_df = df[df.index.duplicated(keep=False)]
    # merge_df = df[df['new_ID'].duplicated(keep=False)].dropna()
    # return split_df, merge_df

    print('number of split genes = %s' % n_split)
    print('number of merge genes = %s' % n_merge)
    return df

def resolve_splits(df):
    '''
    Choose which of the new ids to preserve in old to new mapping.
    This is currently done by matching symbols.
    If no symbols are available (for example because the old version of the genome
    is unknown), then only return the unique ones.
    '''
    df.dropna(subset=['new_ID'], inplace=True)
    unique_df = df[~df.index.duplicated(keep=False)].copy()
    split_df = df[df.index.duplicated(keep=False)].copy()
    if not 'old_sym' in df.columns:
        return unique_df

    #if split_df is empty, nothing to do
    if split_df.empty:
        return df

    #Use startswith instead of ==. Sometimes a gene will be split into a and b.
    split_df['symbol_match'] = split_df.apply(lambda x: x['new_sym'].startswith(x['old_sym']),axis=1)

    accepted_df = split_df[split_df['symbol_match']].copy()
    accepted_df.drop(labels='symbol_match', axis=1, inplace=True)
    final_df = pd.concat([unique_df, accepted_df])
    return final_df

def resolve_merges(df):
    '''
    Choose which of the old IDs is represented by the new ID.
    Not yet implemented.
    '''
    pass
