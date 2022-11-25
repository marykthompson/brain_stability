##### Annotation utilities #####
import pandas as pd
from collections import defaultdict

def parse_ann_file(ID_file):
    '''
    Take in a file of type fbgn_annotation_ID_fb_*_*.tsv from Flybase
    and write a conversion file of old_ID->current_ID.
    In the newer files, there is a column organism abbreviation.
    In the older files, there is no column for organism abbreviation but
    all the symbols start with <Species>\<symbol> where species is
    something like Dmoj. Remove lines with \ in the symbol
    '''
    df = pd.read_csv(ID_file, sep='\t', skiprows=4, keep_default_na=False, na_values=[''])
    cols = df.columns
    if 'organism_abbreviation' in df:
        df = df[df['organism_abbreviation'] == 'Dmel'].copy()
    # Old-style file, remove non Dmel genes
    else:
        # Remove genes with a backslash in the gene symbol, which denotes the species
        import re
        m = re.compile('\\\\')
        df['nonDmel'] = df['##gene_symbol'].apply(lambda x: True if m.search(x) else False)
        df = df.query('~nonDmel').copy()

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

def update_ids(to_anns, from_version='', id_type='FB', genes=None):
    '''
    geneset is an iterable of genes to convert. Otherwise convert all the old IDs
    present in the annotation file to new IDs
    Update an ID set from an old version of Flybase (from_anns) to a new version
    (to_anns). In the case of old ID split to mutliple new IDs, try to match on
    the gene symbol, which hopefully will not have also changed.
    from_anns and to_anns should be the fbgn_annotation_ID_fb_*_*.tsv file.
    First convert symbol to id if id_type == 'symbol'
    Also add option to not have from_anns, useful if it's not known.
    '''
    # In order for this to work properly, the to_anns MUST be more recent than the from anns
    int_v = lambda x: '_'.join(x.strip('.tsv').split('_')[-2:])
    assert int_v(to_anns) >= int_v(from_version)

    # get conversion from any old -> new and also current symbols
    con_dict, new_sym = parse_ann_file(to_anns)
    old_sym = None

    print('to version', to_anns)
    print('from version', from_version)

    if from_version != '':
        # get the old symbols as well
        old_ids, old_sym = parse_ann_file(from_version)

        if id_type == 'symbol':
            sym2id = {v:k for k,v in old_sym.items()}
            genes = [sym2id[i] if i in sym2id else i for i in genes]
        # If no gene group specified for conversion, populate with all genes in the file
        if not genes:
            # Get all the values of the old_ids dictionary
            genes = set().union(*old_ids.values())
    else:
        if id_type == 'symbol':
            raise Exception('symbol mapping cannot be used with unspecified from version')

    # ensure the genes are unique, otherwise they will appear as split annotations
    genes = set(genes)

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
        # report ones not found, i.e. no old symbol.
        # Implies that genome version supplied by the user for the old IDs is incorrect
        n_notfound = len(df[pd.isnull(df['old_sym'])])
        print('number of gene symbols not found in old version = %s' % n_notfound)

    n_split = len(df[df.index.duplicated(keep=False)].dropna(subset=['new_ID']))
    n_merge = len(df[df['new_ID'].duplicated(keep=False)].dropna(subset=['new_ID']))

    # return split and merge df for checking
    # split_df = df[df.index.duplicated(keep=False)]
    # merge_df = df[df['new_ID'].duplicated(keep=False)].dropna()
    # return split_df, merge_df

    print('number of split genes = %s' % n_split)
    print('number of merge genes = %s' % n_merge)
    final_df = resolve_splits(df)
    return final_df

def common_start(s):
    '''
    Return the longest that matches the start of all strings in a series
    https://codereview.stackexchange.com/questions/237340/find-common-substring-that-starts-a-set-of-strings
    '''
    s = s.values
    l = len(s)
    if l == 0:
        return None
    elif l == 1:
        return s[0]
    # This goes through the list and checks if all start with the same string
    # Otherwise remove the last letter and test again
    start = s[0]
    while start != "" and not all(ii.startswith(start) for ii in s[1:]):
        start = start[:-1]
    return start

def resolve_splits(df, old_sym='old_sym', new_sym='new_sym', new_ID='new_ID', common_start_co=3):
    '''
    Choose which of the new ids to preserve in old to new mapping.
    This is currently done by matching symbols.
    If no symbols are available (for example because the old version of the genome
    is unknown), then only return the unique ones.
    common_start_co = the length of a common start prefix to be considered from the same gene
    for the case of matching on non-exact symbols
    '''
    # Remove genes which have no newID. These could be due to withdrawn gene models.
    df = df.dropna(subset=[new_ID]).copy()
    unique_df = df[~df.index.duplicated(keep=False)].copy()
    split_df = df[df.index.duplicated(keep=False)].copy()
    if not old_sym in df.columns:
        return unique_df

    #if split_df is empty, nothing to do
    if split_df.empty:
        return df
    else:
        print('not empty')
    
    # Use startswith instead of ==. Sometimes a gene will be split into a and b.
    # Only use if the symbol is >= common_start_co
    # Make it case insensitive
    split_df['symbol_match'] = split_df.apply(lambda x: x[new_sym].upper().startswith(x[old_sym].upper()),axis=1)
    split_df['symbol_match2'] = split_df['symbol_match'] & split_df[old_sym].apply(lambda x: len(x) > common_start_co)

    # How to check if both values are renamed in a similar way, which suggests a nomenclature change affecting common genes
    # Get the length of the longest common starting string
    split_df['len_common_start'] = split_df.groupby(old_sym)[new_sym].apply(lambda x: len(common_start(x)))
    split_df['common_start'] = split_df['len_common_start'] >= common_start_co
    accepted_df = split_df[(split_df['common_start'] | split_df['symbol_match2'])].copy()
    accepted_df.drop(labels=['symbol_match', 'symbol_match2', 'len_common_start', 'common_start'], axis=1, inplace=True)
    final_df = pd.concat([unique_df, accepted_df])
    return final_df

def txt_2_gene_dict(db):
    '''
    Return dict of transcript -> gene id
    '''
    txt_2_gene = {}
    genes = db.features_of_type('gene')
    for i in genes:
        l = [j for j in db.children(i, featuretype='mRNA')]
        for txt in l:
            txt_id = txt.attributes['transcript_id'][0]
            txt_2_gene[txt_id] = i.id
    return txt_2_gene