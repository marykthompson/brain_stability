#version 2, write functions to work from the individual replicate values
import numpy as np
import math
import pandas as pd
import scipy.stats as stats
from scipy.stats import hypergeom
import os
import gzip
import shutil
from urllib.request import urlretrieve
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt

from plot_helpers import *

class PrettyBox(object):
    '''
    Wrap the Seaborn boxplot in order to allow the box to have lower alpha.
    Also change the lines around the boxplot to black instead of grey.
    The Seaborn box and violin plot lines seem to remain grey despite specifying
    line colors in the rcParams.
    Seaborn boxplot also doesn't seem to respect alpha in colorcycle via RC params.
    '''
    def __new__(self, box_alpha=0.3, *args, **kwargs):
        ax = sns.boxplot(*args, **kwargs)
        plt.setp(ax.lines, color='k')
        boxes = ax.artists
        for box in boxes:
            r, g, b, a = box.get_facecolor()
            box.set_edgecolor('k')
            box.set_facecolor((r, g, b, box_alpha))
        if len(ax.get_legend_handles_labels()[0]) > 0:
            handles, labels = ax.get_legend_handles_labels()
            new_handles = []
            for lh in handles:
                lh.set_ec('k')
                lh.set_fc((*lh.get_fc()[0:3], box_alpha))
                lh.set_ec((0,0,0,1))
                new_handles.append(lh)
            #set new handles to change the linewidth and the alpha
            ax.legend(new_handles, labels, title=kwargs['hue'])
            # #to make the legend outline black, if desired
            # leg = ax.legend()
            # leg.get_frame().set_edgecolor('k')
            # leg.get_frame().set_lw(0.75)

        return ax

def plot_grey_scatter(data=None, x=None, y=None, genegroup=None, grouplabel=None):
    '''
    Plot scatter points in light grey and overlay members of a gene group.
    '''
    fig = plt.figure(figsize=(dfig, dfig), constrained_layout=True)
    ax = fig.add_subplot(111)
    x1 = data[x]
    y1 = data[y]
    ax.scatter(x1, y1, s=5, color='k', alpha=0.3, ec='none')

    if genegroup is not None:
        data['group'] = data.index.isin(genegroup)
        x2 = data[data['group']][x]
        y2 = data[data['group']][y]
        ax.scatter(x2, y2, s=5, color=color_dict['purple'], ec='none', label=grouplabel)

    rval, pval = stats.pearsonr(x1, y1)
    r2_val_av = rval**2
    loc = plticker.MultipleLocator(base=5.0)
    ax.xaxis.set_major_locator(loc)
    ax.yaxis.set_major_locator(loc)
    ax.text(0.1, 0.9, 'r'r'$^2$'' = %1.2f' % r2_val_av, fontsize = 8, transform=ax.transAxes)
    return ax

def get_boxtop(df, **kwargs):
    '''
    Get top of the boxplot whisker, in order to decide where to put the stars
    kwargs can be col1, val1, col2, val2 to denote the subsets in the dataframe
    val_col should be the value actually being ploted
    '''
    #check if two levels of nesting
    if 'col2' in kwargs:
        data = df.loc[(df[kwargs['col1']] == kwargs['val1']) & (df[kwargs['col2']] == kwargs['val2'])][kwargs['val_col']].values
    else:
        data = df.loc[df[kwargs['col1']] == kwargs['val1']][kwargs['val_col']].values

    q3, q1 = np.percentile(data, [75 ,25])
    iqr = q3 - q1
    whisk = q3 + (iqr*1.5)
    h = data[data < whisk].max()
    return h
#     if 'col2' in kwargs:
#         y1 = df.loc[(df[kwargs['col1']] == kwargs['val1']) & (df[kwargs['col2']] == kwargs['val2'])][kwargs['val_col']].quantile(0.25)
#         y2 = df.loc[(df[kwargs['col1']] == kwargs['val1']) & (df[kwargs['col2']] == kwargs['val2'])][kwargs['val_col']].quantile(0.75)
#     else:
#         y1 = df.loc[df[kwargs['col1']] == kwargs['val1']][kwargs['val_col']].quantile(0.25)
#         y2 = df.loc[df[kwargs['col1']] == kwargs['val1']][kwargs['val_col']].quantile(0.75)

#     iqr = y2-y1
#     return y2 + iqr*1.5

def add_stars(x1, x2, y, h, p, ax, siglevels=np.array([0.05, 10e-10, 10e-40]), col='k'):
    '''
    Add stars to the axis denoting significance level.
    Siglevels is the significance levels where p<level to give it a star
    '''
    sigstars = ['ns', '*', '**', '***']
    #this will return 0 if the value is larger than the first one
    j = np.digitize(p, siglevels)
    stars = sigstars[j]
    # j = np.digitize(p, siglevels) - 1
    # if j == -1:
    #     stars = 'ns'
    # else:
    #     stars = sigstars[j]
    #plot first one
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=0.75, c=col, alpha=0.5)
    ax.text((x1+x2)*.5, y+h, stars, ha='center', va='bottom', color=col, alpha=0.5)

def change_from_t0(df):
    '''
    Return a dataframe with the ratio between each timepoint
    and the first column (t0).
    '''
    cols = df.columns
    change_dict = {}
    for i in range(1, len(cols)):
        name = cols[i]
        change = df[cols[i]]/df[cols[0]]
        change_dict[name] = change
    new_df = pd.DataFrame(change_dict)
    return new_df

'''
Timecourse plotting functions
See discussion in NMJ_fig_planning for how to represent error on log scale
1/ln(2) * y_err/y = y_err/(ln(2)*y)
'''
#From F2, save this one
def extract_gene_vals(gene_id, val_df, var_df, tpts, label):
    '''If gene is not found return empty df.'''
    try:
        y = val_df.loc[gene_id]
        y_var = var_df.loc[gene_id]
        y_sd = np.sqrt(y_var)
        y_vals = []
        y_errs = []
        #if the value from the first point is divided by itself, then it's not independent
        #what is the correct way to calculate error in such a case?
        for i in range(0, len(y)):
            y_val = y[i]/y[0]
            y_err = y_val*np.sqrt((y_sd[i]/y[i])**2 + (y_sd[0]/y[0])**2)
            rel_err = y_err/(math.log(2)*y_val)
            if y_val != 0:
                y_vals.append(math.log(y_val, 2))
            else:
                y_vals.append(math.log(0.001, 2))
            y_errs.append(rel_err)

        df = pd.DataFrame({'vals':y_vals, 'err': y_errs, 'tpts':tpts, 'label': label})
        df.set_index('label', inplace = True)
    except KeyError:
        df = pd.DataFrame(columns=['vals', 'err', 'tpts', 'label'])
        df.set_index('label', inplace=True)
        return df
    return df

def extract_gene_vals_reps(gene_id, val_df, label):
    '''If gene is not found return empty df. This version take all the values from the replicates'''
    tpts = [float(i.split('_')[1]) for i in val_df.columns.values]

    try:
        y = val_df.loc[gene_id].mean()
        print('y', y)
        y_sd = val_df.loc[gene_id].std()
        print('ysd', y_sd)
        y_vals = []
        y_errs = []
        #if the value from the first point is divided by itself, then it's not independent
        #what is the correct way to calculate error in such a case?
        for i in range(0, len(y)):
            y_val = y[i]/y[0]
            y_err = y_val*np.sqrt((y_sd[i]/y[i])**2 + (y_sd[0]/y[0])**2)
            print('yerr', y_err)
            rel_err = y_err/(math.log(2)*y_val)
            print('relerr', rel_err)
            if y_val != 0:
                y_vals.append(math.log(y_val, 2))
            else:
                y_vals.append(math.log(0.001, 2))
            y_errs.append(rel_err)
        df = pd.DataFrame({'vals':y_vals, 'err': y_errs, 'tpts':tpts, 'label': label})
        df.set_index('label', inplace = True)
    except KeyError:
        df = pd.DataFrame(columns=['vals', 'err', 'tpts', 'label'])
        df.set_index('label', inplace=True)
        return df
    return df

'''
#From F1
#See discussion in NMJ_fig_planning for how to represent error on log scale
#1/ln(2) * y_err/y = y_err/(ln(2)*y)
def extract_gene_vals(gene_id, val_df, var_df):
    y = val_df.loc[gene_id]
    y_var = var_df.loc[gene_id]
    y_sd = np.sqrt(y_var)
    y_vals = []
    y_errs = []
    #if the value from the first point is divided by itself, then it's not independent
    #what is the correct way to calculate error in such a case?
    for i in range(0, len(y)):
        y_val = y[i]/y[0]
        y_err = y_val*np.sqrt((y_sd[i]/y[i])**2 + (y_sd[0]/y[0])**2)
        rel_err = y_err/(math.log(2)*y_val)
        y_vals.append(math.log(y_val, 2))
        y_errs.append(rel_err)

    return y_vals, y_errs
'''

#For F2
def plot_gene_tc(plot_df, title, ax = '', y_label = '', y_loc = None):
    '''
    Plot values over the timecourse. plot_df should be indexed by series name.
    Plot_df columns = tpts, vals, err.
    https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.errorbar.html
    Remove error bars from the legend:
    https://stackoverflow.com/questions/14297714/matplotlib-dont-show-errorbars-in-legend
    Set y-axis spacing with y_loc, 0.5 for Fig. 1
    '''
    for i in plot_df.index.unique():
        to_plot = plot_df.loc[i]
        ax.errorbar(to_plot['tpts'], to_plot['vals'], yerr = to_plot['err'], label = i)
        #ax.errorbar(to_plot['tpts'], to_plot['vals'], yerr = to_plot['err'], capsize = 5, capthick = 0.75, elinewidth = 0.75, label = i)
        #ax.errorbar(to_plot['tpts'], to_plot['vals'], yerr = to_plot['err'], capsize = 5, elinewidth = matplotlib.rcParams['lines.linewidth'], label = i)
    ax.set_ylabel(y_label)
    ax.set_xlabel('hrs after onset of activation')
    #'$\it{text you want to show in italics}$'
    ax.set_title(title)

    # get handles
    handles, labels = ax.get_legend_handles_labels()
    # remove the errorbars
    handles = [h[0] for h in handles]
    # use them in the legend
    #ax.legend(handles, labels, loc='upper left')
    #If not enough room to center the legend over the plot, then it looks better aligned to right
    #ax.legend(handles, labels, loc = 'lower center', bbox_to_anchor = (0.5, 1.0), handlelength = 0.75, ncol = 2, handletextpad = 0.5, columnspacing = 1)
    ax.legend(handles, labels, loc = 'lower right', bbox_to_anchor = (1.0, 1.03), ncol = 2, borderaxespad = 0)
    #borderaxespad is the arg that says how far away to put from axis. default is 0.5. Set to 0 to align with right
    ax.xaxis.set_major_locator(plticker.MultipleLocator(base=1.0))
    if y_loc is not None:
        ax.yaxis.set_major_locator(plticker.MultipleLocator(base=y_loc))
    return ax

def plot_gene_tc2(plot_df, title, ax = '', y_label = '', y_loc = None):
    '''
    Plot values over the timecourse. plot_df should be indexed by series name.
    Plot_df columns = tpts, vals, err.
    https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.errorbar.html
    Remove error bars from the legend:
    https://stackoverflow.com/questions/14297714/matplotlib-dont-show-errorbars-in-legend
    Set y-axis spacing with y_loc, 0.5 for Fig. 1
    '''
    for i in plot_df.index.unique():
        to_plot = plot_df.loc[i]
        ax.errorbar(to_plot['tpts'], to_plot['vals'], yerr = to_plot['err'], label = i)
        #ax.errorbar(to_plot['tpts'], to_plot['vals'], yerr = to_plot['err'], capsize = 5, capthick = 0.75, elinewidth = 0.75, label = i)
        #ax.errorbar(to_plot['tpts'], to_plot['vals'], yerr = to_plot['err'], capsize = 5, elinewidth = matplotlib.rcParams['lines.linewidth'], label = i)
    ax.set_ylabel(y_label)
    ax.set_xlabel('hrs after onset of activation')
    #'$\it{text you want to show in italics}$'
    ax.set_title(title)

    # get handles
    handles, labels = ax.get_legend_handles_labels()
    # remove the errorbars
    handles = [h[0] for h in handles]
    # use them in the legend
    #ax.legend(handles, labels, loc='upper left')
    #If not enough room to center the legend over the plot, then it looks better aligned to right
    #ax.legend(handles, labels, loc = 'lower center', bbox_to_anchor = (0.5, 1.0), handlelength = 0.75, ncol = 2, handletextpad = 0.5, columnspacing = 1)
    ax.legend(handles, labels, loc = 'lower right', bbox_to_anchor = (1.0, 1.03), ncol = 2, borderaxespad = 0)
    #borderaxespad is the arg that says how far away to put from axis. default is 0.5. Set to 0 to align with right
    ax.xaxis.set_major_locator(plticker.MultipleLocator(base=1.0))
    if y_loc is not None:
        ax.yaxis.set_major_locator(plticker.MultipleLocator(base=y_loc))
    return ax


def update_old_ids(genes, con_dict, as_set=False):
    '''
    Update ID set from old version of Flybase to the version in use.
    If we need to use a dataset with a newer version of Flybase, then it might make sense to update both to the current version
    genes = set of genes to convert
    con_dict = dictionary of old_ID: set of new IDs that can map to it
    if as_set, return genes as a set
    '''
    converted_genes = []
    for i in genes:
        for j in con_dict[i]:
            converted_genes.append(j)
    if as_set:
        converted_genes = set(converted_genes)
    return converted_genes

def update_old_ids_dict(genes, con_dict, discard_split=False, discard_merge=False):
    '''
    Update an ID set from old version of Flybase to the version in use.
    Note that it's possible to have both one-> many and many->one relationships
    Nest the results inside a list so that we can use pandas explode() to convert
    one->many mapping to multiple rows
    Return a dataframe mapping old_ID -> new_ID
    If discard_split == True, remove genes which now map to multiple new genes
    If discard_merge == True, remove genes which map to multiple old genes (does this class exist?)
    '''
    converted_genes = defaultdict(list)
    for i in genes:
        res = []
        for j in con_dict[i]:
            res.append(j)
        converted_genes[i] = [res]
    df = pd.DataFrame.from_dict(converted_genes, orient='index', columns=['new_ID'])
    df = df.explode('new_ID')
    n_split = len(df[df.index.duplicated(keep=False)])
    n_merge = len(df[df['new_ID'].duplicated(keep=False)].dropna())
    print('number of split genes = %s' % n_split)
    print('number of merge genes = %s' % n_merge)
    return df
    #return converted_genes

def update_old_ids_string(genes, con_dict):
    '''
    Update an ID set from old version of Flybase to the version in use.
    Note that it's possible to have both one-> many and many->one relationships
    For this version, make into a string list so that I can use pandas explode()
    '''
    converted_genes = defaultdict()
    for i in genes:
        res = []
        for j in con_dict[i]:
            res.append(j)
        converted_genes[i] = ','.join(res)
    return converted_genes

def calculate_overlap(set1, set2, bg_num):
    '''
    Test overlap of two gene lists with hypergeometric test
    M is the population size (previously N)
    n is the number of successes in the population (previously K)
    N is the sample size (previously n)
    X is still the number of drawn “successes”.
    checked some results here:
    https://systems.crump.ucla.edu/hypergeometric/index.php
    '''
    M = bg_num
    n = len(set1)
    N = len(set2)
    x = len(set1.intersection(set2))
    pval = hypergeom.sf(x-1, M, n, N)
    return pval

def retrieve_file(address, outfile):
    '''
    If url, download. If file, move it to the new location.
    Note the path of a file should be the absolute path because
    snakemake will likely be run in --directory mode.
    '''
    try:
        local_path, headers = urlretrieve(address, outfile)
    #not url
    except ValueError:
        shutil.copy(address, outfile)

def gunzip(infile):
    '''
    Unzip a file with .gz extension. Will remove extension in outfile.
    If the file does not have a .gz extension, do not unzip.
    '''
    if not infile.endswith('.gz'):
        return
    with gzip.open(infile, 'rb') as f_in:
        with open(infile.rstrip('.gz'), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(infile)
'''
#This is for F1, different enough that it makes sense to save separately
def plot_gene_tc(gene_id, gene_symbol, cluster_info, ax):
    #Remove error bars from the legend
    #https://stackoverflow.com/questions/14297714/matplotlib-dont-show-errorbars-in-legend
    syn_vals, syn_err = extract_gene_vals(gene_id, syn_rates, syn_var)
    tot_vals, tot_err = extract_gene_vals(gene_id, tot_levels, tot_var)
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(xvals, syn_vals, y_err = syn_err)
    #transcription = transcription rate
    #total RNA = total RNA level
    ax.errorbar(xvals, syn_vals, yerr = syn_err, capsize = 5, label = 'transcription')
    ax.errorbar(xvals, tot_vals, yerr = tot_err, capsize = 5, label = 'total RNA')
    #ax.plot(xvals, pd.Series(0).append(tot_vals))
    ax.set_ylabel('log'r'$_{2}$' ' fold change')
    ax.set_xlabel('hrs after onset of activation')
    #'$\it{text you want to show in italics}$'
    ax.set_title('$\it%s$ %s' % (gene_symbol, cluster_info))
    ax.xaxis.set_major_locator(plticker.MultipleLocator(base=1.0))
    allvals = syn_vals + tot_vals
    valrange = max(allvals) - min(allvals)
    if valrange < 1:
        ax.yaxis.set_major_locator(plticker.MultipleLocator(base=0.5))
    if valrange > 1:
        ax.yaxis.set_major_locator(plticker.MultipleLocator(base=1.0))
    return ax
'''
