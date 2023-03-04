'''
More complex plotting functions
'''
import numpy as np
import pandas as pd
import scipy.stats as stats
import os
import gzip
import shutil
from urllib.request import urlretrieve
import matplotlib as mpl
from matplotlib import lines
import seaborn as sns
import matplotlib.pyplot as plt
from stats_helpers import calc_fisher_exact
from decimal import Decimal
import matplotlib.colors as colors

from plot_helpers import *

def get_box_coords(ax):
    '''
    Return a list of the middle x-coords, in data coordinates, of boxes. 
    Assumes that all the matplotlib.patches.PathPatch objects are boxes.
    '''
    boxes = ax.patches
    display_to_data = ax.transData.inverted()
    x_pos = []
    for box in boxes:
        if type(box) == mpl.patches.PathPatch:
            x_middle = (box.get_extents().x0 + box.get_extents().x1)/2
            # transform needs 2D arg
            x, _ = display_to_data.transform((x_middle, 0))
            x_pos.append(x)
    return x_pos

def get_letter_height(ax, fig, fontsize=7):
    '''
    Return the height of the letter "R" in data coordinates.
    This can be used to place text and compensate for issues like baseline offsets with asterisks.
    '''
    # Hack to adjust the height of the asterisks to correct baseline
    t = ax.text(0.5, 0.5, 'R', fontsize=fontsize)
    r = fig.canvas.get_renderer()
    bb = t.get_window_extent(renderer=r)
    text_height = bb.height
    text_width = bb.width
    ax_height = ax.get_window_extent(renderer=r).height
    ymin, ymax = ax.get_ylim()
    data_range = abs(ymax-ymin)
    # This is the height of the letter R in the axis units
    R_height = text_height/ax_height
    # To get the height in data units
    R_h_data = R_height*data_range
    del ax.texts[-1]
    return R_h_data

def sc_swarmplot(data=None, all_genes=None, *, x=None, y=None, hue=None, hue_name=None, order=None, y_excluded=[], test_y=None,
                 x_lab=None, y_lab=None, add_n_numbers=True, ax=None, fig=None, s=3, palette=[color_dict['grey'], color_dict['blue']], 
                 enrich_hm=True, cbar_lab_sp=5, hm_lstart=0.65, **kwargs):
    '''
    Plot single-cell sequencing data with specific gene class overlaid.
    all_genes = df containing values for all genes (even unplotted ones). Used to compare to the groups in data.
    test_y = the name of the y_var to use for statistical testing between groups.
    Show p-values for half-lives being different.
    Show p-values for enrichment of gene values in the group.
    * in the function definition means that every following argument must have
    a keyword.
    y_excluded is a list of categories to leave out of the plot. For example, if bg or none is included in the df.
    This version has a hack to get tthe points with hue to plot on top of the background points
    I cannot figure out how to get the axis to update dynamically with jupyter, so just going to plot both and save the second one
    https://stackoverflow.com/questions/39658717/plot-dynamically-changing-graph-using-matplotlib-in-jupyter-notebook
    To get around this, I pass the figure (which has the size I want), copy it and put a new axis on the copied figure
    '''

    if order is None:
        order = data.groupby(y)[x].median().sort_values().index
        # a = (order > 50).values
        # np.argmax(a)
    order2 = [i for i in order if i not in y_excluded]

    if ax is None:
        ax = plt.gca()

    # Drop the categories in y_excluded
    data = data.query(f'{y}!=@y_excluded')

    # First use seaborn swarm to get the jitter positions for the data
    ax1 = sns.swarmplot(data=data, x=x, y=y, order=order2, orient='h', hue=hue, hue_order=[False, True],
                       palette=palette, ax=ax, s=2, **kwargs)

    # hack to get the attributes for the figure and the axis and replot
    figx, figy = fig.bbox.x1/100, fig.bbox.y1/100
    left_x = ax1.bbox.x0/fig.bbox.x1
    bottom_y = ax1.bbox.y0/fig.bbox.y1
    width = ax1.bbox.x1/fig.bbox.x1 - left_x
    height = ax1.bbox.y1/fig.bbox.y1 - bottom_y
    axes = (left_x, bottom_y, width, height)

    fig = plt.figure(figsize=(figx,figy))
    ax = fig.add_axes(axes)

    hue_rgba = colors.to_rgba(palette[1])
    i = 0
    bg_data_all = []
    colored_data_all = []
    for c in ax1.collections:
        # Filter out collections which are not scatterplots
        if c.get_offsets().shape[0] > 0:
            a_colors = c.get_fc()
            a_offsets = c.get_offsets()
            colored_idx = np.where((a_colors == hue_rgba).all(axis=1))[0]
            colored_data = a_offsets[colored_idx]

            mask = np.ones(a_offsets.shape[0], dtype=bool)
            mask[colored_idx] = False
            bg_data = a_offsets[mask]
            assert bg_data.shape[0] + colored_data.shape[0] == a_offsets.shape[0]
            bg_data_all.append(bg_data)
            colored_data_all.append(colored_data)
            i+= 1

    # Replot data with the hue colors on top:
    n_datasets = len(bg_data_all)
    # e.g. (16.7, -0.7), but this stretches the heatmap out a little bit
    # ax.set_ylim(n_datasets - 1 + 0.7, -0.7)
    ax.set_ylim(n_datasets - 1 + 0.5, -0.5)

    ax.yaxis.set_major_locator(plticker.MultipleLocator(base=1.0))
    for i in range(len(bg_data_all)):
        ax.scatter(bg_data_all[i][:,0], bg_data_all[i][:,1], color=palette[0], s=s)
        ax.scatter(colored_data_all[i][:,0], colored_data_all[i][:,1], color=palette[1], s=s)

    # ticks have to be set after plotting, otherwise they get re-set by subsequent plotting
    ax.set_yticks(range(0, n_datasets))

    _ = ax.set_xlabel(x_lab)
    _ = ax.set_ylabel(y_lab)

    ax.text(0.5, 1, hue_name, transform=ax.transAxes, color=palette[-1], ha='center', weight='bold')

    # Add the gene numbers to the cell type genes
    if add_n_numbers:
        gene_nums = data[y].value_counts().loc[order2].values
        new_labels = ['%s (%s)' % (i, j) for i,j in zip(order2, gene_nums)]
        ax.set_yticklabels(new_labels)
    
    else:
        ax.set_yticklabels(order2)

    # Get p-values for the comparison of deg rates between the cell type group and the other genes
    p_vals = []
    for i in order2:
        s_df = data.query(f'{y}==@i')
        group_genes = s_df.index
        g1 = all_genes.loc[group_genes, 'deg_rate']
        g2 = all_genes.loc[~all_genes.index.isin(group_genes), 'deg_rate']

        _, p1 = stats.mannwhitneyu(g1, g2)
        p_vals.append(p1)
    stars = sig_2_stars(p_vals)
    stars = np.where(stars=='ns', '', stars)

    # Put p_vals on the y-axis for having sig. difference to the stability of the other groups
    fontsize = 10
    letter_h = get_letter_height(ax, plt.gcf(), fontsize=fontsize)
    for i in range(0, len(stars)):
        ax.text(110, i+letter_h*0.8, stars[i], ha='center', va='bottom', weight='bold', fontsize=fontsize)

    res_tf = None
    if enrich_hm:
        # prev, lstart=0.65
        res_tf = enrich_heatmap(data=data, all_genes=all_genes, x='stab_percentile', y=y, hue=['TF'], y_lab1='fraction of genes',
                   y_lab2='-log'r'$_{10}$'' p-value', hstart=bottom_y, lstart=hm_lstart, fig=fig, xticklabs1=['counts'], 
                   xticklabs2=['enrichment'], hm_xlab=False, order=order2, cbar_lab_sp=cbar_lab_sp)

    return ax, res_tf, p_vals

def enrich_heatmap(data=None, all_genes=None, x=None, y=None, hue=None, hue_name=None, test_y=None, x_lab=None, y_lab1=None,
                   y_lab2=None, xticklabs1=['counts'], xticklabs2=['enrichment'], fig=None, lstart=0.68, hstart=0.15, hm_width=0.06, 
                   cb_width=0.03, ylabels=False, hm_xlab=True, order=None, cbar_lab_sp=4.5):
    '''
    Plot fraction of group and enrichment stats.
    Generally this will be plotted to the right side of the swarmplot.
    cbar_lab_sp is used to space the cbar label since align_ylabels() doesn't seem to work on the right.
    '''
    if order is None:
        order = data.groupby(y)[x].median().sort_values().index
    # Add the gene numbers to the cell type genes
    gene_nums = data[y].value_counts().loc[order].values
    new_labels = ['%s (%s)' % (i, j) for i,j in zip(order, gene_nums)]

    h = 0.95 - hstart
    cbar_h = (0.95 - hstart - 0.04)/2
    # Now calculate and plot the enrichment heatmaps to the right
    ax1 = fig.add_axes((lstart, hstart, hm_width, h))
    ax2 = fig.add_axes((lstart+hm_width+0.02, hstart, hm_width, h))
    ax3 = fig.add_axes((lstart+hm_width*2+0.04, hstart, cb_width, cbar_h))
    ax4 = fig.add_axes((lstart+hm_width*2+0.04, hstart+cbar_h+0.04, cb_width, cbar_h))

    total_genes = len(all_genes)
    # Adapt this to work with multiple vals for the 'hue' variable, i.e. hue = [TFs, RBPs]
    dfs = []
    for i in hue:
        enrich_pvals = []
        frac_genes = []
        n_genes = []
        for j in order:
            l1 = data.query(f'{y}==@j').index
            l2 = all_genes.query(i).index
            odds_r, log_p, lower, upper, table = calc_fisher_exact(l1, l2, total_genes)
            enrich_pvals.append(log_p)
            n_ol = len(l1.intersection(l2))
            frac_genes.append(n_ol/len(l1))
            n_genes.append(n_ol)
        enrich_df = pd.DataFrame({'frac_genes':frac_genes, 'n_genes':n_genes, 'enrich_p':enrich_pvals})
        dfs.append(enrich_df)

    # Build the frac genes heatmap and the enrichment heatmaps
    frac_hm = pd.concat([i['frac_genes'] for i in dfs], axis=1).values
    count_hm = pd.concat([i['n_genes'] for i in dfs], axis=1).values
    enrich_hm = pd.concat([i['enrich_p'] for i in dfs], axis=1).values

    # Sns heatmap expects a 2D array
    arrays = {'frac':frac_hm, 'count':count_hm, 'enrich':enrich_hm}
    for a in arrays:
        if arrays[a].ndim == 1:
            arrays[a].reshape(-1,1)

    # Replace np.nan with 0 for enrich_pvals in order to plot
    arrays['enrich'] = np.nan_to_num(arrays['enrich'])
    ax1 = sns.heatmap(arrays['frac'], annot=arrays['count'], fmt='d',
                            ax=ax1, cmap='viridis', cbar=False)
    ax2 = sns.heatmap(arrays['enrich'], ax=ax2, cmap='magma', cbar=False)

    # Plot the ylabels if necessary
    # It plots them rotated unless you specify rotation=0
    if ylabels == False:
        ax1.set_yticks([])
    else:
        ax1.set_yticklabels(new_labels, rotation=0)
        ax1.set_ylabel('cell type (num genes)')
    ax2.set_yticks([])

    tick_pos = np.arange(0.5, len(hue), 1)
    ax1.set_xticks(tick_pos)
    ax2.set_xticks(tick_pos)

    ax1.set_xticklabels(xticklabs1, rotation=45)
    ax2.set_xticklabels(xticklabs2, rotation=45)

    if hm_xlab:
        ax1.set_xlabel('counts', loc='right')
        ax2.set_xlabel('enrichment', loc='left')        

    norm1 = plt.Normalize(arrays['frac'].min(), arrays['frac'].max())
    sm1 = plt.cm.ScalarMappable(cmap='viridis', norm=norm1)

    norm2 = plt.Normalize(arrays['enrich'].min(), arrays['enrich'].max())
    sm2 = plt.cm.ScalarMappable(cmap='magma', norm=norm2)

    # https://stackoverflow.com/questions/32462881/add-colorbar-to-existing-axis
    # The top colorbar, fraction of genes 
    fig.colorbar(sm1, cax=ax4, orientation='vertical')
    ax4.set_ylabel(y_lab1)

    # The bottom colorbar, enrichment
    fig.colorbar(sm2, cax=ax3, orientation='vertical')
    ax3.set_ylabel(y_lab2)

    # fig.align_ylabels() # doesn't seem to work with labels on the right

    ax3.yaxis.set_label_coords(cbar_lab_sp, 0.5)
    ax4.yaxis.set_label_coords(cbar_lab_sp, 0.5)
    arrays['order'] = order
    return arrays

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
        lw = 0.75
        # in matplotlib 3.5, now ax.artists is empty
        boxes = ax.patches
        for box in boxes:
            r, g, b, a = box.get_facecolor()
            box.set_facecolor((r, g, b, box_alpha))
            box.set_edgecolor('k')
            box.set_lw(lw)
        for line in ax.lines:
            line.set_color('k')
            line.set_lw(lw)
        # Need to replot legend to add the transparent colors
        # No artists with labels found to put in legend.  
        # Note that artists whose label start with an underscore are ignored when legend() is called with no argument.
        # How can we find out if any of the artists have labels?
        handles, labels = ax.get_legend_handles_labels()
        if handles != []:
            if 'hue' in kwargs:
                ax.legend(title=kwargs['hue'])
            else:
                ax.legend()
        return ax

class Connector:
    '''Add line showing comparison'''
    def __init__(self, ax):
        self.gap = (ax.get_ylim()[1] - ax.get_ylim()[0])*.05
        self.top = ax.get_ylim()[1] + self.gap

    def add_connector(self, ax, xvals):
        '''Pass a list of x vals to use for drawing a horizontal line'''
        line = lines.Line2D(xvals, [self.top, self.top])
        line.set_clip_on(False)
        ax.add_line(line)
        self.top += self.gap

def compare_experiments(df, experiments=None, id_col=None, val_col=None, other_cols=[], pseudo=0, 
                         read_col='summed_est_counts', log=True):
    '''
    df is a "flat" dataframe with only one index level and the id_val in the columns
    experiments is a list of dictionaries with columns specifying the experiment and data labels,
    position 0 will be the x, position 1 will be the y
    i.e. [{'RNAtype':'input', 'condition':'ph'}, {'RNAtype':'input', 'condtion':'mock'}]
    if pseduo='min', then find the minimum value accross the set of experiments to use as the pseudocount.
    For version 2: 
    add other_cols arg to output other columns that can also be used for filtering later.
    '''
    dfs = []
    for i in experiments:
        query_str = '&'.join([f'{k} == "{v}"' if type(v) == str else f'{k} == {v}' for k,v in i.items()])
        # This only works if all the search values are strings, but not if float or int
        # query_str = '&'.join([f'{k} == "{v}"' for k,v in i.items()])
        this_df = df.query(query_str)[[id_col, val_col, *other_cols]].copy()
        dfs.append(this_df)
    cdf = pd.merge(*dfs, left_on=id_col, right_on=id_col, suffixes=('_x', '_y'))
    val_cols = [f'{val_col}_x', f'{val_col}_y']
    if pseudo == 'min':
        a = np.unique(pd.concat([cdf[val_cols[0]], cdf[val_cols[1]]]))
        lowest_nonzero = a[a>0][0]        
        pseudo = lowest_nonzero
    
    cdf[val_cols] += pseudo
    if log:
        cdf[val_cols] = cdf[val_cols].apply(np.log10)
    return cdf

def plot_scatter(df, experiments=None, id_col=None, genegroup=None, grouplabel=None, xlabel=None, ylabel=None, rsquare=True, 
                 loc=None, diagonal=None, histscatter=False, limits=None, same_max=False, ax=None):
    '''
    Plot scatter points in light grey and overlay members of a gene group.
    df is a "flat" dataframe with only one index level and the id_val in the columns
    experiments is a list of column names
    position 0 will be the x, position 1 will be the y
    limits = dictionary with x and y limits, i.e. {'x':{'left':-2, 'right':4}, 'y':{'bottom':-2, 'top':4}}
    same_max means that both axis will have the same max value
    '''
    xname, yname = experiments
    x = df[xname]
    y = df[yname]
    if histscatter:
        ax = sns.histplot(data=df, x=x, y=y, cmap='rocket', ax=ax, zorder=2)
    else:
        ax.scatter(x, y, s=5, color='k', alpha=0.3, ec='none')

    if genegroup:
        df['group'] = df[id_col].isin(genegroup)
        x2 = df.query('group')[xname]
        y2 = df.query('group')[yname]
        ax.scatter(x2, y2, s=5, color=color_dict['purple'], ec='none', label=grouplabel)
    if rsquare:
        rval, pval = stats.pearsonr(x, y)
        r2_val_av = rval**2
        ax.text(0.1, 0.9, 'r'r'$^2$'' = %1.2f' % r2_val_av, fontsize = 8, transform=ax.transAxes)
    if loc:
        ax.xaxis.set_major_locator(loc)
        ax.yaxis.set_major_locator(loc)
    if same_max:
        maxlim = max(ax.get_ylim()[1], ax.get_xlim()[1])
        ax.set_xlim(right=maxlim)
        ax.set_ylim(top=maxlim)
    if limits:
        if 'x' in limits:
            ax.set_xlim(**limits['x'])
        if 'y' in limits:
            ax.set_ylim(**limits['y'])
    if diagonal:
        ax.axline((0, 0), slope=1, color=color_dict['grey'], alpha = 0.8, zorder=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
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

def sig_2_stars(pvals, siglevels=np.array([0.05, 10e-10, 10e-40]), tex=False):
    '''
    Convert array to asterisks. If tex=true, use the mathematical rep. of \ast.
    Problems: math asterisks look frilly and have too much space between them
    '''
    if tex:
        sigstars = np.array(['ns', r'\ast', r'\ast\ast', r'\ast\ast\ast'])
    else:
        sigstars = np.array(['ns', '*', '**', '***'])
    # this will return 0 if the value is larger than the first one
    j = np.digitize(pvals, siglevels)
    stars = sigstars[j]
    return stars

def add_stars(lefti, righti, starty=None, height=0.1, p=None, ax=None, siglevels=np.array([0.05, 10e-10, 10e-40]), col='k', write_p=False):
    '''
    Input the 0-based index of the first and the second groups to draw the line between (lefti, righti).
    Add stars to the axis denoting significance level.
    Siglevels is the significance levels where p<level to give it a star
    '''
    stars = sig_2_stars(p)
    if write_p:
        if p:
            x = Decimal(p)
            dec = '{:.1e}'.format(x)
            stars = f'{stars}p={dec}'
    box_pos = get_box_coords(ax)
    leftx = box_pos[lefti]
    rightx = box_pos[righti]
    xdata, ydata = [[leftx, leftx, rightx, rightx], [starty, starty+height, starty+height, starty]]
    line = lines.Line2D(xdata, ydata, lw=0.75, c=col, alpha=0.5)
    line.set_clip_on(False)
    ax.add_line(line)
    # Asterix already seems to be pretty high, if you give it starty, compensates for that.
    star_vals = ['*', '**', '***', '****']
    if stars not in star_vals:
        starty += (ax.get_ylim()[-1] - ax.get_ylim()[0])*.03
    ax.text((leftx+rightx)*.5, starty, stars, ha='center', va='bottom', color=col, alpha=0.5)

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

def diagonal_cuts(top_ax, bottom_ax, d=0.015):
    '''
    Draw diagonal lines to split the axis.
    https://stackoverflow.com/questions/63726234/how-to-draw-a-broken-y-axis-catplot-graphes-with-seaborn
    '''
    sns.despine(ax=top_ax, bottom=True)
    top_ax.set_xticks([])
    kwargs = dict(transform=top_ax.transAxes, color='k', clip_on=False)
    top_ax.plot((-d, +d), (-d, +d), **kwargs)
    kwargs.update(transform=bottom_ax.transAxes)
    bottom_ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)

def colval_fmt(colval):
    '''
    Pandas needs the column value to have surrounding quotes "" if a string and not if True/False
    '''
    if isinstance(colval, bool):
        return colval
    else:
        return f'"{colval}"'

def colname_fmt(colname):
    if len(colname.split(' ')) > 1:
        return f'`{colname}`'
    else:
        return colname
    
def get_group_Ns(df, x, hue=None, order=[False, True], hue_order=[False, True], ticklabels=None):
    '''
    Return the n-numbers for the groups in the plot
    Specify the order and hue_order if not [False, True]
    '''
    res = []
    if hue is None:
        for i in order:
            l = [(x, i)]
            query_str = '&'.join([f'{colname_fmt(p[0])}=={colval_fmt(p[1])}' for p in l])
            n_i = len(df.query(query_str))
            res.append(n_i)
    else:
        for i in order:
            hue_vals = []
            for j in hue_order:
                l = [(x, i), (hue, j)]
                query_str = '&'.join([f'{colname_fmt(p[0])}=={colval_fmt(p[1])}' for p in l])
                n_i = len(df.query(query_str))
                hue_vals.append(n_i)
            res.append(tuple(hue_vals))

    if ticklabels is None:
        return res
    else:
        new_labels = [f'{i[0]}\nn = {i[1]}' for i in zip(ticklabels, res)]
        return new_labels

def enrich_table(res, pvals):
    show_enrich = pd.DataFrame([res['order'], res['enrich'], -np.log10(pvals)]).transpose()
    show_enrich.columns = ['cell type', '-log10pval_enrich', '-log10pval_stability']
    return show_enrich