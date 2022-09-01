import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from tslearn.preprocessing import TimeSeriesScalerMeanVariance
from tslearn.clustering import TimeSeriesKMeans
from collections import defaultdict
sys.path.append('../scripts/')
from plot_helpers2 import *

#Use Euclidean Kmeans time series clustering
#https://tslearn.readthedocs.io/en/stable/gen_modules/clustering/tslearn.clustering.TimeSeriesKMeans.html#tslearn.clustering.TimeSeriesKMeans
def cluster_ts(df, tpts, nclust, outname):
    seed = 0
    np.random.seed(seed)
    #numpy.random.shuffle(X_train) #Does the data need to be shuffled before clustering?
    #X = obs_pred_sig.values
    X = df.values
    X = TimeSeriesScalerMeanVariance().fit_transform(X)
    km = TimeSeriesKMeans(n_clusters=nclust, verbose=True, random_state=seed)
    y_pred = km.fit_predict(X)

    #Annotate genes with cluster info
    #Count genes per cluster
    genes_per_clust = {}
    for i in range(0, nclust):
        genes_per_clust[i] = (km.labels_ == i).sum()

    #Change cluster number for visualization purposes
    #clust_name_map: convert range 1->n to the original label from clustering result
    clust_name_map = {1:0, 2:1, 3:2, 4:3, 5:4, 6:5, 7:6, 8:7}
    #clust_name_map_r: convert original label to new label
    clust_name_map_r = {v:k for (k,v) in clust_name_map.items()}

    #Annotate genes with cluster assignment and renumbered cluster assignment
    df = df.copy()
    df['cluster_original'] = km.labels_
    df['cluster_renamed'] = df['cluster_original'].map(clust_name_map_r)
    
    genes_to_annotate = {}
    anns_to_add = defaultdict(list)
    for i in genes_to_annotate:
        anns_to_add[df.loc[i, 'cluster_original']].append(genes_to_annotate[i])

    #Specify where to write gene names on the plot
    genename_locs = {1:'topL', 2:'bottomR', 3:'bottomR', 4:'bottomR', 5:'topR', 6:'topR', 7:'bottomL', 8:'topL'}
    gene_pos = {'topL':{'pos':(0.02, 1), 'va':'top', 'ha':'left'}, 'bottomL':{'pos':(0.02, 0), 'va':'bottom', 'ha':'left'},
               'topR':{'pos':(1, 1), 'va':'top', 'ha':'right'}, 'bottomR':{'pos':(1, 0), 'va':'bottom', 'ha':'right'}}

    fig_fmt = 'landscape'
    if fig_fmt == 'portrait':
        fig = plt.figure(figsize = (dfig, dfig*2), constrained_layout = True, 
                         subplotpars=sp)
        gs = fig.add_gridspec(ncols = 2, nrows = 4)
        plot_order = [1, 5, 2, 6, 3, 7, 4, 8]
    if fig_fmt == 'landscape':
        fig = plt.figure(figsize = (dfig*2, dfig))
        gs = fig.add_gridspec(ncols = 4, nrows = 2)
        plot_order = range(1,nclust+1)

    p = 0
    for i in plot_order:
        yi = clust_name_map[i]
        ax = fig.add_subplot(gs[p])
        if anns_to_add[yi] != []:
            ann_string = '\n'.join(anns_to_add[yi])
            pos_dict = gene_pos[genename_locs[i]]
            ax.text(*pos_dict['pos'], ann_string, fontsize=6, va=pos_dict['va'], ha=pos_dict['ha'], style='italic', transform=ax.transAxes)

        for xx in X[y_pred == yi]:
            ax.plot(tpts, xx.ravel(), 'k-', alpha=.2)
        ax.plot(tpts, km.cluster_centers_[yi].ravel(), 'r-')
        ax.set_title('Clust %d (n = %d)' % (i, genes_per_clust[yi]), fontsize=6)
        ax.tick_params(bottom=True, left=False, labelbottom=True, labelleft=False, labelsize=6)
        ax.xaxis.set_major_locator(plticker.MultipleLocator(base=1.0))
        p+=1
    plt.savefig('%s.%s' % (outname, out_fmt), dpi = out_dpi)
    return df, X, km, y_pred, genes_per_clust