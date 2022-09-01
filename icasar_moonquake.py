#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 09:46:35 2020

@author: matthew
"""
import sys

# print('Temporary import of matrix_show for debugging.  ')
# sys.path.append('/home/matthew/university_work/python_stuff/python_scripts/')
# from small_plot_functions import matrix_show

sys.path.append('/Users/TheStuffofAlice/src/ICASAR-master')
from ICASAR_functions import ICASAR, bootstrapped_sources_to_centrotypes
from auxiliary_functions import prepare_point_colours_for_2d, prepare_legends_for_2d, plot_2d_interactive_fig

import numpy as np
import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt

#%% Things to set
# 460 490 worked for 20 events
# 100,595 S16 peaked
ICASAR_settings = {"n_comp" : 4,                                    # number of PCs used
                   "tsne_param" : (30,12),                        # (perplexity, early_exaggeration)
                   "hdbscan_param" : (100,295),                      # (min_cluster_size, min_samples) Discussed in more detail in Mcinnes et al. (2017). min_cluster_size sets the smallest collection of points that can be considered a cluster. min_samples sets how conservative the clustering is. With larger values, more points will be considered noise.
                   "figures" : "png+window"}                       # if png, saved in a folder as .png.  If window, open as interactive matplotlib figures,


#%% load the data and visualise a few

moonquakes_r3 = np.load('U12_MHN_peaked.npy')[:,:,:ICASAR_settings['n_comp']]                              # only load the number of sources that we want, and get rid of last 3 (redundant) dimensions.
n_bootstrap, n_times, n_pcs = moonquakes_r3.shape                                                                   # nice little cube of signals
                                                                                                                    # Alice, am I right in thinking we have 101 bootstrapped runs, each signal has 2650 times, and I've selcted 5 PCs from each of these?
for pc_n in range(ICASAR_settings['n_comp']):
    fig, axes = plt.subplots(3,3)
    fig.canvas.set_window_title("fig_01_PC{pc_n}s")
    for ax in np.ravel(axes):
        ax.plot(np.arange(0,n_times), moonquakes_r3[np.random.randint(0,n_bootstrap),:,pc_n])


#%% Create a single rank 2 array of all the pcs together and plot some at random.
  
moonquakes_r2 = np.zeros((n_bootstrap * n_pcs, n_times))                                                        # which PC is which will be lost in the rank 2 structure.  
pc_n = 0
for i in range(n_bootstrap):                                                                                # loop through, copying signals into a new array.  
    for j in range(n_pcs):
        moonquakes_r2[pc_n,] = moonquakes_r3[i,:,j]
        pc_n += 1

fig1, axes = plt.subplots(5,3)
fig1.canvas.set_window_title('fig_02_random_pcs')
for ax in np.ravel(axes):
    ax.plot(np.arange(0,n_times), moonquakes_r2[np.random.randint(0,moonquakes_r2.shape[0]),])
        
    
#%% Do the cluster and 2d manifold representation

S_best, labels_hdbscan, xy_tsne, clusters_by_max_Iq_no_noise, Iq  = bootstrapped_sources_to_centrotypes(moonquakes_r2, ICASAR_settings['hdbscan_param'], ICASAR_settings['tsne_param'],
                                                                                                        bootstrap_settings = (n_bootstrap, n_pcs), figures = 'window')


#%% Plot the results of previous step

labels_colours = prepare_point_colours_for_2d(labels_hdbscan, clusters_by_max_Iq_no_noise)                                          # make a list of colours so that each point with the same label has the same colour, and all noise points are grey
legend_dict = prepare_legends_for_2d(clusters_by_max_Iq_no_noise, Iq)    
   
plot_2d_labels = {'title' : '03_clustering_and_manifold_results',
                  'xlabel' : 'TSNE dimension 1',
                  'ylabel' : 'TSNE dimension 2'}

temporal_data_S_all = {'tcs_r2' : moonquakes_r2,
                       'xvals'  : np.arange(0, moonquakes_r2.shape[1]) }                                          # make a dictionary of the sources recovered from each run
plot_2d_interactive_fig(xy_tsne.T, colours = labels_colours, temporal_data = temporal_data_S_all,                 # make the 2d interactive plot
                        labels = plot_2d_labels, legend = legend_dict, arrow_length=2, inset_axes_side = {'x':0.3, 'y':0.1})


plt.show()
