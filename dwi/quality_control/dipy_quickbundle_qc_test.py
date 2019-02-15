#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 13:13:48 2019

@author: brian
"""


import numpy as np
from nibabel import trackvis as tv
from dipy.tracking.streamline import Streamlines
from dipy.segment.clustering import QuickBundles
from dipy.io.pickles import save_pickle
from dipy.data import get_data
from dipy.viz import window, actor
from dipy.align.streamlinear import remove_clusters_by_size


# Move to subject directory
subdir = '/media/brian/hpc/data/Jeong-R01-Data/ESM-DTI_study/J0044/PreOP/VN/70deg/fsaverage/'

# Track name to cluster
trk_file = 'VN_70deg_qclen_fsa_18-2.trk'
exemplar_file = 'VN_70deg_qclen_fsa_18-2_exemplar_t5.trk'
streams, hdr = tv.read(subdir + trk_file)
streamlines = [i[0] for i in streams]

# Clustering threshold
thresh = 5.

qb = QuickBundles(threshold=thresh)
clusters = qb.cluster(streamlines)

# Drop small clusters
size_thresh = 10
big_clusters = clusters.get_large_clusters(size_thresh)
streamlines_pruned = [item for sublist in big_clusters for item in sublist]

# Plot original data
interactive = True
ren = window.Renderer()
ren.SetBackground(1, 1, 1)
ren.add(actor.streamtube(streamlines, window.colors.white))
window.record(ren, out_path='test_initial.png', size=(600, 600))
if interactive:
    window.show(ren)
    
# Plot streamlines colored according to centroid    
colormap = actor.create_colormap(np.arange(len(clusters)))
colormap_full = np.ones((len(streamlines), 3))
for cluster, color in zip(clusters, colormap):
    colormap_full[cluster.indices] = color

window.clear(ren)
ren.SetBackground(1, 1, 1)
ren.add(actor.streamtube(streamlines, colormap_full))
window.record(ren, out_path='fornix_clusters.png', size=(600, 600))
if interactive:
    window.show(ren)

# Just plot big clusters
pruned_inds=[]
colormap = actor.create_colormap(np.arange(len(big_clusters)))
colormap_full = np.ones((len(streamlines), 3))
for cluster, color in zip(big_clusters, colormap):
    colormap_full[cluster.indices] = color
    pruned_inds+=cluster.indices
    
window.clear(ren)
ren.SetBackground(1, 1, 1)
ren.add(actor.streamtube(streamlines_pruned, colormap_full[pruned_inds]))
window.record(ren, out_path='fornix_clusters.png', size=(600, 600))
if interactive:
    window.show(ren) 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    