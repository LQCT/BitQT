#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Roy Gonzalez-Aleman                    [rglez.developer@gmail.com]
@author: Daniel Platero Rochart                      [dplatero97@gmail.com]
"""

import time
from collections import deque, Counter

import numpy as np
from bitarray import bitarray as ba
from bitarray import util as bu

import bitfunctions as bf

# =============================================================================
# 0. CLI user parameters
# =============================================================================
args = bf.parse_arguments()

# ====> Debugging <============================================================
# import argparse
# args = argparse.Namespace()
# location = '/opt/0A-PUBLISHING/04-BitClique/bit_benchmark/trajs/'
# location = '/opt/0A-PUBLISHING/04-BitClique/bit_benchmark/trajs/'
# args.first = 0
# args.last = None
# args.stride = 1
# args.minsize = 2
# Change for each testing -----------------------------------------------------
# args.topology = location + 'melvin.pdb'
# args.trajectory = location + 'melvin.dcd'
# args.topology = location + 'aligned_tau.pdb'
# args.trajectory = location + 'aligned_original_tau_6K.dcd'
# args.selection = "backbone"
# args.cutoff = 4

# [6K, 30K, 50K, 100K]
# ["all", "name =~ '[^H.*]'", "name =~ '[^H.*]'", "backbone"]
# [4, 4, 2, 2]
# ====> Debugging <============================================================

#%

# @profile
# def main():
# =========================================================================
# 1. Creating binary matrix (adjacency list)
# =========================================================================
# ++++ Get adjacency matrix of trajectory as list of bitarrays ++++++++++++
trajectory = bf.load_trajectory(args)
init_time = time.time()
matrix = bf.calc_rmsd_matrix(trajectory, args)

# matrix = bf.unpickle_from_file('matrix_tau.pick')
# matrix = bf.unpickle_from_file('matrix_melvin.pick')
# matrix = gnr.unpickle_from_file('matrix_k50.pick')
# matrix = bf.unpickle_from_file('matrix_k100.pick')

matrix_time = time.time() - init_time
# ++++ Tracking clust/unclustered bits to avoid re-computations +++++++++++
N = len(matrix[0])
m = len(matrix)
unclust_bit = ba(N)
unclust_bit.setall(1)
clustered_bit = unclust_bit.copy()
clustered_bit.setall(0)
zeros = np.zeros(N, dtype=np.int32)
# ++++ Save clusters in an array (1 .. N) +++++++++++++++++++++++++++++++++
clusters_array = np.zeros(N, dtype=np.int32)
ncluster = 0
clustered = set()
nmembers = []
# ++++ Coloring ordered vertices (1 .. N) +++++++++++++++++++++++++++++++++
degrees = bf.calc_matrix_degrees(unclust_bit, matrix)
ordered_by_degs = degrees.argsort()[::-1]
colors = bf.colour_matrix(ordered_by_degs, matrix)
# colors[np.frombuffer(clustered_bit.unpack(), dtype=np.bool)] = 0

# =========================================================================
# 2. Main algorithm: BigClique !
# =========================================================================
while True:
    ncluster += 1
    # ++++ Find a big clique early ++++++++++++++++++++++++++++++++++++++++
    big_node = degrees.argmax()
    bit_clique, big_clique = bf.do_bit_cascade(big_node, degrees, colors,
                                               matrix, 0)
    big_clique_size = big_clique.size
    # ++++ Find promising nodes +++++++++++++++++++++++++++++++++++++++++++
    biggers = degrees > big_clique_size
    biggers[big_clique] = False
    cluster_colors = colors[big_clique]
    biggers_colors = colors[biggers]
    promising_colors = np.setdiff1d(biggers_colors, cluster_colors)
    promising_nodes = deque()
    for x in promising_colors:
        promising_nodes.extend(((colors == x) & biggers).nonzero()[0])
    # ++++ Explore all promising nodes ++++++++++++++++++++++++++++++++++++
    cum_found = big_clique
    while promising_nodes:
        node = promising_nodes.popleft()
        try:
            bit_clique, clique = bf.do_bit_cascade(node, degrees, colors,
                                                   matrix, big_clique_size)
            clique_size = len(clique)
        except TypeError:
            clique_size = 0
        # ++++ Cumulative update only if biggers candidates are found +++++
        if clique_size > big_clique_size:
            big_node = node
            big_clique = clique
            big_clique_size = big_clique.size
            # ++++ Repeat previous condition ++++++++++++++++++++++++++++++
            cum_found = np.concatenate((cum_found, big_clique))
            biggers = degrees > big_clique_size
            biggers[cum_found] = False
            cluster_colors = colors[big_clique]
            biggers_colors = colors[biggers]
            promising_colors = np.setdiff1d(biggers_colors, cluster_colors)
            promising_nodes = deque()
            for x in promising_colors:
                promising_nodes.extend(((colors == x) & biggers).nonzero()[0])
    nmembers.append(big_clique_size)

    # ++++ Save new cluster & update ncluster +++++++++++++++++++++++++++++
    clusters_array[big_clique] = ncluster
    # ++++ Update (un)clustered_bit +++++++++++++++++++++++++++++++++++++++
    clustered.update(big_clique)
    clustered_bit = bf.set_to_bitarray(clustered, N)
    unclust_bit = ~clustered_bit

    # ++++ Hard erasing of clustered frames from matrix +++++++++++++++++++
    degrees = zeros.copy()
    for x in unclust_bit[:m].itersearch(ba('1')):
        degrees[x] = matrix[x].count()
        if bu.count_and(matrix[x], clustered_bit):
            matrix[x] &= (matrix[x] ^ clustered_bit)

    # ++++ Stopping the algorithm +++++++++++++++++++++++++++++++++++++++++
    # if ncluster == args.max_clusters:
    #     break
    # if (big_clique_size < args.max_size) or (big_clique_size < 2):
    #     break
    if (big_clique_size < 2):
        break
total_time = time.time() - init_time

# =========================================================================
# 3. Observables
# =========================================================================
fluct = [nmembers[i] - nmembers[i + 1]
         for i, x in enumerate(nmembers) if i < len(nmembers) - 1]
neg_fluct = [x for x in fluct if x < 0]
nflucts = sorted(Counter(neg_fluct).most_common())

print('Fluct: {}. Distributed as: {}'.format(len(neg_fluct), nflucts))
print('Timing {} sec for matrix and {} sec total'.format(matrix_time,
                                                         total_time))
print('NClusters: {}'.format(len(nmembers)))
print('10Firsts: {}'.format(nmembers[:10]))
# bf.to_VMD(args.outdir, clusters_array, 1, args.trajectory)
bf.pickle_to_file(clusters_array, ('TEST.pickle'))


# main()
