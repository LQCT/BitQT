#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Roy Gonzalez-Aleman                        [rglez.developer@gmail.com]
@author: Daniel Platero Rochart                          [dplatero97@gmail.com]
"""

import sys
import pickle
import argparse
from collections import OrderedDict

import numpy as np
import mdtraj as md
from bitarray import bitarray as ba


def parse_arguments():
    """
    Parse user arguments from the command line.

    Returns
    -------
    user_inputs : parser.argparse
        namespace containing user input arguments.
    """
    # Initializing argparse ---------------------------------------------------
    desc = '\nBitQT: QT heuristic clustering for MD trajectories'
    parser = argparse.ArgumentParser(description=desc,
                                     add_help=True,
                                     epilog='As simple as that ;)')
    # Arguments: loading trajectory -------------------------------------------
    parser.add_argument('-top', dest='topology', action='store',
                        help='path to topology file (psf/pdb)',
                        type=str, required=False)
    parser.add_argument('-traj', dest='trajectory', action='store',
                        help='path to trajectory file',
                        type=str)
    parser.add_argument('-first', dest='first', action='store',
                        help='first frame to analyze (starting from 0)',
                        type=int, required=False, default=0)
    parser.add_argument('-last', dest='last', action='store',
                        help='last frame to analyze (starting from 0)',
                        type=int, required=False, default=None)
    parser.add_argument('-stride', dest='stride', action='store',
                        help='stride of frames to analyze',
                        type=int, required=False, default=1)
    parser.add_argument('-sel', dest='selection', action='store',
                        help='atom selection (MDTraj syntax)',
                        type=str, required=False, default='all')
    parser.add_argument('-rmwat', dest='remove_waters', action='store',
                        help='remove waters from trajectory?',
                        type=bool, required=False, default=0,
                        choices=[True, False])
    # Arguments: clustering ---------------------------------------------------
    parser.add_argument('-cutoff', action='store', dest='cutoff',
                        help='RMSD cutoff for pairwise comparisons in A',
                        type=float, required=False, default=1.0)
    parser.add_argument('-minsize', action='store', dest='minsize',
                        help='minim number of frames inside returned clusters',
                        type=int, required=False, default=2)
    parser.add_argument('-ref', action='store', dest='reference',
                        help='reference frame to align trajectory',
                        type=int, required=False, default=0)
    # Arguments: analysis -----------------------------------------------------
    parser.add_argument('-odir', action='store', dest='outdir',
                        help='output directory to store analysis',
                        type=str, required=False, default='./')
    user_inputs = parser.parse_args()
    return user_inputs


# @profile
def load_trajectory(args):
    """
    Load trajectory file using MDTraj. If trajectory format is h5, lh5 or
    pdb, a topology file is not required. Otherwise, you should specify a
    topology file.

    Parameters
    ----------
    args : argparse.Namespace
        user input parameters parsed by argparse (CLI).

    Returns
    -------
    traj : mdtraj.Trajectory
        MDTraj trajectory object.
    """
    traj_file = args.trajectory
    traj_ext = traj_file.split('.')[-1]
    # Does trajectory file format need topology ? -----------------------------
    if traj_ext in ['h5', 'lh5', 'pdb']:
        traj = md.load(traj_file)
    else:
        traj = md.load(traj_file, top=args.topology)

    # Reduce RAM consumption by loading selected atoms only -------------------
    if args.selection != 'all':
        try:
            sel_indx = traj.topology.select(args.selection)
        except ValueError:
            print('Specified selection is invalid')
            sys.exit()
        if sel_indx.size == 0:
            print('Specified selection in your system corresponds to no atoms')
            sys.exit()
        traj = traj.atom_slice(sel_indx)[args.first:args.last:args.stride]
    else:
        traj = traj[0:-1:1]
    return traj


# @profile
def calc_rmsd_matrix(trajectory, args):
    """
    Calculate optimal RMSD binary-encoded square matrix using MDTraj. Pairwise
    similarity is saved in RAM as bits (dict of bitarrays), not floats.

    Parameters
    ----------
    trajectory : mdtraj.Trajectory
        MDTraj trajectory object.
    args : argparse.Namespace
        user input parameters parsed by argparse (CLI)..

    Returns
    -------
    matrix : collections.OrderedDict
        dict of bitarrays.

    """
    trajectory.center_coordinates()
    cutoff = np.full(trajectory.n_frames, args.cutoff / 10, dtype=np.float32)
    matrix = OrderedDict()
    to_explore = range(trajectory.n_frames)
    for i in to_explore:
        rmsd_ = md.rmsd(trajectory, trajectory, i, precentered=True)
        vector_np = np.less_equal(rmsd_, cutoff)
        bitarr = ba()
        bitarr.pack(vector_np.tobytes())
        bitarr.fill()
        matrix.update({i: bitarr})
    return matrix


# @profile
def calc_matrix_degrees(unclustered_bit, matrix):
    """
    Calculate number of neighbors (degree) of unclustered nodes in matrix.

    Parameters
    ----------
    unclustered_bit : bitarray.bitarray
        bitarray with indices of unclustered nodes turned on.
    matrix : collections.OrderedDict
        dict of bitarrays.

    Returns
    -------
    degrees : numpy.ndarray
        array containing each node degree. Clustered nodes have degree = 0.

    """
    one = ba('1')
    degrees = np.zeros(len(unclustered_bit), dtype=np.int32)
    for node in unclustered_bit.itersearch(one):
        try:
            degrees[node] = matrix[node].count()
        except KeyError:
            pass
    return degrees


# @profile
def colour_matrix(degrees, matrix):
    """
    Greedy coloring of bit-encoded RMSD matrix.

    Parameters
    ----------
    degrees : numpy.ndarray
        array containing each node degree. Clustered nodes have degree = 0.
    matrix : collections.OrderedDict
        dict of bitarrays.

    Returns
    -------
    colors : numpy.ndarray
        array of colors assigned to each node of the matrix.
    """
    # Constants ---------------------------------------------------------------
    N = degrees.size
    m = len(matrix)
    one = ba('1')
    xcolor = 0
    # Initialize containers ---------------------------------------------------
    ordered_by_degrees = iter((-degrees[:m]).argsort())
    colors = np.zeros(N, dtype=np.int32)
    colored = ba(N)
    colored.setall(0)
    seen = set()
    while True:
        # Retrieve the max-degree node ----------------------------------------
        max_node = next(ordered_by_degrees)
        if max_node in seen:
            continue
        else:
            seen.add(max_node)
            xcolor += 1
            not_neighbors = ~ matrix[max_node]
            not_colored = ~colored
            candidates = not_neighbors & not_colored
            # Nodes passing conditions (not-neighb, not-colored, not-neighb) --
            passed = [max_node]
            for candidate in candidates.itersearch(one):
                passed.append(candidate)
                try:
                    candidates &= ~matrix[candidate]
                except KeyError:
                    continue
                if not candidates.any():
                    break
            seen.update(passed)
            # Deliver a color class to passed nodes ---------------------------
            colors[passed] = xcolor
            colored = ba()
            colored.pack(colors.astype(np.bool).tobytes())
            if colored.count(0) == 0:
                break
    return colors


# @profile
def bitarray_to_np(bitarr):
    """
    Convert from bitarray.bitarray to numpy.ndarray efficiently.

    Parameters
    ----------
    bitarr : bitarray.bitarray
        a bitarray.

    Returns
    -------
    numpy.ndarray
        boolean bitarray equivalent to the binary bitarray input object.
    """
    return np.unpackbits(bitarr).astype(np.bool)


# @profile
def do_bit_cascade(big_node, degrees, colors, matrix, max_):
    """
    Perform succesive AND operations between an initial bitarray and subsequent
    bitarray candidates to search for a clique.

    Parameters
    ----------
    big_node : int
        node whose bitarray will start the operations.
    degrees : numpy.ndarray
        array containing each node degree. Clustered nodes have degree = 0.
    colors : numpy.ndarray
        array of colors assigned to each node of the matrix.
    clustered_bit : bitarray.bitarray
        bitarray with indices of clustered nodes turned on.
    matrix : collections.OrderedDict
        dict of bitarrays.
    max_ : int
        Stop iterative AND operations after the initial bitarray has max_
        bits turned on.

    Returns
    -------
    init_cascade : bitarray.bitarray
        initial bitarray before any AND operation.
    ar : numpy.ndarray
        array of nodes forming a clique.
    """
    init_cascade = matrix[big_node]
    # .... recovering neighbors and their information .........................
    neighb = bitarray_to_np(init_cascade).nonzero()[0]
    neighb_colors = colors[neighb]
    if len(set(neighb_colors.tolist())) <= max_:
        return None
    neighb_degrees = degrees[neighb]
    g = np.bincount(neighb_colors)
    neighb_g = g[neighb_colors]
    # .... ordering neighbors by g ---> colors ---> degrees ...................
    idx = np.lexsort([-neighb_degrees, neighb_colors, neighb_g])
    candidates_info = zip(neighb[idx], neighb_colors[idx])

    # .... BitCascade considering divergence ..................................
    counter = 0
    seen = set()
    for candidate, color in candidates_info:
        if (color in seen) or (not init_cascade[candidate]):
            continue
        seen.add(color)
        init_cascade = matrix[candidate] & init_cascade
        counter += 1
        COUNT = init_cascade.count()
        if (COUNT <= max_):
            return None
            break
        elif counter >= COUNT:
            break
    ar = np.nonzero(np.unpackbits(init_cascade).astype(np.bool))[0]
    return init_cascade, ar


# @profile
def set_to_bitarray(set_, N):
    """
    Convert from python set to bitarray.bitarray.

    Parameters
    ----------
    set_ : set
        a python set.
    N : int
        lenght of the desired bitarray. It must be greater than the maximum
        value of indices present in set.

    Returns
    -------
    bitarr : bitarray.bitarray
        bitarray of lenght N with indices present in set turned on.
    """
    zero_arr = np.zeros(N, dtype=np.bool)
    zero_arr[list(set_)] = 1
    bitarr = ba()
    bitarr.pack(zero_arr.tobytes())
    return bitarr


def pickle_to_file(data, file_name):
    """
    Serialize data using the pickle library.

    Parameters
    ----------
    data : serializable object
        variable name of the object to serialize.
    file_name : str
        name of the pickle file to be created.

    Returns
    -------
    file_name : str
        name of the serialized file.
    """
    with open(file_name, 'wb') as file:
        pickle.dump(data, file)
    return file_name


def unpickle_from_file(file_name):
    """
    Unserialize data using the pickle library.

    Parameters
    ----------
    file_name : str
        name of the serialized file.

    Returns
    -------
    data : serializable object
        unserialized object.
    """
    with open(file_name, 'rb') as file:
        data = pickle.load(file)
    return data


def to_VMD(output_name, clusters_array, stride, trajectory):
    with open(output_name, 'wt') as clq:
        for num in np.unique(clusters_array):
            clq.write('{}:\n'.format(num))
            numcluster = np.where(clusters_array == num)[0]
            frames = [str(x * stride) for x in numcluster]
            members = ' '.join(frames)
            clq.write('Members: ' + members + '\n\n')
