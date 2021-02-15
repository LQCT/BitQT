#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Roy Gonzalez-Aleman                            [roy_gonzalez@fq.uh.cu]
@author: Daniel Platero Rochart                      [daniel.platero@gmail.com]
"""
import os
import pickle
import argparse
from collections import deque, OrderedDict

import numpy as np
import pandas as pd
import mdtraj as md
from bitarray import util as bu
from bitarray import bitarray as ba


valid_tops = set(['pdb', 'pdb.gz', 'h5', 'lh5', 'prmtop', 'parm7', 'prm7',
                  'psf', 'mol2', 'hoomdxml', 'gro', 'arc', 'hdf5', 'gsd'])


valid_trajs = set(['arc', 'dcd', 'binpos', 'xtc', 'trr', 'hdf5', 'h5', 'ncdf',
                   'netcdf', 'nc', 'pdb.gz', 'pdb', 'lh5', 'crd', 'mdcrd',
                   'inpcrd', 'restrt', 'rst7', 'ncrst', 'lammpstrj', 'dtr',
                   'stk', 'gro', 'xyz.gz', 'xyz', 'tng', 'xml', 'mol2',
                   'hoomdxml', 'gsd'])


def parse_arguments():
    """
    Parse arguments from the cli.

    Returns
    -------
    user_inputs : argparse.Namespace
        Namespace with all required arguments.

    """
    # Initializing argparse ---------------------------------------------------
    desc = '\nBitQT: A Graph-based Approach to the  Quality Threshold\
        Clustering of Molecular Dynamics'
    parser = argparse.ArgumentParser(prog='bitqt',
                                     description=desc,
                                     add_help=True,
                                     epilog='As simple as that ;)',
                                     allow_abbrev=False,
                                     usage='%(prog)s -traj trajectory [options]')
    # Arguments: loading trajectory -------------------------------------------
    traj = parser.add_argument_group(title='Trajectory options')
    traj.add_argument('-traj', dest='trajectory', action='store',
                      help='Path to trajectory file [required]', type=str,
                      metavar='trajectory', required=True)
    traj.add_argument('-top', dest='topology', action='store',
                      help='Path to the topology file', type=str,
                      required=False, metavar='topology', default=None)
    traj.add_argument('-first', dest='first', action='store',
                      help='First frame to analyze (counting from 0)\
                      [default: %(default)s]', type=int, required=False,
                      default=0, metavar='first_frame')
    traj.add_argument('-last', dest='last', action='store',
                      help='Last frame to analyze (counting from 0)\
                      [default: last frame]', type=int, required=False,
                      default=None, metavar='last_frame')
    traj.add_argument('-stride', dest='stride', action='store',
                      help='Stride of frames to analyze\
                      [default: %(default)s]', type=int, required=False,
                      default=1, metavar='stride')
    traj.add_argument('-sel', dest='selection', action='store',
                      help='Atom selection (MDTraj syntax)\
                      [default: %(default)s]', type=str, required=False,
                      default='all', metavar='selection')
    # Arguments: clustering parameters ----------------------------------------
    clust = parser.add_argument_group(title='Clustering options')
    clust.add_argument('-cutoff', action='store', dest='cutoff',
                       help='RMSD cutoff [default: %(default)s]',
                       type=int, required=False, default=2, metavar='k')
    clust.add_argument('-min_clust_size', action='store',
                       dest='min_clust_size',
                       help='Minimum size of returned clusters\
                           [default: %(default)s]',
                       type=int, required=False, default=2, metavar='m')
    clust.add_argument('-nclust', action='store',
                       dest='nclust',
                       help='Number of clusters to retrieve\
                       [default: %(default)s]',
                       type=int, required=False, default=np.inf, metavar='n')
    # Arguments: analysis -----------------------------------------------------
    out = parser.add_argument_group(title='Output options')
    out.add_argument('-odir', action='store', dest='outdir',
                     help='Output directory to store analysis\
                     [default: %(default)s]',
                     type=str, required=False, default='bitQT_outputs',
                     metavar='bitQT_outputs')
    user_inputs = parser.parse_args()
    return user_inputs


def is_valid_traj(traj, valid_trajs):
    """
    Check if the trajectory extension is supported by MDTraj engine.

    Parameters
    ----------
    traj : str
        Path to the trajectory file.
    valid_trajs : set
        Set of supported trajectory extensions.

    Raises
    ------
    ValueError
        If trajectory extension is not supported by MDTraj.

    Returns
    -------
    bool
        True if trajectory extension is supported.

    """
    traj_ext = traj.split('.')[-1]
    if traj_ext not in valid_trajs:
        raise ValueError('The trajectory format "{}" '.format(traj_ext) +
                         'is not available. Valid trajectory formats '
                         'are: {}'.format(valid_trajs))
    return True


def traj_needs_top(traj):
    """
    Determine if trajectory extension does not contain topological information.

    Parameters
    ----------
    traj : str
        Path to the trajectory file.

    Returns
    -------
    bool
        True if trajectory needs topological information.

    """
    traj_ext = traj.split('.')[-1]
    if traj_ext in ['h5', 'lh5', 'pdb']:
        return False
    return True


def is_valid_top(topology, valid_tops):
    """
    Check if the topology extension is supported by MDTraj engine.

    Parameters
    ----------
    topology : str
        Path to the trajectory file.
    valid_tops : set
        Set of supported topology extensions.

    Raises
    ------
    ValueError
        If topology extension is not supported by MDTraj.

    Returns
    -------
    bool
        DESCRIPTION.

    """
    try:
        top_ext = topology.split('.')[-1]
    except AttributeError:
        raise ValueError('You should pass a topology object. '
                         'Valid topology formats are: {}'.format(valid_tops))

    if top_ext not in valid_tops:
        raise ValueError('The topology format "{}"'.format(top_ext) +
                         'is not available. Valid topology formats'
                         'are: {}'.format(valid_tops))
    return True


def load_raw_traj(traj, valid_trajs, topology=None):
    """
    Load the whole trajectory without any modifications.

    Parameters
    ----------
    traj : str
        Path to the trajectory file.
    valid_trajs : set
        Set of supported trajectory extensions.
    topology : str, optional
        Path to the trajectory file. The default is None.

    Returns
    -------
    mdtraj.Trajectory
        Raw trajectory.

    """
    if is_valid_traj(traj, valid_trajs) and traj_needs_top(traj):
        if is_valid_top(topology, valid_tops):
            return md.load(traj, top=topology)

    if is_valid_traj(traj, valid_trajs) and not traj_needs_top(traj):
        return md.load(traj)


def shrink_traj_selection(traj, selection):
    """
    Select a subset of atoms from the trajectory.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        Trajectory object to which selection will be applied.
    selection : str
        Any MDTraj valid selection.

    Raises
    ------
    ValueError
        If specified selection is not valid.
        If specified selection corresponds to no atoms.

    Returns
    -------
    traj : mdtraj.Trajectory
        Trajectory containing the subset of specified atoms.

    """
    if selection != 'all':
        try:
            sel_indx = traj.topology.select(selection)
        except Exception:
            raise ValueError('Specified selection is invalid')
        if sel_indx.size == 0:
            raise ValueError('Specified selection corresponds to no atoms')
        traj = traj.atom_slice(sel_indx, inplace=True)
    return traj


def shrink_traj_range(first, last, stride, traj):
    """
    Select a subset of frames from the trajectory.

    Parameters
    ----------
    first : int
        First frame to consider (0-based indexing).
    last : TYPE
        Last frame to consider (0-based indexing).
    stride : TYPE
        Stride (step).
    traj : mdtraj.Trajectory
        Trajectory object to which slicing will be applied.

    Raises
    ------
    ValueError
        If first, last or stride are falling out of their valid ranges.

    Returns
    -------
    mdtraj.Trajectory
        Trajectory containing the subset of specified frames.

    """
    # Calculate range of available intervals ----------------------------------
    N = traj.n_frames
    first_range = range(0, N - 1)
    last_range = range(first + 1, N)
    try:
        delta = last - first
    except TypeError:
        delta = N - first
    stride_range = range(1, delta)
    # Raising if violations ---------------------------------------------------
    if first not in first_range:
        raise ValueError('"first" parameter should be in the interval [{},{}]'
                         .format(first_range.start, first_range.stop))
    if last and (last not in last_range):
        raise ValueError('"last" parameter should be in the interval [{},{}]'
                         .format(last_range.start, last_range.stop))
    if stride not in stride_range:
        raise ValueError('"stride" parameter should be in the interval [{},{}]'
                         .format(stride_range.start, stride_range.stop))
    # Slicing trajectory ------------------------------------------------------
    sliced = slice(first, last, stride)
    if sliced not in [slice(0, N, 1), slice(0, None, 1)]:
        return traj.slice(sliced)
    return traj


def calc_rmsd_matrix(trajectory, args):
    """
    Calculate optimal RMSD binary-encoded square matrix using MDTraj. Pairwise
    similarity is saved in RAM as bits (dict of bitarrays), not floats.

    Parameters
    ----------
    trajectory : mdtraj.Trajectory
        MDTraj trajectory object.
    args : argparse.Namespace
        user input parameters parsed by argparse (CLI).

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
        seen.add(max_node)
        xcolor += 1
        not_neighbors = ~ matrix[max_node]
        not_colored = ~colored
        candidates = not_neighbors & not_colored
        # Nodes passing conditions (not-neighb, not-colored, not-neighb) ------
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
        # Deliver a color class to passed nodes -------------------------------
        colors[passed] = xcolor
        colored = ba()
        colored.pack(colors.astype(np.bool).tobytes())
        if colored.count(0) == 0:
            break
    return colors


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
        if counter >= COUNT:
            break
    ar = np.nonzero(np.unpackbits(init_cascade).astype(np.bool))[0]
    return init_cascade, ar


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


def top_has_coords(topology):
    """
    Check if topology has cartesian coordinates information.

    Parameters
    ----------
    topology : str
        Path to the topology file.

    Returns
    -------
    int
        Number of cartesian frames if topology contains cartesians.
        False otherwise.
    """
    try:
        tt = md.load(topology)
    except OSError:
        return False
    return tt.xyz.shape[0]


def to_VMD(outdir, topology, first, N1, last, stride, final_array):
    """
    Create a .log file for visualization of clusters in VMD through a
    third-party plugin.

    Parameters
    ----------
    outdir : str
        Path where to create the VMD visualization .log.
    topology : str
        Path to the topology file.
    first : int
        First frame to consider (0-based indexing).
    N1 : int
        default value when last == None.
    last : TYPE
        Last frame to consider (0-based indexing).
    stride : TYPE
        Stride (step).
    final_array : numpy.ndarray
        Final labeling of the selected clusters ordered by size (descending).

    Returns
    -------
    logname : str
        Log file to be used with VMD.
    """
    basename = os.path.basename(topology).split('.')[0]
    logname = os.path.join(outdir, '{}.log'.format(basename))
    vmd_offset = top_has_coords(topology)
    start = first
    if not last:
        stop = N1
    else:
        stop = last
    slice_frames = np.arange(start, stop, stride, dtype=np.int32)
    nmr_offset = 1
    with open(logname, 'wt') as clq:
        for num in np.unique(final_array):
            if num != 0:
                clq.write('{}:\n'.format(num))
                cframes = np.where(final_array == num)[0]
                if vmd_offset:
                    real_frames = slice_frames[cframes] + nmr_offset + vmd_offset
                else:
                    real_frames = slice_frames[cframes] + nmr_offset
                str_frames = [str(x) for x in real_frames]
                members = ' '.join(str_frames)
                clq.write('Members: ' + members + '\n\n')
        if 0 in np.unique(final_array):
            clq.write('{}:\n'.format(0))
            cframes = np.where(final_array == 0)[0]
            if vmd_offset:
                real_frames = slice_frames[cframes] + nmr_offset + vmd_offset
            else:
                real_frames = slice_frames[cframes] + nmr_offset
            str_frames = [str(x) for x in real_frames]
            members = ' '.join(str_frames)
            clq.write('Members: ' + members + '\n\n')
    return logname


def get_frames_stats(N1, first, last, stride, clusters, outdir):
    """
    Get "frames_statistics.txt" containing frameID, clusterID.

    Parameters
    ----------
    N1 : int
        default value when last == None.
    first : int
        First frame to consider (0-based indexing).

    last : TYPE
        Last frame to consider (0-based indexing).
    stride : TYPE
        Stride (step).
    clusters : numpy.ndarray
        array of clusters ID.
    outdir : str
        Path where to create the VMD visualization .log.

    Returns
    -------
    frames_df : pandas.DataFrame
        dataframe with frames_statistics info.
    """
    start = first
    if not last:
        stop = N1
    else:
        stop = last
    slice_frames = np.arange(start, stop, stride, dtype=np.int32)
    frames_df = pd.DataFrame(columns=['frame', 'cluster_id'])
    frames_df['frame'] = range(N1)
    frames_df['cluster_id'].loc[slice_frames] = clusters
    with open(os.path.join(outdir, 'frames_statistics.txt'), 'wt') as on:
        frames_df.to_string(buf=on, index=False)
    return frames_df


def get_cluster_stats(clusters, outdir):
    """
    Get "cluster_statistics.txt" containing clusterID, cluster_size, and
    cluster percentage from trajectory.

    Parameters
    ----------
    clusters : numpy.ndarray
        array of clusters ID.
    outdir : str
        Path where to create the VMD visualization .log.

    Returns
    -------
    clusters_df : pandas.DataFrame
        dataframe with cluster_statistics info.
    """
    clusters_df = pd.DataFrame(columns=['cluster_id', 'size', 'percent'])
    clusters_df['cluster_id'] = list(range(0, clusters.max() + 1))
    sizes = []
    for x in clusters_df.cluster_id:
        sizes.append(len(np.where(clusters == x)[0]))
    clusters_df['size'] = sizes

    sum_ = clusters_df['size'].sum()
    percents = [round(x / sum_ * 100, 4) for x in clusters_df['size']]
    clusters_df['percent'] = percents

    with open(os.path.join(outdir, 'cluster_statistics.txt'), 'wt') as on:
        clusters_df.to_string(buf=on, index=False)
    return clusters_df


if __name__ == '__main__':
    # >>>> Debugging <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # import argparse
    # args = argparse.Namespace()
    # args.trajectory = '../trajs/aligned_tau.dcd'
    # args.topology = '../trajs/aligned_tau.pdb'
    # args.nclust = np.inf
    # args.min_clust_size = 2
    # args.first = 1000
    # args.last = 6000
    # args.stride = 3
    # args.selection = 'all'
    # args.cutoff = 4
    # args.outdir = 'bitQT_outputs'
    # =========================================================================
    # 1. Creating binary matrix (adjacency list)
    # =========================================================================
    # ++++ Get adjacency matrix of trajectory as list of bitarrays ++++++++++++
    args = parse_arguments()

    try:
        os.makedirs(args.outdir)
    except FileExistsError:
        raise Exception('{} directory already exists.'.format(args.outdir) +
                        'Please specify another location or rename it.')

    trajectory = load_raw_traj(args.trajectory, valid_trajs, args.topology)
    trajectory = shrink_traj_selection(trajectory, args.selection)
    N1 = trajectory.n_frames
    trajectory = shrink_traj_range(args.first, args.last, args.stride, trajectory)
    trajectory.center_coordinates()
    matrix = calc_rmsd_matrix(trajectory, args)
    # ++++ Tracking clust/uNCLUSTERed bits to avoid re-computations +++++++++++
    N = len(matrix[0])
    m = len(matrix)
    unclust_bit = ba(N)
    unclust_bit.setall(1)
    clustered_bit = unclust_bit.copy()
    clustered_bit.setall(0)
    zeros = np.zeros(N, dtype=np.int32)
    # ++++ Save clusters in an array (1 .. N) +++++++++++++++++++++++++++++++++
    clusters_array = np.zeros(N, dtype=np.int32)
    NCLUSTER = 0
    clustered = set()
    nmembers = []
    # ++++ Coloring ordered vertices (1 .. N) +++++++++++++++++++++++++++++++++
    degrees = calc_matrix_degrees(unclust_bit, matrix)
    ordered_by_degs = degrees.argsort()[::-1]
    colors = colour_matrix(ordered_by_degs, matrix)
    # colors[np.frombuffer(clustered_bit.unpack(), dtype=np.bool)] = 0

    # =========================================================================
    # 2. Main algorithm: BitQT !
    # =========================================================================
    while True:
        NCLUSTER += 1
        # ++++ Find a big clique early ++++++++++++++++++++++++++++++++++++++++
        big_node = degrees.argmax()
        bit_clique, big_clique = do_bit_cascade(big_node, degrees, colors,
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
                bit_clique, clique = do_bit_cascade(node, degrees, colors,
                                                    matrix, big_clique_size)
                CLIQUE_SIZE = len(clique)
            except TypeError:
                CLIQUE_SIZE = 0
            # ++++ Cumulative update only if biggers candidates are found +++++
            if CLIQUE_SIZE > big_clique_size:
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

        if (big_clique_size < args.min_clust_size) or (NCLUSTER == args.nclust):
            break

        # ++++ Save new cluster & update NCLUSTER +++++++++++++++++++++++++++++
        clusters_array[big_clique] = NCLUSTER
        # ++++ Update (un)clustered_bit +++++++++++++++++++++++++++++++++++++++
        clustered.update(big_clique)
        clustered_bit = set_to_bitarray(clustered, N)
        unclust_bit = ~clustered_bit
        # ++++ Hard erasing of clustered frames from matrix +++++++++++++++++++
        degrees = zeros.copy()
        for x in unclust_bit[:m].itersearch(ba('1')):
            degrees[x] = matrix[x].count()
            if bu.count_and(matrix[x], clustered_bit):
                matrix[x] &= (matrix[x] ^ clustered_bit)

    # =========================================================================
    # 3. Output
    # =========================================================================
    # saving pickle for api debugging tests
    outname = os.path.basename(args.topology).split('.')[0]
    pickle_to_file(clusters_array, os.path.join(args.outdir,
                                                '{}.pick'.format(outname)))
    # saving VMD visualization script
    to_VMD(args.outdir, args.topology, args.first, args.last, N1, args.stride,
           clusters_array[:m])
    # saving clustering info  files
    frames_stats = get_frames_stats(N1, args.first, args.last, args.stride,
                                    clusters_array[:m], args.outdir)
    cluster_stats = get_cluster_stats(clusters_array[:m], args.outdir)
    print('\n\nNormal Termination of BitQT :)')
