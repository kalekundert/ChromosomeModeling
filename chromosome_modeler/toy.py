#!/usr/bin/env python2

from __future__ import division
import numpy as np, scipy as sp
if __name__ != '__main__': from . import data
else: import data

def make_hilbert_reference(num_particles):
    """
    Create a reference model by placing the given number of particles along the 
    vertices of a Hilbert curve.  This makes a dense but knot-free arrangement
    particles in 3-dimensional space.
    
    Parameters
    ----------
    num_particles : int

    Return
    ------
    toy_coords : coords
    """

    import hilbert

    num_particles = parse_num_particles(num_particles)
    coords = np.zeros((num_particles, 3))

    for index in range(num_particles):
        coords[index] = hilbert.int_to_hilbert(index, 3)

    return coords

def make_nagano_reference():
    from h5py import File
    from numpy import array, concatenate
    from glob import glob
    from os.path import join, dirname

    ensemble = []
    data = join(dirname(__file__), 'nagano_models', '*.hdf5')

    for path in glob(data):
        file = File(path)
        coords = file['/structures/0/coords/X']
        ensemble.append(coords)

    return concatenate(ensemble)

def make_interpolated_trajectory():
    raise NotImplementedError

def make_md_trajectory():
    raise NotImplementedError

def make_xyz_restraints(reference, num_restraints):
    """
    Restraint an arbitrary set of well-spaced particles to their reference 
    coordinates.  This mimics the data from a microscopy experiment that has 
    labeled the given number of known points along a chromosome.
    """
    
    import collections, random; random.seed(0)

    data.require_coords(reference)
    num_particles = len(reference)
    num_restraints = parse_num_restraints(num_restraints, num_particles)

    restraints = data.xyz_restraints()
    restraint_stride = num_particles // num_restraints
    possible_indices = range(0, num_particles, restraint_stride)
    shuffled_indices = sorted(possible_indices, key=lambda x: random.random())
    restraint_indices = shuffled_indices[:num_restraints]

    for index in restraint_indices:
        restraints[index] = reference[index]

    return restraints

def make_pair_restraints(ensemble, num_restraints, cutoff_distance=None):
    """
    Pick an arbitrary set of neighboring particles to use as pair restraints.  
    This mimics population Hi-C data, which identifies chromatin contacts 
    throughout the genome.  Because the same set of neighbors may be picked 
    more than once, the return value is a histogram.
    """

    # We will use the random module with a fixed seed to generate arbitrary 
    # (i.e. irregular) set of pair restraints.

    import random; random.seed(0)
    from scipy.spatial.distance import cdist

    restraints = data.pair_restraints()
    data.require_ensemble(ensemble)

    if cutoff_distance is None:
        cutoff_distance = 1.1 * np.sqrt(2)

    for i in range(num_restraints):

        # Pick an arbitrary set of coords from the ensemble.

        coords = ensemble[random.randrange(len(ensemble))]

        # Pick an arbitrary particle from those coords.

        index_1 = random.randrange(len(coords))

        # Find all the neighboring particles that are within a certain distance 
        # but not adjacent in sequence.  Arbitrarily pick one to restrain.
        
        pair_distances = cdist([coords[index_1]], coords)[0]
        pair_distances[max(index_1-1,0):index_1+2] = np.inf
        neighbors = np.where(pair_distances < cutoff_distance)[0]

        if len(neighbors):
            index_2 = random.choice(neighbors)
            restraints[index_1, index_2] += 1

    return restraints

def find_average_bond_length(coords):
    diff = coords[:-1] - coords[1:]
    diff_sqr = diff**2
    lengths_sqr = np.sum(diff_sqr, axis=1)
    lengths = np.sqrt(sum_diff_sqr)
    return np.mean(lengths)

def superimpose_model(model, reference):
    from numpy import eye, sign, mean, matrix
    from numpy.linalg import svd, det

    # Reorient the matrices, because this algorithm expects them to be 3xN.  
    # At the same time, convert them into matrices to make some of the 
    # following operations more convenient.

    raw_target = matrix(reference.T)
    raw_mobile = matrix(model.T)

    # Translate both structures to the origin, so that they only differ by a 
    # simple rotation.  Save the translation vector that can be used to 
    # position the mobile structure on top of the target one after the 
    # rotation has been applied.

    target = raw_target - mean(raw_target, 1)
    mobile = raw_mobile - mean(raw_mobile, 1)

    translation = mean(raw_target, 1)

    # Calculate the rotation that transforms the mobile structure into the 
    # best possible alignment with the target structure.  I don't really 
    # understand how this works, although the website I took it from had a 
    # very complete description. 

    C = mobile * target.T
    V, S, W = svd(C)

    I = eye(3)
    I[2,2] = sign(det(C))

    rotation = W.T * I * V.T

    # Save the aligned coordinates.
    
    return np.array(rotation * mobile + translation).T

def evaluate_model(model, reference, mode='rmsd'):
    """
    Score the similarity between a model and its true reference structure.  
    Different metrics can be used to calculate scores.

    Parameters
    ----------
    model : array
        An Nx3 array containing all the coordinates for the model.

    reference : array
        An Nx3 array containing all the coordinates for the reference structure 
        that the model was trying to recapitulate.

    mode : string
        The score metric to use.  The default is 'rmsd', which calculates the 
        root-mean-squared distance between each particle in the two structures.  
        The alternative is 'pairwise', which calculates the average absolute 
        difference between the pairwise distance matrices for the two 
        structures.

    Returns
    -------
    score : float
        A scalar value reflecting how similar the two sets of coordinates are.
    """

    from scipy.spatial.distance import pdist

    data.require_coords(model);
    data.require_coords(reference)
    assert model.shape == reference.shape

    if mode == 'rmsd':
        diff = model - reference
        return np.sqrt(np.sum(diff**2) / len(model))

    elif mode == 'pairwise':
        model_dists = pdist(model)
        reference_dists = pdist(reference)
        dist_errors = np.abs(model_dists - reference_dists) / reference_dists
        num_dists = dist_errors.shape[0]
        return sum(dist_errors) / num_dists

    else:
        raise ValueError("unknown mode: '{}'".format(mode))


def parse_num_particles(num_particles):
    return int(num_particles)

def parse_num_restraints(num_restraints, num_particles):
    if isinstance(num_restraints, str) and num_restraints.endswith('%'):
        num_restraints = int(0.01 * float(num_restraints[:-1]) * num_particles)
    else:
        num_restraints = int(num_restraints)

    return min(num_restraints, num_particles)


if __name__ == '__main__':

    ## Test make_hilbert_reference()

    hilbert_coords = make_hilbert_reference(8)

    assert tuple(hilbert_coords[0]) == (0, 0, 0)
    assert tuple(hilbert_coords[1]) == (1, 0, 0)
    assert tuple(hilbert_coords[2]) == (1, 0, 1)
    assert tuple(hilbert_coords[3]) == (0, 0, 1)
    assert tuple(hilbert_coords[4]) == (0, 1, 1)
    assert tuple(hilbert_coords[5]) == (1, 1, 1)
    assert tuple(hilbert_coords[6]) == (1, 1, 0)
    assert tuple(hilbert_coords[7]) == (0, 1, 0)

    ## Test make_nagano_reference()

    nagano_ensemble = make_nagano_reference()
    data.require_ensemble(nagano_ensemble)

    ## Test make_xyz_restraints()

    xyz_restraints = make_xyz_restraints(hilbert_coords, 2)
    assert all(np.abs(np.diff(xyz_restraints.keys())) >= 8 // 2)

    ## Test make_pair_restraints()

    hilbert_ensemble = np.array([hilbert_coords])
    pair_restraints = make_pair_restraints(hilbert_ensemble, 10)
    assert all(np.abs(np.diff(pair_restraints.keys())) >= 2)

    ## Test evaluate_model()

    reference = np.array([
        [1.0, 0.0, 1.0],
        [2.0, 0.0, 1.0],
        [2.0, 0.0, 2.0],
        [1.0, 0.0, 2.0],
    ])

    model = np.array([
        [0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [3.0, 0.0, 3.0],
        [0.0, 0.0, 3.0],
    ])

    score = evaluate_model(model, reference)

    assert score == np.sqrt(2)


    print 'All tests passed!'

