#!/usr/bin/env python2

from __future__ import division
import numpy as np, scipy as sp
from . import data

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

def make_pair_restraints(reference, max_distance=10.0, noise_weight=2.0):
    """
    Return a set of restraints based on the distance between all pairs of 
    particles in each member of the ensemble.  These restraints are meant to 
    mimic Hi-C data.  The restraints are extracted from an ensemble of 
    structures to mimic the fact that Hi-C data represents many structures.
    """

    from scipy.spatial.distance import pdist, squareform

    # Calculate the pairwise distance between every bead.
    restraints = squareform(pdist(reference))

    # A negative distance means the beads are at least that far apart.
    restraints[restraints > max_distance] = -max_distance

    # Add a user-controlled amount of noise to the data.
    restraints += noise_weight * np.random.normal(0.0, 1.0, restraints.shape)

    return restraints

def find_average_bond_length(coords):
    diff = coords[:-1] - coords[1:]
    diff_sqr = diff**2
    lengths_sqr = np.sum(diff_sqr, axis=1)
    lengths = np.sqrt(lengths_sqr)
    return np.mean(lengths) / 2

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

