#!/usr/bin/env python2

"""\
A collection of functions capable of building both chromosome models and 
trajectories of chromosome models from fluorescence microscopy data.  The basic 
modeling approach is to convert the experimental data into a restraint-based 
score function and to search for low-scoring models.  In principle, any data 
can be used as a restraint, so it is easy to add new sources of data to this 
framework.  

Data Types
==========
The functions in this framework operate on a few different data types, which 
are documented below.  Understanding these types is essential to understanding 
what the code is supposed to do.

coords : num_particles x 3 array
    Store X, Y, and Z coordinates for each particle in a model.

ensemble/trajectory : num_models x num_particles x 3 array
    Store coordinates for every member of an ensemble.  Ensembles can be built 
    using either MD or MC and typically attempt to fit some experimental data.

ensemble_trajectory : list of num_frames ensembles
    Store ensembles generated for each frame of a trajectory.  An ensemble 
    trajectory can be sampled to get at a set of reasonable trajectories.

xyz_restraints : {int : array}
    Map particle indices (int) to 3D coordinates (array) for each position that 
    was observed in a microscopy experiment.

pair_restraints : {set(int) : int}
    Map pairs of particle indices (set) to the count (int) of how many times 
    those indices were observed to be cross-linked in a Hi-C experiment.
"""

from __future__ import division

from . import data, utils
import IMP, IMP.core, IMP.atom, IMP.container
from IMP.algebra import Vector3D
import sys, numpy as np, scipy as sp

def load_microscopy_restraints():
    raise NotImplementedError

def load_hic_restraints():
    raise NotImplementedError

def define_system(num_particles, **kwargs):
    """
    Return a data structure representing the entire system to be simulated.  
    This system includes restraints and particles with radii, masses, and 
    initial positions.  Other functions are provided to search this system for 
    low-scoring conformations.
    """

    ## Define simulation parameters

    recursive_kwargs = kwargs.copy()
    recursive_kwargs.pop('initial_coords', None)

    xyz_restraints = kwargs.pop('xyz_restraints', None)
    pair_restraints = kwargs.pop('pair_restraints', None)
    particle_radius = float(kwargs.pop('radius', 0.5))
    particle_mass = float(kwargs.pop('mass', 1.0))
    initial_coords = kwargs.pop('initial_coords', 'random')
    nonbonded = kwargs.pop('nonbonded', 'lennard-jones')
    k_bonded = float(kwargs.pop('k_bonded', 1.0))
    k_nonbonded = float(kwargs.pop('k_nonbonded', 1.0))
    k_xyz = float(kwargs.pop('k_xyz', 1.0))
    k_pair = float(kwargs.pop('k_pair', 1.0))
    lj_well_depth = float(kwargs.pop('lj_depth', 1e-2))
    nblist_cutoff = float(kwargs.pop('nblist_cutoff', 3.0 * particle_radius))
    imp_logging = kwargs.pop('imp_logging', False)

    utils.check_all_kwargs_used(kwargs)

    ## Create the particles

    if imp_logging:
        if imp_logging is True: imp_logging = IMP.base.TERSE
        IMP.base.set_log_level(imp_logging)

    system = IMP.kernel.Model()
    particle_list = IMP.core.create_xyzr_particles(
            system, num_particles, particle_radius)

    for p in particle_list:
        IMP.atom.Mass.setup_particle(p, particle_mass)
        IMP.atom.LennardJones.setup_particle(p, lj_well_depth)

    if initial_coords == 'random':
        pass

    elif initial_coords == 'circle':
        from math import pi, sin, cos
        num_particles = len(particle_list)
        for i, p in enumerate(particle_list):
            r = num_particles * particle_radius / pi
            th = 2 * pi * i / num_particles
            xyz = r * cos(th), r * sin(th), 0
            IMP.core.XYZ(p).set_coordinates(xyz)

    elif initial_coords == 'interpolated':
        coords = build_interpolated_model(num_particles, **recursive_kwargs)
        for xyz, p in zip(coords, particle_list):
            IMP.core.XYZ(p).set_coordinates(xyz)

    elif initial_coords == 'minimized':
        coords = build_minimized_model(num_particles, **recursive_kwargs)
        for xyz, p in zip(coords, particle_list):
            IMP.core.XYZ(p).set_coordinates(xyz)

    elif data.are_coords(initial_coords):
        for xyz, p in zip(initial_coords, particle_list):
            IMP.core.XYZ(p).set_coordinates(xyz)

    else:
        raise ValueError("unknown config '{}'".format(config))

    particles = IMP.container.ListSingletonContainer(particle_list)

    ## Bonded restraints

    hdps = IMP.core.HarmonicDistancePairScore(1.5 * particle_radius, k_bonded)
    bonded = IMP.container.ExclusiveConsecutivePairContainer(particle_list)
    res = IMP.container.PairsRestraint(hdps, bonded, "Bonded Restraints")
    system.add_restraint(res)

    ## Non-bonded restraints

    if nonbonded == 'lennard-jones':
        nbl = IMP.container.ClosePairContainer(particles, nblist_cutoff)
        nbl.add_pair_filter(IMP.container.ExclusiveConsecutivePairFilter())
        fs = IMP.atom.ForceSwitch(0.8 * nblist_cutoff, nblist_cutoff)
        ljs = IMP.atom.LennardJonesPairScore(fs)
        res = IMP.container.PairsRestraint(ljs, nbl, "Lennard-Jones Restraint")

    elif nonbonded == 'excluded-volume':
        res = IMP.core.ExcludedVolumeRestraint(particles, k_nonbonded)
        res.set_name("Excluded Volume Restraint")

    else:
        raise ValueError("unknown non-bonded restraint: '{}'".format(nonbonded))

    system.add_restraint(res)

    ## Coordinate restraints (from microscopy data)

    if xyz_restraints is not None:
        res_set = IMP.kernel.RestraintSet(system, "Coordinate Restraints")

        for i, coord in xyz_restraints.items():
            particle = particle_list[i]
            func = IMP.core.Harmonic(0, k_xyz)
            score = IMP.core.DistanceToSingletonScore(func, coord)
            res = IMP.core.SingletonRestraint(score, particle)
            res_set.add_restraint(res)

        system.add_restraint(res_set)

    ## Pair restraints (from Hi-C data)

    if pair_restraints is not None:
        hdps = IMP.core.HarmonicDistancePairScore(1.5 * particle_radius, k_pair)
        pairs = IMP.container.ListPairContainer(
             [(particle_list[i], particle_list[j]) for i,j in pair_restraints])
        res = IMP.container.PairsRestraint(hdps, pairs, "Pair Restraints")
        system.add_restraint(res)


    return system

def run_minimization(system):
    """
    Use a gradient-minimization algorithm to search for a low-scoring 
    conformation in the given system.  Because gradient-minimization can only 
    find local minima, the coordinates it returns may not be that good.
    """

    minimizer = IMP.core.ConjugateGradients(system)
    minimizer.optimize(1000)
    return convert_system_to_coords(system)

def run_molecular_dynamics(system, **kwargs):
    """
    Use molecular dynamics to search for a low-scoring conformation of the 
    given system.  The dynamics of this system are not physical, but using an 
    MD integrator allows the system to overcome energy barriers between 
    low-scoring states.  The coordinates returned by this function depend 
    heavily on the parameters given, so expect to have to optimize a bit.
    """

    # The in(float(...)) idiom makes it easier for values like '1e5' to be 
    # specified on the command line and passed into this function directly.

    steps = int(float(kwargs.pop('steps', 1e4)))
    movie_path = kwargs.pop('movie_path', None)
    frames = int(float(kwargs.pop('frames', 100 if movie_path else 1)))
    temperature = float(kwargs.pop('temperature', 1e-2))
    progress_bar = bool(kwargs.pop('progress_bar', bool(movie_path)))

    utils.check_all_kwargs_used(kwargs)

    # Optimize using MD.

    num_particles = system.get_number_of_particles()
    trajectory = np.zeros((frames, num_particles, 3))

    md = IMP.atom.MolecularDynamics(system)
    md.add_optimizer_state(
            IMP.atom.VelocityScalingOptimizerState(
                system, system.get_particles(), temperature))

    if movie_path:
        writer = IMP.display.PymolWriter(movie_path)
        convert_system_to_movie_frame(writer, system)

    for frame in range(frames):
        md.optimize(steps // frames)
        trajectory[frame] = convert_system_to_coords(system)

        if movie_path:
            convert_system_to_movie_frame(writer, system)

        if progress_bar:
            sys.stdout.write('\r[{}/{}]'.format(frame+1, frames))
            sys.stdout.flush()

    if progress_bar:
        print

    return trajectory

def run_monte_carlo(system, **kwargs):
    pass

def build_interpolated_model(num_particles, xyz_restraints):
    """
    Build a model by simple linearly interpolating between the particles with 
    known positions.  The purpose of this model is to show whether or not the 
    fancier methods are really doing anything.  I might also want to use a 
    naive model like this to initialize one of the fancier methods.
    """

    from .utils import pairwise

    model = np.zeros((num_particles, 3))

    # Place the restrained particles themselves.

    for i, coord in xyz_restraints.items():
        model[i] = coord

    # Place all the particles between any two known coordinates.

    restraint_indices = sorted(xyz_restraints.keys())

    for i, j in pairwise(restraint_indices):
        start = xyz_restraints[i]
        end = xyz_restraints[j]
        steps = abs(j - i) - 1

        if steps == 0:
            continue

        disp_vector = end - start
        step_vector = disp_vector / (steps + 1)

        for k in range(i + 1, j):
            model[k] = model[k-1] + step_vector

    # Place the particles before the first restraint.

    start = restraint_indices[0] - 1

    for i in range(start, -1, -1):
        model[i] = 2 * model[i+1] - model[i+2]

    # Place the particles after the last restraint.
    
    start = restraint_indices[-1] + 1

    for i in range(start, num_particles):
        model[i] = 2 * model[i-1] - model[i-2]

    return model

def build_minimized_model(num_particles, **kwargs):
    system = define_system(num_particles, **kwargs)
    return run_minimization(system)

def build_md_model(num_particles, **kwargs):
    return build_md_trajectory(num_particles, **kwargs)[0]

def build_md_trajectory(num_particles, **kwargs):
    system = define_system(num_particles, check_args=False, **kwargs)
    return run_molecular_dynamics(system, check_args=False, **kwargs)

def build_md_trajectories(
        num_particles, xyz_restraint_traj=None, pair_restraint_traj=None,
        **kwargs):

    assert xyz_restraint_traj or pair_restraint_traj

    trajectories = []
    restraints = zip(
            xyz_restraint_traj or itertools.repeat(None),
            pair_restraint_traj or itertools.repeat(None))

    for xyz_restraints, pair_restraints in restraints:
        trajectory = build_md_trajectory(
                num_particles,
                xyz_restraints=xyz_restraints,
                pair_restraints=pair_restraints,
                **kwargs)
        trajectories.append(trajectory)

    return trajectories

def build_initial_trajectory():
    raise NotImplementedError

def search_for_best_trajectory():
    raise NotImplementedError


def convert_system_to_coords(system):
    particles = system.get_particles()
    num_particles = len(particles)

    coords = np.zeros((num_particles, 3))
    for i, p in enumerate(particles):
        coords[i] = IMP.core.XYZR(p).get_coordinates()

    return coords

def convert_system_to_movie_frame(writer, system):
    particle_list = system.get_particles()
    particles = IMP.container.ListSingletonContainer(particle_list)
    neighbors = IMP.container.ConsecutivePairContainer(particle_list)

    frame = writer.get_frame() + 1
    writer.set_frame(frame)

    g = IMP.core.XYZRsGeometry(particles)
    g.set_name('particles')
    writer.add_geometry(g)

    g = IMP.core.EdgePairsGeometry(neighbors)
    g.set_name('bonds')
    writer.add_geometry(g)

def print_restraint_satisfaction(system):
    for restraint in system.get_root_restraint_set().get_restraints():
        print '{}: {}'.format(restraint.get_name(), restraint.get_last_score())

def export_to_pdb(pdb_path, model, reference, xyz_restraints):
    """
    Export the whole scene to a pdb file.  The model and reference structures 
    will be put in separate chains.  Another chain will be created containing 
    just the atoms that were restrained, so that these can be highlighted 
    separately.
    """

    import string

    atom_lines = []
    conect_lines = []

    def coords_to_pdb(coords, conect=True):                     # (no fold)
        atom_attrs = {
                'atom_name': 'C', 'residue_name': 'CHR'}
        atom_template = (
                'HETATM{atom_id:5d}  '
                '{atom_name:3s} {residue_name:3s} '
                '{chain_name:1s}{residue_id:4d}    '
                '{coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}  '
                '1.00  1.00\n')
        
        for index, coord in enumerate(coords):
            atom_id = coords_to_pdb.atom_id
            chain_id = coords_to_pdb.chain_id

            atom_attrs['atom_id'] = atom_id
            atom_attrs['residue_id'] = index + 1
            atom_attrs['chain_name'] = string.uppercase[chain_id]
            atom_attrs['coords'] = coord
            atom_lines.append(atom_template.format(**atom_attrs))

            if conect:
                neighbors = []
                if index > 0: neighbors.append(atom_id - 1)
                if index < len(coords) - 1: neighbors.append(atom_id + 1)
                conect_line = 'CONECT {:4d} '.format(atom_id)
                conect_line += ' '.join('{:4d}'.format(x) for x in neighbors)
                conect_lines.append(conect_line + '\n')

            coords_to_pdb.atom_id += 1
        coords_to_pdb.chain_id += 1
        
    coords_to_pdb.atom_id = 1
    coords_to_pdb.chain_id = 0

    coords_to_pdb(reference)
    coords_to_pdb(xyz_restraints.values(), conect=False)
    coords_to_pdb(model)
    
    with open(pdb_path, 'w') as file:
        file.writelines(atom_lines + conect_lines)


