#!/usr/bin/env python2

"""
Generate toy data showing that combined Hi-C and imaging data gives better 
models than either type of data alone.  This is a fairly intuitive result, 
because the two types of data are quite orthogonal.

Usage:
    run_md.py (xyz|pair|both) [-h]
"""

import sys, docopt
import chromosome_modeler as cmod

args = docopt.docopt(__doc__)
num_particles = 512
reference = cmod.toy.make_hilbert_reference(num_particles)
protocol = cmod.define_logarithmic_protocol(3e4, 1e2, 1e-3, 100)
xyz_restraints, pair_restraints, initial_coords = None, None, 'random'

def output_path(ext):
    if args['xyz']: base = 'xyz'
    if args['pair']: base = 'pair'
    if args['both']: base = 'both'
    return base + ext


# Setup the toy system.

if args['xyz'] or args['both']:
    initial_coords = 'interpolated'
    xyz_restraints = cmod.toy.make_xyz_restraints(
            reference, num_restraints=10)

if args['pair'] or args['both']:
    pair_restraints = cmod.toy.make_pair_restraints(
            reference, max_distance=2.0, signal_to_noise=2.0)

system = cmod.define_system(
        num_particles,
        xyz_restraints=xyz_restraints,
        pair_restraints=pair_restraints,
        initial_coords=initial_coords,
)

# Build a model to fit the toy system.

if args['xyz']:
    model = cmod.build_interpolated_model(
            num_particles,
            xyz_restraints,
    )
else:
    model = cmod.run_molecular_dynamics(
            system,
            protocol=protocol,
            frames=100,
            progress_bar=True,
            movie_path=output_path('.pym'),
    )[-1]

if args['pair']:
    cmod.export_to_pdb(output_path('_no_super.pdb'), model, reference)
    model = cmod.superimpose_model(model, reference)

# Report how good the model is.

rmsd = cmod.toy.evaluate_model(model, reference)
cmod.print_restraint_satisfaction(system)
print 'Reference RMSD:', rmsd
cmod.export_to_pdb(output_path('.pdb'), model, reference, xyz_restraints)
with open(output_path('.rmsd'), 'a') as file:
    file.write('{}\n'.format(rmsd))

