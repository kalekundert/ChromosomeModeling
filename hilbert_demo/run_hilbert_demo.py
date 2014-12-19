#!/usr/bin/env python2

import sys, time, json; sys.path.append('.')
import chromosome_modeler as cmod
from libraries import utils

params = utils.read_params()
num_particles = params['num_particles']
num_xyz_restraints = params['num_xyz_restraints']
num_pair_restraints = params['num_pair_restraints']

models = {}
reference = cmod.toy.make_hilbert_reference(num_particles)
xyz_restraints = cmod.toy.make_xyz_restraints(reference, num_xyz_restraints)
pair_restraints = cmod.toy.make_pair_restraints(reference, num_pair_restraints)
annealing_protocol = cmod.define_protocol(3e4, [10.0, 1.0, 0.1])

print 'num_particles, num_xyz_restraints, num_pair_restraints = {}, {}, {}'.format(
        num_particles, num_xyz_restraints, num_pair_restraints)

json_path = 'jsons/{:02}_{:03}.json'.format(
        num_xyz_restraints, num_pair_restraints)
pdb_path_template = 'pdbs/{:02}_{:03}_{{}}.pdb'.format(
        num_xyz_restraints, num_pair_restraints)
movie_path_template = 'movies/{:02}_{:03}_{{}}.pym'.format(
        num_xyz_restraints, num_pair_restraints)

now = lambda: time.strftime("%X")

print ' ', now(), 'iterp...'
models['interp'] = cmod.build_interpolated_model(
        num_particles,
        xyz_restraints,
)
print ' ', now(), 'xyz_min...'
models['xyz_min'] = cmod.build_minimized_model(
        num_particles,
        xyz_restraints=xyz_restraints,
)
print ' ', now(), 'xyz_md...'
models['xyz_md'] = cmod.build_md_model(
        num_particles,
        xyz_restraints=xyz_restraints,
        protocol=annealing_protocol,
        movie_path=movie_path_template.format('xyz_md'),
)
print ' ', now(), 'pair_md...'
models['pair_md'] = cmod.build_md_model(
        num_particles,
        pair_restraints=pair_restraints,
        protocol=annealing_protocol,
        movie_path=movie_path_template.format('pair_md'),
)
print ' ', now(), 'xyz_pair_min...'
models['xyz_pair_min'] = cmod.build_minimized_model(
        num_particles,
        xyz_restraints=xyz_restraints,
        pair_restraints=pair_restraints,
)
print ' ', now(), 'xyz_pair_md...'
models['xyz_pair_md'] = cmod.build_md_model(
        num_particles,
        xyz_restraints=xyz_restraints,
        pair_restraints=pair_restraints,
        protocol=annealing_protocol,
        movie_path=movie_path_template.format('xyz_pair_md'),
)
print ' ', now(), 'done.'
print

models['pair_md'] = cmod.toy.superimpose_model(
        models['pair_md'], reference)

for key, model in models.items():
    path = pdb_path_template.format(key)
    cmod.export_to_pdb(path, model, reference, xyz_restraints)

results = {
        'n_xyz': num_xyz_restraints,
        'n_pair': num_pair_restraints,
}
for key, model in models.items():
    results[key] = cmod.toy.evaluate_model(model, reference)

with open(json_path, 'w') as file:
    json.dump(results, file)

