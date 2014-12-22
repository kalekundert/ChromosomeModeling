#!/usr/bin/env python2

import sys, time, json
import chromosome_modeler as cmod
from libraries import utils

params = utils.read_params()
num_particles = params['num_particles']
annealing_protocol = params['annealing_protocol']
num_xyz_restraints = 30
num_pair_restraints = 500

name = '_'.join('{0}@{1}'.format(i,T) for i,T in annealing_protocol)
json_path = 'jsons/{0}.json'.format(name)
pdb_path = 'pdbs/{0}.pdb'.format(name)
movie_path = 'movies/{0}.pym'.format(name)

reference = cmod.toy.make_hilbert_reference(num_particles)
xyz_restraints = cmod.toy.make_xyz_restraints(reference, num_xyz_restraints)
pair_restraints = cmod.toy.make_pair_restraints(reference, num_pair_restraints)

interp_model = cmod.build_interpolated_model(
        num_particles,
        xyz_restraints,
)
md_model = cmod.build_md_model(
        num_particles,
        xyz_restraints=xyz_restraints,
        pair_restraints=pair_restraints,
        protocol=annealing_protocol,
        movie_path=movie_path,
)

cmod.export_to_pdb(pdb_path, md_model, reference, xyz_restraints)

results = {
        'n_xyz': num_xyz_restraints,
        'n_pair': num_pair_restraints,
        'interp': cmod.toy.evaluate_model(interp_model, reference),
        'md': cmod.toy.evaluate_model(md_model, reference),
}
with open(json_path, 'w') as file:
    json.dump(results, file)


