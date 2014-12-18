#!/usr/bin/env python2

import chromosome_modeler as cmod
import time, json

def run_algorithms(num_particles, num_xyz_restraints, num_pair_restraints):
    reference = cmod.toy.make_hilbert_reference(num_particles)
    xyz_restraints = cmod.toy.make_xyz_restraints(reference, num_xyz_restraints)
    pair_restraints = cmod.toy.make_pair_restraints(reference, num_pair_restraints)

    print 'num_xyz_restraints, num_pair_restraints = {}, {}'.format(
            num_xyz_restraints, num_pair_restraints)

    models = {}
    now = lambda: time.strftime("%X")

    print ' ', now(), 'iterp...'
    models['interp'] = cmod.build_interpolated_model(
            num_particles, xyz_restraints)

    print ' ', now(), 'xyz_min...'
    models['xyz_min'] = cmod.build_minimized_model(
            num_particles,
            xyz_restraints=xyz_restraints)

    print ' ', now(), 'xyz_md...'
    models['xyz_md'] = cmod.build_md_model(
            num_particles,
            xyz_restraints=xyz_restraints)

    print ' ', now(), 'pair_md...'
    models['pair_md'] = cmod.build_md_model(
            num_particles,
            pair_restraints=pair_restraints)

    print ' ', now(), 'xyz_pair_min...'
    models['xyz_pair_min'] = cmod.build_minimized_model(
            num_particles,
            xyz_restraints=xyz_restraints,
            pair_restraints=pair_restraints)

    print ' ', now(), 'xyz_pair_md...'
    models['xyz_pair_md'] = cmod.build_md_model(
            num_particles,
            xyz_restraints=xyz_restraints,
            pair_restraints=pair_restraints)

    print ' ', now(), 'done.'
    print

    json_path = 'hilbert_demo/jsons/{:02}_{:03}.json'.format(
            num_xyz_restraints, num_pair_restraints)
    pdb_path_template = 'hilbert_demo/pdbs/{:02}_{:03}_{{}}.pdb'.format(
            num_xyz_restraints, num_pair_restraints)

    for key, model in models.items():
        path = pdb_path_template.format(key)
        cmod.export_to_pdb(path, model, reference, xyz_restraints)

    results = {
            key: cmod.toy.evaluate_model(model, reference)
            for key, model in models.items()
    }
    results['n_xyz'] = num_xyz_restraints
    results['n_pair'] = num_pair_restraints

    with open(json_path, 'w') as file:
        json.dump(results, file)

if __name__ == '__main__':
    import numpy as np

    num_particles = 512 #4096
    num_xyz_restraints = np.linspace(10, 50, num=3, dtype=int)
    num_pair_restraints = np.linspace(150, 500, num=3, dtype=int)

    for n_xyz in num_xyz_restraints:
        for n_pair in num_pair_restraints:
            run_algorithms(num_particles, n_xyz, n_pair)

