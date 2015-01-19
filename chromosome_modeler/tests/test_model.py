#!/usr/bin/env python2

import chromosome_modeler as cmod

num_particles = 64
num_xyz_restraints = 10

reference = cmod.toy.make_hilbert_reference(num_particles)
xyz_restraints = cmod.toy.make_xyz_restraints(reference, num_xyz_restraints)
pair_restraints = cmod.toy.make_pair_restraints(reference, max_distance=1.5, noise_weight=0)

system = cmod.define_system(
        num_particles,
        xyz_restraints=xyz_restraints,
        pair_restraints=pair_restraints,
        initial_coords='interpolated',
        nonbonded='excluded-volume',
)
model = cmod.run_molecular_dynamics(
        system,
        frames=100,
        progress_bar=True,
        movie_path='md_test.pym',
)[-1]

cmod.print_restraint_satisfaction(system)
print 'Reference RMSD:', cmod.toy.evaluate_model(model, reference)
cmod.export_to_pdb('md_test.pdb', model, reference, xyz_restraints)



