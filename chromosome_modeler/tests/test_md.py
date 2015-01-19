#!/usr/bin/env python2

import chromosome_modeler as cmod

num_particles = 64
num_xyz_restraints = 10
reference = cmod.toy.make_hilbert_reference(num_particles)
protocol = cmod.define_logarithmic_protocol(3e4, 1e2, 1e-2, 100)
xyz_restraints = cmod.toy.make_xyz_restraints(reference, num_xyz_restraints)
pair_restraints = cmod.toy.make_pair_restraints(
        reference, max_distance=2.0, signal_to_noise=1.0)

system = cmod.define_system(
        num_particles,
        xyz_restraints=xyz_restraints,
        pair_restraints=pair_restraints,
        initial_coords='random',
        nonbonded='excluded-volume',
)
model = cmod.run_molecular_dynamics(
        system,
        protocol=protocol,
        frames=100,
        progress_bar=True,
        movie_path='md_test.pym',
)[-1]

cmod.print_restraint_satisfaction(system)
print 'Reference RMSD:', cmod.toy.evaluate_model(model, reference)
cmod.export_to_pdb('md_test.pdb', model, reference, xyz_restraints)



