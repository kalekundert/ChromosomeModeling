#!/usr/bin/env python2

"""\
Usage: hilbert_demo.py [--steps=STEPS] [--test-run]
"""

import subprocess, numpy as np
from libraries import docopt, utils

args = docopt.docopt(__doc__)

num_particles = 512 if args['--test-run'] else 4096
num_steps = int(args['--steps'] or 
    (1 if args['--test-run'] else (20 if utils.running_on_chef() else 3)))
num_xyz_restraints = np.linspace(10, 50, num=num_steps, dtype=int)
num_pair_restraints = np.linspace(150, 500, num=num_steps, dtype=int)

def submit_job(num_xyz_restraints, num_pair_restraints):
    if utils.running_on_chef():
        submit_cluster_job(num_xyz_restraints, num_pair_restraints)
    else:
        submit_local_job(num_xyz_restraints, num_pair_restraints)

def submit_cluster_job(num_xyz_restraints, num_pair_restraints):
    qsub_command = (
            'qsub',
            'run_hilbert_demo.py',
            num_particles,
            num_xyz_restraints,
            num_pair_restraints,
    )
    subprocess.call([str(x) for x in qsub_command])

def submit_local_job(num_xyz_restraints, num_pair_restraints):
    local_command = (
            './run_hilbert_demo.py',
            num_particles,
            num_xyz_restraints,
            num_pair_restraints,
    )
    subprocess.call([str(x) for x in local_command])


utils.clear_directory('jsons')
utils.clear_directory('pdbs')
utils.clear_directory('movies')

for n_xyz in num_xyz_restraints:
    for n_pair in num_pair_restraints:
        submit_job(n_xyz, n_pair)
