#!/usr/bin/env python2

"""\
Usage: hilbert_demo.py [--steps=STEPS] [--test-run]
"""

import subprocess, numpy as np
from libraries import docopt, utils

args = docopt.docopt(__doc__)
num_particles = 512 if args['--test-run'] else 512
num_steps = int(args['--steps'] or 
    (1 if args['--test-run'] else (20 if utils.running_on_cluster() else 3)))
num_xyz_restraints = np.linspace(10, 50, num=num_steps).astype(int)
num_pair_restraints = np.linspace(150, 500, num=num_steps).astype(int)

params = []
for n_xyz in num_xyz_restraints:
    for n_pair in num_pair_restraints:
        params.append({
            'num_particles': num_particles,
            'num_xyz_restraints': n_xyz,
            'num_pair_restraints': n_pair,
        })

utils.clear_directories('jsons', 'pdbs', 'movies')
utils.submit_job('run_hilbert_demo.py', params)

