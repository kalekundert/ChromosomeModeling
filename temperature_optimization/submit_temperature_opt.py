#!/usr/bin/env python

"""\
Usage: submit_temperature_opt.py [--test-run]
"""

import numpy as np
import chromosome_modeler as cmod
from libraries import docopt, utils

args = docopt.docopt(__doc__)

params = []
step_params = np.linspace(3, 100, 10).astype(int)
start_temp_params = 10**np.linspace(1, 4, num=10)
stop_temp_params = 10**np.linspace(-2, 2, num=10)

for steps in step_params:
    for start_temp in start_temp_params:
        for stop_temp in stop_temp_params:
            if start_temp >= stop_temp:
                params.append({
                    'num_particles': 4096,
                    'annealing_protocol': cmod.define_linear_protocol(
                        3e4, start_temp, stop_temp, steps),
                })

if args['--test-run']:
    params = [params[0]]
    params[0]['num_particles'] = 512

utils.clear_directories('jsons', 'pdbs', 'movies')
utils.submit_job('run_temperature_opt.py', params, args['--test-run'])

