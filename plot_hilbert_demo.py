#!/usr/bin/env python2

import json, glob
import numpy as np
import pandas as pd
import pylab as pp

json_paths = glob.glob('hilbert_demo/jsons/*.json')

data = []

for path in json_paths:
    with open(path) as file:
        data.append(json.load(file))

data = pd.DataFrame(data)

pair_slice = data[data['n_pair'] == np.median(data['n_pair'])]
pair_slice.drop(['n_pair'], axis=1)

xyz_slice = data[data['n_xyz'] == np.median(data['n_xyz'])]
xyz_slice.drop(['n_xyz'], axis=1)

#pp.subplot(121)
pair_slice.plot(x='n_xyz')
#pp.subplot(122)
#xyz_slice.plot(x='n_pair')
pp.show()

