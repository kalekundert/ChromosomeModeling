#!/usr/bin/env python

import glob, json, operator
from pprint import pprint

data = []

for path in glob.glob('jsons/*.json'):
    with open(path) as file:
        datum = json.load(file)
        datum['path'] = path
        data.append(datum)

data = sorted(data, key=operator.itemgetter('md'))[::-1]
pprint(data)
