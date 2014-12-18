#!/usr/bin/env python2

from __future__ import division
import collections, numpy as np, scipy as sp

class xyz_restraints (collections.OrderedDict):

    def __init__(self):
        super(xyz_restraints, self).__init__()

    def __str__(self):
        from pprint import pformat
        return pformat(dict(self))

    def __setitem__(self, key, value):
        require_coord(value)
        return super(xyz_restraints, self).__setitem__(key, value)
    

class pair_restraints (collections.Counter):

    def __init__(self):
        super(pair_restraints, self).__init__()

    def __str__(self):
        from pprint import pformat
        return pformat(dict(self))

    def __getitem__(self, key):
        return super(pair_restraints, self).__getitem__(self._key(key))

    def __setitem__(self, key, value):
        return super(pair_restraints, self).__setitem__(self._key(key), value)
    
    def _key(self, key):
        assert len(key) == 2
        return tuple(sorted(key))



def are_coords(coords):
    try: require_coords(coords)
    except: return False
    else: return True


def require_coord(coord):
    assert coord.shape == (3,)

def require_coords(coords):
    assert coords.ndim == 2
    assert coords.shape[1] == 3, coords.shape

def require_ensemble(ensemble):
    assert ensemble.ndim == 3
    assert ensemble.shape[2] == 3

def require_trajectory(trajectory):
    if trajectory: require_coords(trajectory[0])

def require_ensemble_trajectory(ens_traj):
    if ens_traj: require_ensemble(ens_traj[0])
