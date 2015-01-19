#!/usr/bin/env python2

import finalexam
import numpy as np
import chromosome_modeler as cmod

np.set_printoptions(linewidth=1000)

@finalexam.test
def test_make_hilbert_reference():
    reference = cmod.toy.make_hilbert_reference(8)

    print repr(reference)

    assert tuple(reference[0]) == (0, 0, 0)
    assert tuple(reference[1]) == (1, 0, 0)
    assert tuple(reference[2]) == (1, 0, 1)
    assert tuple(reference[3]) == (0, 0, 1)
    assert tuple(reference[4]) == (0, 1, 1)
    assert tuple(reference[5]) == (1, 1, 1)
    assert tuple(reference[6]) == (1, 1, 0)
    assert tuple(reference[7]) == (0, 1, 0)

@finalexam.test
def test_make_nagano_refernce():
    nagano_ensemble = cmod.toy.make_nagano_reference()
    cmod.data.require_ensemble(nagano_ensemble)

@finalexam.test
def test_make_xyz_restraints():
    reference = cmod.toy.make_hilbert_reference(8)
    restraints = cmod.toy.make_xyz_restraints(reference, 2)
    assert all(np.abs(np.diff(restraints.keys())) >= 8 // 2)

@finalexam.test
def test_make_pair_restraints():
    reference = cmod.toy.make_hilbert_reference(8)
    restraints = cmod.toy.make_pair_restraints(
            reference, max_distance = 1.5, noise_weight=0)

    print repr(reference)
    print repr(restraints)

    expected_restraints = np.array([
        [ 0.        ,  1.        ,  1.41421356,  1.        ,  1.41421356, -1.5       ,  1.41421356,  1.        ],
        [ 1.        ,  0.        ,  1.        ,  1.41421356, -1.5       ,  1.41421356,  1.        ,  1.41421356],
        [ 1.41421356,  1.        ,  0.        ,  1.        ,  1.41421356,  1.        ,  1.41421356, -1.5       ],
        [ 1.        ,  1.41421356,  1.        ,  0.        ,  1.        ,  1.41421356, -1.5       ,  1.41421356],
        [ 1.41421356, -1.5       ,  1.41421356,  1.        ,  0.        ,  1.        ,  1.41421356,  1.        ],
        [-1.5       ,  1.41421356,  1.        ,  1.41421356,  1.        ,  0.        ,  1.        ,  1.41421356],
        [ 1.41421356,  1.        ,  1.41421356, -1.5       ,  1.41421356,  1.        ,  0.        ,  1.        ],
        [ 1.        ,  1.41421356, -1.5       ,  1.41421356,  1.        ,  1.41421356,  1.        ,  0.        ]])

    assert np.allclose(restraints, expected_restraints)

def test_evaluate_model():
    reference = np.array([
        [1.0, 0.0, 1.0],
        [2.0, 0.0, 1.0],
        [2.0, 0.0, 2.0],
        [1.0, 0.0, 2.0],
    ])

    model = np.array([
        [0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [3.0, 0.0, 3.0],
        [0.0, 0.0, 3.0],
    ])

    score = cmod.toy.evaluate_model(model, reference)
    assert score == np.sqrt(2), score


if __name__ == '__main__':
    finalexam.title("Testing the toy model...")
    finalexam.run()



