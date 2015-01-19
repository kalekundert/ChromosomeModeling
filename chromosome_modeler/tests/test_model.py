#!/usr/bin/env python2

import finalexam
import numpy as np
import chromosome_modeler as cmod

@finalexam.test
def test_define_linear_protocol():
    protocol = cmod.define_linear_protocol(3e4, 1e3, 1e-2, 5)
    expected_protocol = [
            (6000.0, 1000.0),
            (6000.0, 750.0025), 
            (6000.0, 500.005),
            (6000.0, 250.0075),
            (6000.0, 0.01),
    ]
    assert np.allclose(protocol, expected_protocol), protocol

@finalexam.test
def test_define_logarithmic_protocol():
    protocol = cmod.define_logarithmic_protocol(3e4, 1e3, 1e-2, 5)
    expected_protocol = [
            (6000.0, 1000.0),
            (6000.0, 56.234132519034908),
            (6000.0, 3.1622776601683795),
            (6000.0, 0.17782794100389229),
            (6000.0, 0.01),
    ]
    assert np.allclose(protocol, expected_protocol), protocol


if __name__ == '__main__':
    finalexam.title("Testing the chromosome model...")
    finalexam.run()



