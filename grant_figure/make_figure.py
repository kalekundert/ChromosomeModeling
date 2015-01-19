#!/usr/bin/env python

from pylab import *
import collections

Methods = collections.namedtuple('Methods', 'pair xyz both')

rmsds = Methods(4.00562362662, 2.63649056823, 0.952165881042)

width = 0.6
offsets = arange(len(rmsds)) + 0.5 * (1 - width)

bar(offsets, rmsds, width, color='k')
xticks(offsets + 0.5 * width,
        ('Hi-C', 'Microscopy', 'Hi-C + Microscopy'), ha='center')

ylabel('RMSD')

show()
