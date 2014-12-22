#!/usr/bin/env python2

import sys, json, glob
import numpy as np
import pandas as pd
import pylab as pp


import gtk
import matplotlib
import matplotlib.pyplot
import pango

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg

cols, rows = pd.util.terminal.get_terminal_size()
pd.set_option('display.max_rows', rows)
pd.set_option('display.max_columns', cols)
pd.set_option('display.width', cols)

class DemoViewer (gtk.Window):

    def __init__(self, data):
        gtk.Window.__init__(self)

        self.data = data
        self.current_axis = 'n_pair'
        self.current_slice = 0

        self.setup_plot()
        self.update_plot()

    def setup_plot(self):
        self.figure = Figure(facecolor='#edecea')
        self.axes = self.figure.add_axes((0.15, 0.15, 0.75, 0.75))

        box = self.axes.get_position()
        self.axes.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        canvas = FigureCanvasGTKAgg(self.figure)
        toolbar = NavigationToolbar2GTKAgg(canvas, self)

        vbox = gtk.VBox()
        vbox.pack_start(canvas)
        vbox.pack_start(toolbar, expand=False)

        self.add(vbox)
        self.connect('destroy', gtk.main_quit)
        self.connect('scroll-event', self.on_scroll)
        self.set_default_size(644, 365)
        self.show_all()

    def update_plot(self):
        self.axes.clear()

        off_axis = 'n_pair' if self.current_axis == 'n_xyz' else 'n_xyz'
        off_axis_data = self.data[off_axis]
        off_axis_slices = sorted(off_axis_data.unique())

        self.current_slice = np.clip(
                self.current_slice, 0, len(off_axis_slices) - 1)

        off_axis_slice = off_axis_slices[self.current_slice]

        slice = self.data[off_axis_data == off_axis_slice]
        slice = slice.drop([off_axis], axis=1)
        slice = slice.sort([self.current_axis])
        slice.plot(x=self.current_axis, ax=self.axes)

        x_axis = '{} ({}={})'.format(
                self.current_axis, off_axis, off_axis_slice)

        self.axes.set_xlabel(x_axis)
        self.axes.set_ylabel('RMSD')
        self.axes.set_ylim(0, 12)
        self.axes.legend(loc=2, bbox_to_anchor=(1.02, 1.0))

        self.figure.canvas.draw()

    def on_scroll(self, widget, event):
        if event.direction == gtk.gdk.SCROLL_UP:
            self.current_slice += 1
        elif event.direction == gtk.gdk.SCROLL_DOWN:
            self.current_slice -= 1
        self.update_plot()



json_paths = glob.glob('jsons/*.json')

data = []

for i, path in enumerate(json_paths):
    sys.stdout.write('\r[{}/{}]'.format(i+1, len(json_paths)))
    sys.stdout.flush()

    with open(path) as file:
        data.append(json.load(file))

print

data = pd.DataFrame(data)

DemoViewer(data)
gtk.main()
