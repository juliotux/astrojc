import os
import sys

import matplotlib
matplotlib.use('Qt5Agg')
from qtpy import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm, Normalize, PowerNorm
import numpy as np

from ..reducing.psf_fitting import (extract_data, xy2r,
                                    fit_moffat_spatial, fit_moffat_radial,
                                    fit_gaussian_spatial, fit_gaussian_radial)
from photutils

from .qt_helper import new_spacer
from ..mpl_helper import ZoomPan

default_im_rect = [0.0, 0.0, 1.0, 1.0]
picked_im_rect = [0.0, 0.1, 0.8, 0.8]
pick_rect = [[0.82, 0.1, 0.18, 0.25],
             [0.82, 0.4, 0.18, 0.25],
             [0.82, 0.7, 0.18, 0.25]]

class Picker():
    #TODO: add a method to create external plots out of the figure (dock widgets?)
    def __init__(self, ax, hdu, boxsize=30, model='gaussian', n_contour=10,
                 surface_wire = False, mode='radial'):
        self.hdu = hdu
        self.ax = ax
        self.fig = ax.get_figure()

        self.boxsize = boxsize
        self.model = model
        self.n_contour = n_contour
        self.surface_wire = surface_wire

        self.pick = None
        self.naxes = None
        self.connected = None

    def onPress(self, event):
        if event.inaxes == self.ax:
            if event.button == 1:
                self.pick = (event.xdata, event.ydata)
                self.do_pick(self.pick)

    def do_pick(self, pick):
        if self.naxes == None or len(self.naxes) == 0:
            self.naxes = [self.fig.add_axes(i) for i in pick_rect]
            self.ax.set_position(picked_im_rect)

    def connect(self):
        if self.connected == None:
            self.connected = self.fig.canvas.mpl_connect('button_press_event',
                                                         self.onPress)

    def remove_auxiliary_axes(self):
        if self.naxes is not None:
            for i in self.naxes:
                i.remove()
            self.naxes = None
            self.ax.set_position(default_im_rect)

    def disconnect(self):
        if self.connected is not None:
            self.fig.canvas.mpl_disconnect(self.connected)
            self.connected = None
            self.remove_auxiliary_axes()
            self.fig.canvas.draw()

class HDUFigureCanvas2D(QtWidgets.QWidget):
    def __init__(self, hdu, parent=None):
        self.hdu = hdu
        super(QtWidgets.QWidget, self).__init__(parent)

        self.properties = {'percentage' : 0.95}

        self.create_toolbar()
        self.create_plot()

        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.toolbar)
        self.layout.addWidget(self.figCanvas)

    def create_toolbar(self) :
        self.toolbar = QtWidgets.QToolBar(self)

        self.action_reset = QtWidgets.QAction(QtGui.QIcon.fromTheme('reset'),
                                             'Reset', self)

        self.action_config_plot = QtWidgets.QAction(QtGui.QIcon.fromTheme('document-properties'),
                                             'Configure', self)

        self.config_pick = QtWidgets.QAction(QtGui.QIcon.fromTheme('document-properties'),
                                             'Configure Picker', self)

        self.action_pick = QtWidgets.QAction(QtGui.QIcon.fromTheme('find-location-symbolic'),
                                             'Pick Star', self)
        self.action_pick.setCheckable(True)
        self.action_pick.changed.connect(self.toogle_pick)

        self.toolbar.addAction(self.action_reset)
        self.toolbar.addWidget(new_spacer())
        self.toolbar.addAction(self.action_pick)
        self.toolbar.addWidget(new_spacer())
        self.toolbar.addAction(self.action_config_plot)

    def create_plot(self):
        self.fig = Figure(figsize=(8, 8), dpi=90)
        self.ax = self.fig.add_axes((0.0, 0.0, 1.0, 1.0))
        self.figCanvas = FigureCanvas(self.fig)
        FigureCanvas.setSizePolicy(self.figCanvas,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)

        self._zp = ZoomPan()
        self._zpzoom = self._zp.zoom_factory(self.ax)
        self._zppan = self._zp.pan_factory(self.ax)

        self._pick = Picker(self.ax, self.hdu)

        self.plot_config = {'norm' : Normalize(*self._get_min_max()),
                            'cmap' : 'viridis'}
        self.plot_image()

    def plot_image(self):
        im = self.ax.imshow(self.hdu.data, **self.plot_config, origin='lower')
        #cbar = self.fig.colorbar(im)

    def toogle_pick(self):
        if self._pick.connected is None:
            self._pick.connect()
        else:
            self._pick.disconnect()

    def _get_min_max(self):
        vmin = np.max([0, np.min(self.hdu.data)])
        vmax = np.max(np.sort(self.hdu.data, axis=None)[:int(self.properties['percentage']*len(self.hdu.data.flatten()) - 1)])
        return vmin, vmax

    def reset_fig(self):
        for i in self._pick.naxes:
            i.remove()
