import os
import sys

import matplotlib
matplotlib.use('Qt5Agg')
from qtpy import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm, Normalize, PowerNorm
import numpy as np

from .qt_helper import new_spacer
from ..mpl_helper import ZoomPan

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

        self.action_config_plot = QtWidgets.QAction(QtGui.QIcon.fromTheme('document-properties'),
                                             'Configure', self)
        self.config_pick = QtWidgets.QAction(QtGui.QIcon.fromTheme('document-properties'),
                                             'Configure Picker', self)
        self.action_pick = QtWidgets.QAction(QtGui.QIcon.fromTheme('find-location-symbolic'),
                                             'Pick Star', self)

        self.toolbar.addAction(self.action_pick)
        self.spacer = new_spacer()
        self.toolbar.addWidget(self.spacer)
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

        self.plot_config = {'norm' : Normalize(*self._get_min_max()),
                            'cmap' : 'viridis'}
        self.plot_image()

    def plot_image(self):
        im = self.ax.imshow(self.hdu.data, **self.plot_config, origin='lower')
        cbar = self.fig.colorbar(im)

    def _get_min_max(self):
        vmin = np.max([0, np.min(self.hdu.data)])
        vmax = np.max(np.sort(self.hdu.data, axis=None)[:int(self.properties['percentage']*len(self.hdu.data.flatten()) - 1)])
        return vmin, vmax
