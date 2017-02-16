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

from .qt_helper import new_spacer
from ..mpl.zoompan import ZoomPan
from ..mpl.imexam import Imexam
from ..mpl.ds9norm import DS9Normalize, DS9Interacter
from ..signal import MySignal

class FigConfigPanel(QtWidgets.QWidget):
    def __init__(self, hdu, ax, ds9normalize, parent=None):
        super(QtWidgets.QWidget, self).__init__(parent)
        self.hdu = hdu
        self.ax = ax
        self.ds9 = ds9normalize

        self.setFixedWidth(200)

class HDUFigureCanvas2D(QtWidgets.QWidget):
    def __init__(self, hdu, parent=None):
        self.hdu = hdu
        super(QtWidgets.QWidget, self).__init__(parent)

        self.layout = QtWidgets.QGridLayout(self)

        self.sidebar = QtWidgets.QWidget()
        self.sidebar_layout = QtWidgets.QGridLayout(self.sidebar)

        self.create_toolbar()
        self.create_plot()

        self.setMinimumSize(600, 400)

        self.layout.addWidget(self.toolbar, 0, 0, 1, 2)
        self.layout.addWidget(self.figCanvas, 1, 0, 1, 1)
        self.layout.addWidget(self.sidebar, 1, 1, 1, 1)

    def create_toolbar(self) :
        self.toolbar = QtWidgets.QToolBar(self)

        self.action_reset = QtWidgets.QAction(QtGui.QIcon.fromTheme('reset'),
                                             'Reset', self)

        self.action_config_plot = QtWidgets.QAction(QtGui.QIcon.fromTheme('document-properties'),
                                             'Configure', self)

        self.toolbar.addAction(self.action_reset)
        self.toolbar.addWidget(new_spacer())
        self.toolbar.addAction(self.action_config_plot)

    def create_plot(self):
        self.fig = Figure(figsize=(8, 8), dpi=90)
        self.ax = self.fig.add_axes((0.05, 0.05, 0.9, 0.9))
        self.ax.set_axis_off()
        self.figCanvas = FigureCanvas(self.fig)
        FigureCanvas.setSizePolicy(self.figCanvas,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)

        self._zp = ZoomPan()
        self._zpzoom = self._zp.zoom_factory(self.ax)
        self._zppan = self._zp.pan_factory(self.ax)

        self._ds9 = DS9Normalize(clip_lo = 30, clip_hi=99)
        self._ds9interact = DS9Interacter()
        self._ds9interact.ds9norm_factory(self.ax, self._ds9)

        self.plot_image()

        self.config_panel = FigConfigPanel(self.hdu, self.ax, self._ds9)

        #self.sidebar_layout.addWidget(self.config_panel)

    def plot_image(self):
        im = self.ax.imshow(self.hdu.data, norm=self._ds9, origin='lower')
