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

class ImOutputPanel(QtWidgets.QWidget):
    def __init__(self, hdu, parent=None):
        super(QtWidgets.QWidget, self).__init__(parent)

        self.hdu = hdu
        self.create_layout()
        self.setFixedWidth(200)

    def create_layout(self):
        self.layout = QtWidgets.QVBoxLayout(self)

        #self.mapfig = Figure(figsize=(2,2), dpi=90)
        #self.mapax = self.mapfig.gca()
        #self.immap = FigureCanvas(self.mapfig)
        #self.immap.setFixedSize(180,180)

        self.text = QtWidgets.QTextEdit(self)
        self.text.setReadOnly(True)
        self.text.setFixedSize(180,100)

        self.scroll = QtWidgets.QScrollArea(self)
        self.scroll.setFixedSize(180,180)

        #self.layout.addWidget(self.immap)
        self.layout.addWidget(self.text)
        self.layout.addWidget(self.scroll)
        self.fig_layout = QtWidgets.QVBoxLayout(self.scroll)

    def add_axes(self):
        fig = Figure(figsize=(2,2), dpi=90)
        ax = fig.gca()
        figcanvas = FigureCanvas(fig)
        FigureCanvas.setSizePolicy(figcanvas,
                                   QtWidgets.QSizePolicy.Fixed,
                                   QtWidgets.QSizePolicy.Fixed)
        self.fig_layout.addWidget(figcanvas)
        return ax

    def remove_axes(self, ax):
        self.layout.removeWidget(ax.get_figure().canvas)
        ax.get_figure().canvas.deleteLater()

    def print_text(self, message):
        self.text.append(message)

    def clear(self):
        for i in reversed(range(self.fig_layout.count())):
            self.fig_layout.itemAt(i).widget().deleteLater()
        self.text.clear()

class HDUFigureCanvas2D(QtWidgets.QWidget):
    def __init__(self, hdu, parent=None):
        self.hdu = hdu
        super(QtWidgets.QWidget, self).__init__(parent)

        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout2 = QtWidgets.QHBoxLayout()

        self.properties = {'percentage' : 0.95}

        self.create_toolbar()
        self.create_plot()
        self.create_output_panel()

        self.setMinimumSize(400, 300)

        self.layout.addWidget(self.toolbar)
        self.layout2.addWidget(self.figCanvas)
        self.layout2.addWidget(self.output_panel)
        self.layout.addLayout(self.layout2)

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
        self.ax = self.fig.add_axes((0.0, 0.0, 1.0, 1.0))
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

    def create_output_panel(self):
        self.output_panel = ImOutputPanel(self.hdu)

        self._imexam = Imexam()
        self._imexam.imexam_connect(self.ax)
        self._imexam.create_plot_axes.connect(self.output_panel.add_axes)
        self._imexam.print_text.connect(self.output_panel.print_text)
        self._imexam.set_hdu(self.hdu)

    def _clear_output_panel(self, panel):
        panel.clear()

    def _create_plot_ax(self, panel):
        return panel.add_axes()

    def _write_on_panel(self, panel, message):
        panel.print_text(message)

    def plot_image(self):
        im = self.ax.imshow(self.hdu.data, norm=self._ds9, origin='lower')
        #cbar = self.fig.colorbar(im)

    def _get_min_max(self):
        vmin = np.max([0, np.min(self.hdu.data)])
        vmax = np.max(np.sort(self.hdu.data, axis=None)[:int(self.properties['percentage']*len(self.hdu.data.flatten()) - 1)])
        return vmin, vmax
