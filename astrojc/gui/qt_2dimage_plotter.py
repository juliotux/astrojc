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
from ..mpl_helper import ZoomPan
from .imexam import Imexam
from .my_signal import MySignal

default_im_rect = [0.0, 0.0, 1.0, 1.0]
picked_im_rect = [0.0, 0.1, 0.8, 0.8]
pick_rect = [[0.82, 0.1, 0.18, 0.25],
             [0.82, 0.4, 0.18, 0.25],
             [0.82, 0.7, 0.18, 0.25]]

class Picker():
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
        self.connected = None

        self.output_panel = None

        self.imexam = Imexam()
        self.imexam.imexam_factory(self.ax, None)
        self.imexam.set_hdu(self.hdu)

        self.create_auxiliar_axes = MySignal() #args: output_panel
        self.create_output_panel = MySignal()  #No Arguments
        self.print_on_panel = MySignal()       #args: output_panel, message
        self.clear_output_panel = MySignal()   #args: output_panel
        self.imexam.create_plot_axes.connect(self._create_imexam_axes)
        self.imexam.print_text.connect(self.print_text)

    def _create_imexam_axes(self):
        if self.output_panel is None:
            self.output_panel = self.create_output_panel.emit()
        return self.create_auxiliar_axes(self.output_panel)

    def print_text(self, message):
        if self.output_panel is None:
            self.output_panel = self.create_output_panel.emit()
        self.print_on_panel.emit(self.output_panel, message)

    def onPress(self, event):
        if event.inaxes == self.ax:
            if event.button == 1:
                self.pick = (event.xdata, event.ydata)
                self.do_pick(self.pick)

    def do_pick(self, pick):
        if self.output_panel is None:
            self.output_panel = self.create_output_panel.emit()
        self.clear_output_panel.emit(self.output_panel)

        ax = self.create_auxiliar_axes.emit(self.output_panel)
        self.imexam.plot_ax = ax
        self.imexam.do_option('r', pick[0], pick[1])
        ax = self.create_auxiliar_axes.emit(self.output_panel)
        self.imexam.plot_ax = ax
        self.imexam.do_option('e', pick[0], pick[1])
        ax = self.create_auxiliar_axes.emit(self.output_panel)
        self.imexam.plot_ax = ax
        self.imexam.do_option('s', pick[0], pick[1])
        self.imexam.do_option('a', pick[0], pick[1])
        self.imexam.plot_ax = None

    def connect(self):
        if self.connected == None:
            self.connected = self.fig.canvas.mpl_connect('button_press_event',
                                                         self.onPress)

    def disconnect(self):
        if self.connected is not None:
            self.fig.canvas.mpl_disconnect(self.connected)
            self.connected = None
            self.fig.canvas.draw()

class ImOutputPanel(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super(QtWidgets.QWidget, self).__init__(parent)

        self.create_layout()

    def create_layout(self):
        self.layout = QtWidgets.QVBoxLayout(self)
        self.setMinimumWidth(150)

        self.text = QtWidgets.QTextEdit(self)
        self.text.setReadOnly(True)
        self.text.setMinimumSize(150,100)
        self.scroll = QtWidgets.QScrollArea(self)
        self.layout.addWidget(self.text)
        self.layout.addWidget(self.scroll)
        self.fig_layout = QtWidgets.QVBoxLayout(self.scroll)

    def add_axes(self):
        fig = Figure(dpi=90)
        ax = fig.gca()
        figcanvas = FigureCanvas(fig)
        FigureCanvas.setSizePolicy(figcanvas,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
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

        self.setMinimumSize(400, 300)

        self.layout.addWidget(self.toolbar)
        self.layout2.addWidget(self.figCanvas)
        self.layout.addLayout(self.layout2)

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
        self._pick.create_output_panel.connect(self._create_output_panel)
        self._pick.create_auxiliar_axes.connect(self._create_plot_ax)
        self._pick.print_on_panel.connect(self._write_on_panel)
        self._pick.clear_output_panel.connect(self._clear_output_panel)

        self.plot_config = {'norm' : Normalize(*self._get_min_max()),
                            'cmap' : 'viridis'}
        self.plot_image()

    def _create_output_panel(self):
        self.output_panel = ImOutputPanel()
        self.layout2.addWidget(self.output_panel)
        return self.output_panel

    def _clear_output_panel(self, panel):
        panel.clear()

    def _create_plot_ax(self, panel):
        return panel.add_axes()

    def _write_on_panel(self, panel, message):
        panel.print_text(message)

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
