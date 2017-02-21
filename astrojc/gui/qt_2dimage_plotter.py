import os
import sys

import matplotlib
matplotlib.use('Qt5Agg')
from qtpy import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm, Normalize, PowerNorm
import numpy as np
from functools import partial

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

        self.setFixedWidth(250)

        self.start = MySignal(raise_error=False)
        self.stop = MySignal(raise_error=False)

class ImexamPanel(QtWidgets.QWidget):
    def __init__(self, hdu, ax, imexam, parent=None):
        super(QtWidgets.QWidget, self).__init__(parent)
        self.hdu = hdu
        self.ax = ax
        self.imexam = imexam

        self.plot_fig = Figure(dpi=60)
        self.plot_ax = self.plot_fig.add_subplot(111)

        self.setFixedWidth(250)

        self.clear_button = QtWidgets.QPushButton('Clear All')
        self.clear_button.clicked.connect(self._on_clear)

        self.layout = QtWidgets.QVBoxLayout(self)
        self.text = QtWidgets.QTextEdit()
        self.fig_canvas = FigureCanvas(self.plot_fig)
        self.fig_canvas.setMinimumSize(200,200)
        policy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred,
                                       QtWidgets.QSizePolicy.Preferred)
        policy.setHeightForWidth(True)
        self.fig_canvas.setSizePolicy(policy)

        self.layout.addWidget(self.fig_canvas)
        self.layout.addWidget(self.text)
        self.layout.addWidget(self.clear_button)

        self.imexam.print_text.connect(self.print_text)

        self.start = MySignal(raise_error=False)
        self.stop = MySignal(raise_error=False)

    def connect_imexam(self, ax):
        self.imexam.connect(ax, self.plot_ax)

    def disconnect_imexam(self):
        self.imexam.disconnect()
        self._on_clear()
        self.plot_fig.canvas.draw()

    def _on_clear(self):
        self.text.setText('')
        self.plot_ax.clear()

    def print_text(self, message):
        self.text.append(message)

class HDUFigureCanvas2D(QtWidgets.QWidget):
    def __init__(self, hdu, parent=None):
        self.hdu = hdu
        super(QtWidgets.QWidget, self).__init__(parent)

        self.layout = QtWidgets.QGridLayout(self)

        self.sidebar = QtWidgets.QTabWidget()
        self.sidebar.setTabsClosable(True)
        self.sidebar.setFixedWidth(250)
        self.sidebar.tabCloseRequested.connect(self.close_tab)
        self.sidebar_widget = {}

        self.create_toolbar()
        self.create_plot()

        self.setMinimumSize(700, 500)

        self.layout.addWidget(self.toolbar, 0, 0, 1, 2)
        self.layout.addWidget(self.figCanvas, 1, 0, 1, 1)

    def update_sidebar(self):
        if self.sidebar.count() > 0 and self.layout.indexOf(self.sidebar) == -1:
            self.layout.addWidget(self.sidebar, 1, 1, 1, 1)
        elif self.sidebar.count() == 0 and self.layout.indexOf(self.sidebar) != -1:
            self.layout.removeWidget(self.sidebar)
            self.sidebar.setParent(None)

    def create_toolbar(self) :
        self.toolbar = QtWidgets.QToolBar(self)
        self.toolbar_buttons = {}

        self.toolbar_buttons['config'] = QtWidgets.QAction(QtGui.QIcon.fromTheme('document-properties'),
                                             'Configure', self)
        self.toolbar_buttons['config'].setCheckable(True)
        self.toolbar_buttons['config'].triggered[bool].connect(self._on_config_changed)

        self.toolbar_buttons['imexam'] = QtWidgets.QAction(QtGui.QIcon.fromTheme('document-properties'),
                                             'imexam', self)
        self.toolbar_buttons['imexam'].setCheckable(True)
        self.toolbar_buttons['imexam'].triggered[bool].connect(self._on_imexam_changed)

        #self.toolbar.addWidget(new_spacer())
        self.toolbar.addAction(self.toolbar_buttons['config'])
        self.toolbar.addAction(self.toolbar_buttons['imexam'])

    def create_plot(self):
        self.fig = Figure(figsize=(8, 8), dpi=90)
        self.ax = self.fig.add_axes((0.0, 0.0, 1.0, 1.0))
        self.ax.set_axis_off()
        self.figCanvas = FigureCanvas(self.fig)
        FigureCanvas.setSizePolicy(self.figCanvas,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)

        self.figCanvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.figCanvas.setFocus()

        self._zp = ZoomPan()
        self._zpzoom = self._zp.zoom_factory(self.ax)
        self._zppan = self._zp.pan_factory(self.ax)

        self._ds9 = DS9Normalize(clip_lo = 30, clip_hi=99)
        self._ds9interact = DS9Interacter()
        self._ds9interact.ds9norm_factory(self.ax, self._ds9)

        self._imexam = Imexam()
        self._imexam.set_hdu(self.hdu)

        self.plot_image()

        self.sidebar_widget['config'] = FigConfigPanel(self.hdu, self.ax, self._ds9)
        self.sidebar_widget['imexam'] = ImexamPanel(self.hdu, self.ax, self._imexam)
        self.sidebar_widget['imexam'].start.connect(partial(self.sidebar_widget['imexam'].connect_imexam, self.ax))
        self.sidebar_widget['imexam'].stop.connect(self.sidebar_widget['imexam'].disconnect_imexam)

        self.sidebar_names = {'config': 'Configure',
                              'imexam': 'Examine'}

    def add_sidebar_mode(self, mode):
        if self.sidebar.indexOf(self.sidebar_widget[mode]) == -1:
            self.sidebar.addTab(self.sidebar_widget[mode], self.sidebar_names[mode])
            self.sidebar_widget[mode].start()
        if not self.toolbar_buttons[mode].isChecked():
            self.toolbar_buttons[mode].setChecked(True)
        self.update_sidebar()

    def remove_sidebar_mode(self, mode):
        if self.sidebar.indexOf(self.sidebar_widget[mode]) != -1:
            self.sidebar.removeTab(self.sidebar.indexOf(self.sidebar_widget[mode]))
            self.sidebar_widget[mode].stop()
        if self.toolbar_buttons[mode].isChecked():
            self.toolbar_buttons[mode].setChecked(False)
        self.update_sidebar()

    def close_tab(self, index):
        for i in self.sidebar_widget.keys():
            if self.sidebar_widget[i] == self.sidebar.widget(index):
                self.remove_sidebar_mode(i)
        self.fig.canvas.draw()

    def _on_config_changed(self, checked=False):
        if not checked:
            self.remove_sidebar_mode('config')
        else:
            self.add_sidebar_mode('config')

    def _on_imexam_changed(self, checked=False):
        if not checked:
            self.remove_sidebar_mode('imexam')
        else:
            self.add_sidebar_mode('imexam')

    def plot_image(self):
        im = self.ax.imshow(self.hdu.data, norm=self._ds9, origin='lower')
