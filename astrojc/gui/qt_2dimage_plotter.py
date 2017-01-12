import os
import sys

import matplotlib
matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np

from ..mpl_helper import ZoomPan

class HDUFigureCanvas2D(QtWidgets.QWidget):
    def __init__(self, hdu, parent=None):
        self.hdu = hdu
        super(QtWidgets.QWidget, self).__init__(parent)

        self.fig = Figure(figsize=(4, 4), dpi=100)
        self.ax = self.fig.add_axes((0.05, 0.05, 0.9, 0.9))
        self.figCanvas = FigureCanvas(self.fig)
        FigureCanvas.setSizePolicy(self.figCanvas,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)

        self._zp = ZoomPan()
        self._zpzoom = self._zp.zoom_factory(self.ax)
        self._zppan = self._zp.pan_factory(self.ax)

        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.figCanvas)
        self.plot_config = {'vmin' : np.nanmin(self.hdu.data),
                            'vmax' : np.nanmax(self.hdu.data),
                            'cmap' : 'viridis'}

        self.plot_image()

    def plot_image(self):
        self.ax.imshow(self.hdu.data, **self.plot_config)
