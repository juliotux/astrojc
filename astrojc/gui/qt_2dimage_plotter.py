import os
import sys

import matplotlib
matplotlib.use('Qt5Agg')
from qtpy import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm, Normalize, PowerNorm
import numpy as np

class ZoomPan:
    '''
    Activates zoom and pan with the scrollwhell in an axes.

    To use it, just enable it following:
    zp = ZoomPan()
    zoom = zp.zoom_factory(ax, base_scale=1.1)
    pan = zp.pan_factory(ax)
    '''
    def __init__(self):
        self.press = None
        self.cur_xlim = None
        self.cur_ylim = None
        self.x0 = None
        self.y0 = None
        self.x1 = None
        self.y1 = None
        self.xpress = None
        self.ypress = None

    def zoom_factory(self, ax, base_scale = 1.1):
        def zoom(event):
            cur_xlim = ax.get_xlim()
            cur_ylim = ax.get_ylim()

            xdata = event.xdata # get event x location
            ydata = event.ydata # get event y location

            if event.button == 'down':
                # deal with zoom in
                scale_factor = 1 / base_scale
            elif event.button == 'up':
                # deal with zoom out
                scale_factor = base_scale
            else:
                # deal with something that should never happen
                scale_factor = 1
                print(event.button)

            new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
            new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor

            relx = (cur_xlim[1] - xdata)/(cur_xlim[1] - cur_xlim[0])
            rely = (cur_ylim[1] - ydata)/(cur_ylim[1] - cur_ylim[0])

            ax.set_xlim([xdata - new_width * (1-relx), xdata + new_width * (relx)])
            ax.set_ylim([ydata - new_height * (1-rely), ydata + new_height * (rely)])
            ax.figure.canvas.draw()

        fig = ax.get_figure() # get the figure of interest
        fig.canvas.mpl_connect('scroll_event', zoom)

        return zoom

    def pan_factory(self, ax):
        #TODO: Only mid button for pan. Right button for brigthness/contrast
        def onPress(event):
            if event.inaxes != ax: return
            self.cur_xlim = ax.get_xlim()
            self.cur_ylim = ax.get_ylim()
            self.press = self.x0, self.y0, event.xdata, event.ydata
            self.x0, self.y0, self.xpress, self.ypress = self.press

        def onRelease(event):
            self.press = None
            ax.figure.canvas.draw()

        def onMotion(event):
            if self.press is None: return
            if event.inaxes != ax: return
            dx = event.xdata - self.xpress
            dy = event.ydata - self.ypress
            self.cur_xlim -= dx
            self.cur_ylim -= dy
            ax.set_xlim(self.cur_xlim)
            ax.set_ylim(self.cur_ylim)

            ax.figure.canvas.draw()

        fig = ax.get_figure() # get the figure of interest

        # attach the call back
        fig.canvas.mpl_connect('button_press_event',onPress)
        fig.canvas.mpl_connect('button_release_event',onRelease)
        fig.canvas.mpl_connect('motion_notify_event',onMotion)

        #return the function
        return onMotion

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
        self.spacer = QtWidgets.QWidget()
        self.spacer.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.toolbar.addWidget(self.spacer)
        self.toolbar.addAction(self.action_config_plot)

    def create_plot(self):
        self.fig = Figure(figsize=(8, 8), dpi=100)
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
        im = self.ax.imshow(self.hdu.data, **self.plot_config)
        cbar = self.fig.colorbar(im)

    def _get_min_max(self):
        vmin = np.max([0, np.min(self.hdu.data)])
        vmax = np.max(np.sort(self.hdu.data, axis=None)[:int(self.properties['percentage']*len(self.hdu.data.flatten()) - 1)])
        return vmin, vmax
