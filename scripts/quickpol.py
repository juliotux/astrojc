import os
import sys
sys.path.append('build/lib/')

import matplotlib
matplotlib.use('Qt5Agg')
from qtpy import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import colors
import numpy as np
import glob

from astropy.io import fits

from astrojc.gui.qt_helper import new_spacer
from astrojc.mpl.ds9norm import DS9Normalize, DS9Interacter
from astrojc.mpl.zoompan import ZoomPan
from astrojc.logging import log


class CamPolMW(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(QtWidgets.QMainWindow, self).__init__(parent)
        self.setWindowTitle('Quick Polarimetry')
        self.setMinimumSize(800, 600)

        self.files = []

        self.main_widget = QtWidgets.QWidget()
        self.main_grid = QtWidgets.QGridLayout(self.main_widget)

        #Main Menu
        menu = QtWidgets.QToolBar()
        menu.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                           QtWidgets.QSizePolicy.Fixed)
        menu_open = menu.addAction('Open')
        menu_open.triggered.connect(self._open_files)
        menu_clear = menu.addAction('Clear')
        menu_clear.triggered.connect(self._clear_list)
        menu.addSeparator()
        menu_auto = menu.addAction('Auto Pol.')
        menu_seto = menu.addAction('Set Star')
        menu.addWidget(new_spacer())
        menu_seto = menu.addAction('Config.')

        #Display image figure
        self.im_fig = Figure()
        self.im_ax = self.im_fig.add_subplot(111)
        image_widget = FigureCanvas(self.im_fig)
        self._ds9 = DS9Normalize(clip_lo = 1, clip_hi=99)
        self._ds9interact = DS9Interacter()
        self._ds9interact.ds9norm_factory(self.im_ax, self._ds9)
        self._zp = ZoomPan()
        self._zpzoom = self._zp.zoom_factory(self.im_ax)
        self._zppan = self._zp.pan_factory(self.im_ax)
        self.im_fig.canvas.draw()

        #File List View
        list_view = QtWidgets.QListWidget()
        #list_view.currentItemChanged.connect(self._on_item_clicked)
        list_view.setSizePolicy(QtWidgets.QSizePolicy.Minimum,
                                QtWidgets.QSizePolicy.Expanding)
        list_widget = QtWidgets.QWidget()
        list_buttons = QtWidgets.QHBoxLayout(list_widget)
        add_button = QtWidgets.QPushButton()
        add_button.setText('+')
        add_button.clicked.connect(self._open_files)
        remove_button = QtWidgets.QPushButton()
        remove_button.setText('-')
        remove_button.clicked.connect(self._remove_files)
        clear_button = QtWidgets.QPushButton()
        clear_button.setText('clear')
        clear_button.clicked.connect(self._clear_list)
        list_buttons.addWidget(add_button)
        list_buttons.addWidget(remove_button)
        list_buttons.addWidget(clear_button)

        #Plot Results
        self.res_fig = Figure(figsize=(4, 3))
        self.res_ax = self.res_fig.add_subplot(111)
        self.res_ax.set_xlim(0, 360)
        self.res_ax.set_xlabel('Retarder Position')
        self.res_ax.set_ylabel('Flux Ratio')
        self.res_fig.tight_layout()
        res_widget = QtWidgets.QWidget()
        res_widget.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                 QtWidgets.QSizePolicy.Minimum)
        res_layout = QtWidgets.QVBoxLayout(res_widget)
        res_fig_widget = FigureCanvas(self.res_fig)
        res_fig_widget.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                     QtWidgets.QSizePolicy.Minimum)
        res_table_widget = QtWidgets.QTableWidget()
        res_layout.addWidget(res_fig_widget)
        res_layout.addWidget(res_table_widget)

        # Grid positioning
        self.main_grid.addWidget(menu, 0, 0, 1, 3)
        self.main_grid.addWidget(list_view, 1, 0, 1, 1)
        self.main_grid.addWidget(list_widget, 2, 0, 1, 1)
        self.main_grid.addWidget(image_widget, 1, 1, 2, 1)
        self.main_grid.addWidget(res_widget, 1, 2, 2, 1)
        self.setCentralWidget(self.main_widget)

    def _clear_list(self):
        '''Clear the images list'''

    def _open_files(self):
        '''Open the list of files'''

    def _remove_files(self):
        '''Remove the selected files'''



def main(argv):
    app = QtWidgets.QApplication(argv)
    mw = CamPolMW()
    mw.show()

    sys.exit(app.exec_())

if __name__ == '__main__':
    main(sys.argv)
