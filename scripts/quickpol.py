import os
import sys
sys.path.append('build/lib/')

import matplotlib
matplotlib.use('Qt5Agg')
from qtpy import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import colors
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy as np
import glob

from astropy.io import fits

from astrojc.gui.qt_helper import new_spacer
from astrojc.mpl.ds9norm import DS9Normalize, DS9Interacter
from astrojc.mpl.zoompan import ZoomPan
from astrojc.logging import log




def clear_ax_by_label(ax, label):
    '''Clear all the artists with a given label from an axes.'''
    for i in ax.get_lines():
        if i.get_label() == label:
            i.remove()
    for i in ax.patches:
        if i.get_label() == label:
            i.remove()


class CamPolMW(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(QtWidgets.QMainWindow, self).__init__(parent)
        self.setWindowTitle('Quick Polarimetry')
        self.setMinimumSize(1000, 600)

        self._x1, self._y1 = (None, None)
        self._x2, self._y2 = (None, None)
        self.opened = None

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
        menu_auto.triggered.connect(self._on_autoprocess_clicked)
        self.menu_seto = menu.addAction('Set Star')
        self.menu_seto.toggled.connect(self._on_set_star_clicked)
        self.menu_seto.setCheckable(True)
        menu.addWidget(new_spacer())
        menu_seto = menu.addAction('Config.')

        #Display image figure
        self.im_fig = Figure()
        self.im_ax = self.im_fig.add_subplot(111)
        image_widget = FigureCanvas(self.im_fig)
        self._ds9 = DS9Normalize(clip_lo = 1, clip_hi = 99)
        self._ds9interact = DS9Interacter()
        self._ds9interact.ds9norm_factory(self.im_ax, self._ds9)
        self._zp = ZoomPan()
        self._zpzoom = self._zp.zoom_factory(self.im_ax)
        self._zppan = self._zp.pan_factory(self.im_ax)
        self.im_fig.canvas.draw()

        #File List View
        self.list_view = QtWidgets.QListWidget()
        #list_view.currentItemChanged.connect(self._on_item_clicked)
        self.list_view.setSizePolicy(QtWidgets.QSizePolicy.Minimum,
                                     QtWidgets.QSizePolicy.Expanding)
        list_model = QtGui.QStandardItemModel(self.list_view)
        self.list_view.currentItemChanged.connect(self._on_item_clicked)
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
        self.main_grid.addWidget(self.list_view, 1, 0, 1, 1)
        self.main_grid.addWidget(list_widget, 2, 0, 1, 1)
        self.main_grid.addWidget(image_widget, 1, 1, 2, 1)
        self.main_grid.addWidget(res_widget, 1, 2, 2, 1)
        self.setCentralWidget(self.main_widget)

    def _clear_list(self):
        '''Clear the images list'''
        self.list_view.clear()

    def _open_files(self):
        '''Open the list of files'''
        files = QtWidgets.QFileDialog.getOpenFileNames(self, 'Select the fits files')[0]
        for i in files:
            self.list_view.addItem(i)

    def _remove_files(self):
        '''Remove the selected files'''
        for i in self.list_view.selectedItems():
            self.list_view.takeItem(self.list_view.row(i))

    def _on_item_clicked(self, item, previous):
        self._try_plot(item.text())

    def _on_set_star_clicked(self, pushed):
        '''Click on set star button'''
        canvas = self.im_fig.canvas
        if pushed:
            self._x1, self._y1, self._x2, self._y2 = (None, None, None, None)
            clear_ax_by_label(self.im_ax, 'picked')
            clear_ax_by_label(self.im_ax, 'pair')
            self.im_fig.canvas.draw()
            self._event_pair = canvas.mpl_connect('button_press_event',
                                                  self._mpl_select_points)
        else:
            try:
                canvas.mpl_disconnect(self._event_pair)
            except:
                raise

    def _mpl_select_points(self, event):
        '''mpl event to select stars'''
        if event.button != 1 or event.inaxes != self.im_ax:
            return
        if self.opened is None:
            self.menu_seto.setChecked(False)
            return
        if self._x1 is None or self._y1 is None:
            self._x1, self._y1 = event.xdata, event.ydata
            self.im_ax.plot(self._x1, self._y1, 'yo', ms=15, label='picked')
            self.im_fig.canvas.draw()
        elif self._x2 is None or self._y2 is None:
            self._x2, self._y2 = event.xdata, event.ydata
            self.menu_seto.setChecked(False)
            self.im_ax.plot(self._x2, self._y2, 'yo', ms=15, label='picked')
            self.im_ax.plot((self._x1, self._x2), (self._y1, self._y2),
                            'y-', lw=2.0, label='pair', picker=5)
            self.im_fig.canvas.draw()

    def _on_autoprocess_clicked(self):
        '''Click on autoprocess button'''
        clear_ax_by_label(self.im_ax, 'picked')
        clear_ax_by_label(self.im_ax, 'pair')

    def _try_plot(self, file):
        try:
            self.opened = fits.open(file)[0]
            self.im_ax.imshow(self.opened.data, origin='lower', norm=self._ds9,
                              cmap='gray')
            self._ds9.update_clip(self.opened.data)
        except:
            self.opened = None
            self.im_ax.cla()
        self.im_fig.canvas.draw()

    def _auto_process(self):
        '''Automatic process all stars in the field.'''

    def _process_star(self, x1, y1, x2, y2):
        '''Process one star in (x1, y1) and (x2, y2)'''


def main(argv):
    app = QtWidgets.QApplication(argv)
    mw = CamPolMW()
    mw.show()

    sys.exit(app.exec_())

if __name__ == '__main__':
    main(sys.argv)
