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
#from astrojc.gui.qt_fitsheadertable import HeaderTableModel
from astrojc.mpl.ds9norm import DS9Normalize, DS9Interacter
from astrojc.mpl.zoompan import ZoomPan
from astrojc.logging import log

class FitsViewerMW(QtWidgets.QMainWindow):
    def __init__(self, dir=None, parent=None):
        super(QtWidgets.QMainWindow, self).__init__(parent)
        self.setWindowTitle('Quick Fits Viewer')
        self.setMinimumSize(400, 400)
        self.actual_index = 0
        self.files = []
        self.hdu = 0
        self.opened = None

        self.create_figure()
        self.create_layout()

        self._ds9 = DS9Normalize(clip_lo = 1, clip_hi=99)
        self._ds9interact = DS9Interacter()
        self._ds9interact.ds9norm_factory(self.ax, self._ds9)

        self._zp = ZoomPan()
        self._zpzoom = self._zp.zoom_factory(self.ax)
        self._zppan = self._zp.pan_factory(self.ax)

        if dir is not None:
            self._read_folder(dir)

    def create_layout(self):
        self.main_widget = QtWidgets.QWidget()
        self.main_layout = QtWidgets.QGridLayout(self.main_widget)

        self.open_folder = QtWidgets.QPushButton()
        self.open_folder.setText('Open\nFolder')
        self.open_folder.clicked.connect(self._on_folder_clicked)
        self.open_file = QtWidgets.QPushButton()
        self.open_file.setText('Open\nFile')
        self.open_file.clicked.connect(self._on_file_clicked)

        self.list_view = QtWidgets.QListWidget()
        self.list_view.addItems(self.files)
        self.list_view.currentItemChanged.connect(self._on_item_clicked)
        self.list_view.setSizePolicy(QtWidgets.QSizePolicy.Minimum,
                                     QtWidgets.QSizePolicy.Minimum)

        self.list_widget = QtWidgets.QWidget()
        self.list_buttons = QtWidgets.QHBoxLayout(self.list_widget)
        self.list_buttons.addWidget(self.open_file)
        self.list_buttons.addWidget(self.open_folder)

        self.image_view = FigureCanvas(self.fig)
        FigureCanvas.setSizePolicy(self.image_view,
                                   QtWidgets.QSizePolicy.Preferred,
                                   QtWidgets.QSizePolicy.Preferred)
        self.image_view.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.image_view.setFocus()

        self.previous = QtWidgets.QPushButton()
        self.previous.setText('Previous\nFile')
        self.previous.clicked.connect(self._on_previous_clicked)
        self.next = QtWidgets.QPushButton()
        self.next.setText('Next\nFile')
        self.next.clicked.connect(self._on_next_clicked)

        self.image_buttons = QtWidgets.QWidget()
        self.image_buttons.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                         QtWidgets.QSizePolicy.Fixed)
        self.image_toolbar = QtWidgets.QHBoxLayout(self.image_buttons)
        self.image_toolbar.addWidget(self.previous)
        self.image_toolbar.addWidget(new_spacer())
        self.image_toolbar.addWidget(self.next)
        self.sidebar = QtWidgets.QTabWidget()
        self.header_view = QtWidgets.QTextEdit()
        self.header_view.setReadOnly(True)
        self.header_view.setSizePolicy(QtWidgets.QSizePolicy.Minimum,
                                       QtWidgets.QSizePolicy.Minimum)
        self.sidebar.addTab(self.header_view, 'Header')

        self.main_layout.addWidget(self.list_widget, 0, 0, 1, 1)
        self.main_layout.addWidget(self.list_view, 1, 0, 1, 1)
        self.main_layout.addWidget(self.image_buttons, 0, 1, 1, 1)
        self.main_layout.addWidget(self.image_view, 1, 1, 1, 1)
        self.main_layout.addWidget(self.sidebar, 0, 2, 2, 1)

        self.main_widget.setLayout(self.main_layout)
        self.setCentralWidget(self.main_widget)

    def create_figure(self):
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111)

    def _on_next_clicked(self):
        row = self.actual_index + 1
        self._try_plot(self.files[row])
        self.actual_index = row
        self.list_view.setCurrentRow(row)

    def _on_previous_clicked(self):
        row = self.actual_index - 1
        self._try_plot(self.files[row])
        self.actual_index = row
        self.list_view.setCurrentRow(row)

    def _on_item_clicked(self, item, previous):
        row = self.list_view.row(item)
        self._try_plot(self.files[row])
        self.actual_index = row

    def _on_folder_clicked(self):
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select a folder to open')
        self._read_folder(folder)

    def _on_file_clicked(self):
        file = QtWidgets.QFileDialog.getOpenFileName(self, 'Select a Fits file')[0]
        self._read_folder(os.path.dirname(os.path.abspath(file)))
        self._try_plot(file)

    def _on_view_header(self):
        try:
            header = self.opened[self.hdu].header.tostring(sep='\n')
            self.header_view.setText(str(header))
        except:
            pass

    def _try_plot(self, file):
        try:
            self.opened = fits.open(file)
            self.ax.imshow(self.opened[self.hdu].data, origin='lower', norm=self._ds9, cmap='viridis')
            self._ds9.update_clip(self.opened[self.hdu].data)
            self.setWindowTitle(file)
            self._on_view_header()
        except:
            self.opened = None
        self.fig.canvas.draw()

    def _read_folder(self, folder):
        self.files = sorted(glob.glob(os.path.join(folder, '*')))
        self.list_view.clear()
        for i in range(len(self.files)):
            if not os.path.isfile(self.files[i]):
                self.files[i] = None
        self.files = list(filter(None.__ne__, self.files))
        files = [None]*len(self.files)
        for i in range(len(self.files)):
            files[i] = os.path.basename(self.files[i])
        self.list_view.addItems(files)

def main(argv):
    app = QtWidgets.QApplication(argv)
    dir = sys.argv[1]
    mw = FitsViewerMW(dir)
    mw.show()

    sys.exit(app.exec_())

if __name__ == '__main__':
    main(sys.argv)
