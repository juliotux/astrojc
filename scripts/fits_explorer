#!/bin/env python3

import os
import sys
sys.path.append('./build/lib/')

from qtpy import QtCore, QtWidgets
from astrojc.gui.qt_2dimage_plotter import HDUFigureCanvas2D
from astrojc.gui.qt_fitsfiletree import HDUFileSystemModel, HDUListView
from astrojc.gui.qt_fitstable import HDUTableWidget

from astropy.io import fits

class ViewWidget(QtWidgets.QTabWidget):
    def __init__(self, parent=None):
        super(QtWidgets.QTabWidget, self).__init__(parent)
        self.setTabsClosable(True)

    def add_hdu_tab(self, hdu, name):
        if isinstance(hdu, (fits.PrimaryHDU, fits.ImageHDU)):
            self.addTab(HDUFigureCanvas2D(hdu, self), name)
        elif isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
            self.addTab(HDUTableWidget(hdu, self), name)

class FitsExplorerMW(QtWidgets.QMainWindow):
    def __init__(self, root='', parent=None):
        super(QtWidgets.QMainWindow, self).__init__(parent)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle('FitsExplorer')

        self.setMinimumSize(800, 600)

        self.create_fits_viewer()
        self.create_file_tree(root)

    def create_file_tree(self, root=''):
        self.file_tree_dock = QtWidgets.QDockWidget('File Tree', self)
        self.file_tree_dock.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea | QtCore.Qt.RightDockWidgetArea)
        self.hdulist_dock = QtWidgets.QDockWidget('HDUList', self)
        self.hdulist_dock.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea | QtCore.Qt.RightDockWidgetArea)

        self.hdulist_view = HDUListView(self.hdulist_dock)
        self.hdulist_view.set_open_data(self.open_data)
        self.file_tree = HDUFileSystemModel(root, self.file_tree_dock)
        self.file_tree.setHDUListView(self.hdulist_view)

        self.file_tree_dock.setWidget(self.file_tree)
        self.hdulist_dock.setWidget(self.hdulist_view)

        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.file_tree_dock)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, self.hdulist_dock)

    def create_fits_viewer(self):
        self.fits_viewer = ViewWidget(self)
        self.setCentralWidget(self.fits_viewer)

    def open_data(self, hdu):
        self.fits_viewer.add_hdu_tab(hdu, 'testing data')

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName('FitsExplorer')

    mw = FitsExplorerMW()

    mw.show()

    sys.exit(app.exec_())
