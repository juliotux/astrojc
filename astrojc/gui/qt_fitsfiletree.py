from PyQt5.QtWidgets import (QWidget, QTreeView, QApplication, QFileSystemModel,
                             QLabel, QLineEdit, QVBoxLayout, QListView)
from PyQt5.QtCore import pyqtSlot, QModelIndex
from PyQt5.QtGui import QIcon, QStandardItemModel, QStandardItem
import os
import sys

from astropy.io import fits

class HDUFileSystemModel(QWidget):
    def __init__(self, path, parent=None):
        super(QWidget, self).__init__(parent)

        self.pathRoot = path
        self.dirmodel = QFileSystemModel(self)
        self.dirmodel.setRootPath(self.pathRoot)
        self.indexRoot = self.dirmodel.index(self.dirmodel.rootPath())

        self.filetree = QTreeView(self)
        self.filetree.setModel(self.dirmodel)
        self.filetree.setRootIndex(self.indexRoot)
        self.filetree.clicked.connect(self.on_filetree_clicked)

        self.hdulist = QListView(self)
        self.hdulistmodel = QStandardItemModel(self)
        self.hdulist.setModel(self.hdulistmodel)
        self.hdulist.clicked.connect(self.on_hdulist_clicked)

        self.layout = QVBoxLayout(self)
        self.layout.addWidget(self.filetree)
        self.layout.addWidget(self.hdulist)

    @pyqtSlot(QModelIndex)
    def on_filetree_clicked(self, index):
        indexItem = self.dirmodel.index(index.row(), 0, index.parent())
        self.hdulistmodel.clear()
        if not self.dirmodel.isDir(indexItem):
            try:
                f = fits.open(self.dirmodel.filePath(indexItem))
                for i in f:
                    self.hdulistmodel.appendRow(self._get_hdu_item(i))
            except:
                pass

    @pyqtSlot(QModelIndex)
    def on_hdulist_clicked(self, index):
        indexItem = self.hdulistmodel.index(index.row(), 0, index.parent())

    def _get_hdu_item(self, hdu):
        name = ('... ' if hdu.name == '' else hdu.name + ' ')
        if isinstance(hdu, fits.PrimaryHDU):
            icon = QIcon.fromTheme('image-x-generic')
            data = name + 'PrimaryHDU (%i axis)' % len(hdu.shape)
        elif isinstance(hdu, fits.PrimaryHDU):
            icon = QIcon.fromTheme('image-x-generic')
            data = name + 'ImageHDU (%i axis)' % len(hdu.shape)
        elif isinstance(hdu, fits.BinTableHDU):
            icon = QIcon.fromTheme('x-office-spreadsheet')
            data = name + 'BinTableHDU'
        elif isinstance(hdu, fits.TableHDU):
            icon = QIcon.fromTheme('x-office-spreadsheet')
            data = name + 'TableHDU'
        else:
            icon = None
            data = 'Unkowed HDU'
        item = QStandardItem(icon, data)
        item.hdu = hdu
        return item

    def set_root_path(self, path):
        self.pathRoot = path
        self.dirmodel.setRootPath(path)
        self.indexRoot = self.dirmodel.index(self.dirmodel.rootPath())
