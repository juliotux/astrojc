from qtpy import QtCore, QtWidgets, QtGui
import os
import sys

from astropy.io import fits

from .qt_helper import new_spacer

hdu_icon_map = {'image': 'image-x-generic',
                'table': 'x-office-spreadsheet',
                'unkown': ''}

class HDUTreeItem(QtWidgets.QWidget):
    def __init__(self, hdu, parent=None):
        super(QtWidgets.QWidget, self).__init__(parent)
        self.hdu = hdu
        self.name = hdu.name
        if self.hdu.__class__.__name__ in ['PrimaryHDU', 'ImageHDU']:
            self.type = 'image'
            try:
                self.shape = self.hdu.data.shape
                self.size_str = str(self.shape[0])
                for i in self.shape[1:]:
                    self.size_str = self.size_str + 'x' + str(i)
            except:
                self.shape = None
                self.size_str = 'empty'
        elif self.hdu.__class__.__name__ in ['TableHDU', 'BinTableHDU']:
            self.type = 'table'
            try:
                self.shape = (len(self.hdu.data), len(self.hdu.data.columns))
                self.size_str = '%i cols x %i rows' % self.shape
            except:
                self.shape = None
                self.size_str = 'empty'
        else:
            self.type = 'unkown'
            self.size_str = 'unkown'
            self.shape = None

        self.create_widget()
        self.plotter = None

    def create_widget(self):
        '''
        Creates the widget layout.
        '''
        l1 = QtWidgets.QHBoxLayout()
        l2 = QtWidgets.QVBoxLayout()
        l3 = QtWidgets.QHBoxLayout()

        name_label = QtWidgets.QLabel(self.name)
        size_label = QtWidgets.QLabel(self.size_str)

        l3.addWidget(name_label)
        l3.addWidget(size_label)

        class_name_label = QtWidgets.QLabel(self.hdu.__class__.__name__)
        l2.addWidget(class_name_label)
        l2.addLayout(l3)

        icon = QtWidgets.QLabel()
        icon.setPixmap(QtGui.QIcon.fromTheme(hdu_icon_map[self.type]).pixmap(48))
        l1.addWidget(icon)
        l1.addLayout(l2)
        l1.addWidget(new_spacer())

        self.open_data_button = QtWidgets.QPushButton()
        self.open_data_button.setText('Open\nData')
        self.open_header_button = QtWidgets.QPushButton()
        self.open_header_button.setText('Open\nHeader')

        l1.addWidget(self.open_data_button)
        l1.addWidget(self.open_header_button)

        self.setLayout(l1)

    def set_plotter(self, plotter):
        '''
        Set the plotter used to show the data.
        '''
        self.plotter = plotter

class HDUFileSystemModel(QtWidgets.QWidget):
    def __init__(self, path, parent=None):
        super(QtWidgets.QWidget, self).__init__(parent)

        self.pathRoot = path
        self.dirmodel = QtWidgets.QFileSystemModel(self)
        self.dirmodel.setRootPath(self.pathRoot)
        self.indexRoot = self.dirmodel.index(self.dirmodel.rootPath())

        self.filetree = QtWidgets.QTreeView(self)
        self.filetree.setModel(self.dirmodel)
        self.filetree.setRootIndex(self.indexRoot)
        self.filetree.clicked.connect(self.on_filetree_clicked)

        self.hdulist = None

        self.layout = QtWidgets.QVBoxLayout(self)
        self.layout.addWidget(self.filetree)
        self.layout.addWidget(self.hdulist)

    def on_filetree_clicked(self, index):
        indexItem = self.dirmodel.index(index.row(), 0, index.parent())
        self.hdulist.clear()
        if not self.dirmodel.isDir(indexItem):
            if hdulist is None:
                raise ValueError('No HDUListView is connected. Use `setHDUListView`'
                                 'To connect one.')
            self.hdulist.add_fitsfile(self.dirmodel.filePath(indexItem))

    def set_root_path(self, path):
        self.pathRoot = path
        self.dirmodel.setRootPath(path)
        self.indexRoot = self.dirmodel.index(self.dirmodel.rootPath())

    def setHDUListView(self, hdulist):
        self.hdulist = hdulist

class HDUListView(QtWidgets.QListWidget):
    def __init__(self, parent=None):
        super(QtWidgets.QListWidget, self).__init__(parent)
        self.plotter = None

    def add_fitsfile(self, hdulist):
        self.clear()
        try:
            f = fits.open()
                for i in f:
                    self.add_hdu(i)
            except Exception as e:
                print(e)
                pass

    def add_hdu(self, hdu):
        newitem = QtWidgets.QListWidgetItem(self)
        hduitem = HDUTreeItem(hdu)
        newitem.setSizeHint(hduitem.sizeHint())
        self.addItem(newitem)
        self.setItemWidget(newitem, hduitem)
        if plotter is not None:
            self.hduitem.set_plotter(self.plotter)

    def set_plotter(self, plotter):
        self.plotter = plotter
