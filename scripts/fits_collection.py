#!/bin/env python3

from qtpy.QtCore import Qt, QAbstractTableModel, QTimer
from qtpy.QtGui import QIcon
from qtpy.QtWidgets import (QApplication, QDialog, QGridLayout,
                            QPushButton, QHeaderView,
                            QWidget, QFileDialog, QTableView,
                            QSizePolicy, QSpacerItem)
from functools import partial
from os import path
from collections import OrderedDict
import json

from astropy.logger import log as logger

from ccdproc import ImageFileCollection


default_keywords = ['OBJECT', 'GAIN', 'LAMINA', 'FILTER', 'CAMGAIN',
                    'OUTPTAMP', 'EXPTIME', 'READTIME', 'IMGRECT', 'HBIN',
                    'VBIN', 'RDNOISE']


class CollectionTableModel(QAbstractTableModel):
    """
    Class to populate a table view with a Fits tables
    """
    def __init__(self, summary=None, parent=None):
        super(QAbstractTableModel, self).__init__(parent=parent)
        if summary is not None:
            self.summ = summary
            self.c = len(self.summ.colnames)
        else:
            self.summ = None

        self.previous_sort = []

    def rowCount(self, parent=None):
        if self.summ is None:
            return 0
        else:
            return len(self.summ)

    def columnCount(self, parent=None):
        if self.summ is None:
            return 0
        else:
            return len(self.summ.colnames)

    def data(self, index, role=Qt.DisplayRole):
        if self.summ is None:
            return None
        if index.isValid():
            if role == Qt.DisplayRole:
                key = self.summ.colnames[index.column()]
                return str(self.summ[key][index.row()])
        return None

    def headerData(self, p_int, orientation, role):
        if self.summ is None:
            return None
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return str(self.summ.colnames[p_int])
            elif orientation == Qt.Vertical:
                return p_int
        return None

    def sort(self, col, order):
        self.layoutAboutToBeChanged.emit()
        key = self.summ.colnames[col]
        if key in self.previous_sort:
            self.previous_sort.remove(key)
        self.previous_sort.insert(0, key)
        self.summ.sort(self.previous_sort)
        if order == Qt.DescendingOrder:
            self.summ.reverse()
        self.layoutChanged.emit()

    def set_data(self, new_summary):
        self.layoutAboutToBeChanged.emit()
        self.summ = new_summary
        self.layoutChanged.emit()


class CollectionTableWidget(QWidget):
    # TODO: Move columns
    # TODO: set default keyowrds
    def __init__(self, summary=None, parent=None):
        super(QWidget, self).__init__(parent=parent)
        self.layout = QGridLayout(self)

        self.table_model = CollectionTableModel(summary, self)
        self.table_view = QTableView()
        self.table_view.setSortingEnabled(True)
        hhead = self.table_view.horizontalHeader()
        hhead.setSectionResizeMode(QHeaderView.Interactive)
        # hhead.setSectionsMovable(True)
        hhead.setStretchLastSection(True)

        self.table_view.setModel(self.table_model)
        self.table_view.setMinimumSize(400, 300)
        self.table_view.resizeColumnsToContents()

        self.layout.addWidget(self.table_view, 0, 0)


class MainWindow(QDialog):
    def __init__(self, parent=None):
        super(QDialog, self).__init__(parent=parent)
        self.init_layout()

    def init_layout(self):
        self._layout = QGridLayout(self)

        self.button = QPushButton(QIcon.fromTheme('folder'), 'open')
        self.button.clicked.connect(self.open_folder)
        self._layout.addWidget(self.button, 0, 0, 1, 1)
        self._layout.addItem(QSpacerItem(200, 1, QSizePolicy.Expanding,
                                         QSizePolicy.Minimum), 0, 1, 1, 1)

        self.table = CollectionTableWidget()
        self._layout.addWidget(self.table, 1, 0, 1, 3)

        self.setLayout(self._layout)

    def open_folder(self):
        folder = QFileDialog.getExistingDirectory()
        coll = ImageFileCollection(folder, keywords=default_keywords)
        self.table.table_model.set_data(coll.summary)


if __name__ == '__main__':
    import sys
    app = QApplication(sys.argv)
    dialog = MainWindow()
    sys.exit(dialog.exec_())
