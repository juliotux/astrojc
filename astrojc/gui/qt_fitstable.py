from qtpy import QtCore, QtWidgets, QtGui
import numpy as np

class FITSTableModel(QtCore.QAbstractTableModel):
    """
    Class to populate a table view with a Fits tables
    """
    def __init__(self, hdu, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self.hdu = hdu
        self.r = len(hdu.data)
        self.c = len(hdu.data.columns)

    def rowCount(self, parent=None):
        return self.r

    def columnCount(self, parent=None):
        return self.c

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if index.isValid():
            if role == QtCore.Qt.DisplayRole:
                return str(self.hdu.data[index.row()][index.column()])
        return None

    def headerData(self, p_int, orientation, role):
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                return str(self.hdu.data.columns.names[p_int])
            elif orientation == QtCore.Qt.Vertical:
                return p_int
        return None

class HDUTableWidget(QtWidgets.QTableView):
    def __init__(self, hdu, parent):
        super(QtWidgets.QTableView, self).__init__(parent)

        self.hdu = hdu
        self.table_model = FITSTableModel(hdu, self)

        self.setModel(self.table_model)
        self.setMinimumSize(400, 300)
        self.resizeColumnsToContents()

        #TODO: plot, edit, save ascii, etc...
