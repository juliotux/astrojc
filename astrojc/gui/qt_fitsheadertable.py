from qtpy import QtCore, QtWidgets, QtGui
from astropy.io import fits

class HeaderTableModel(QtCore.QAbstractTableModel):
    """
    Class to populate a table view with fits header
    """
    def __init__(self, hdu, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self.hdu = hdu
        self.r = len(hdu.data)
        self.c = 3

    def set_hdu(self, hdu):
        self.beginResetModel()
        self.hdu = hdu
        self.endResetModel()

    def rowCount(self, parent=None):
        return self.r

    def columnCount(self, parent=None):
        return self.c

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if index.isValid():
            if role == QtCore.Qt.DisplayRole:
                key = self.hdu.header.keys()[index.row()]
                if index.column() == 0:
                    return key
                elif index.column() == 1:
                    return self.hdu.header[key]
                elif index.column() == 2:
                    return self.hdu.header.comments[key]
        return None

    def headerData(self, p_int, orientation, role):
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                if p_int == 0:
                    return 'Key'
                if p_int == 1:
                    return 'Value'
                if p_int == 2:
                    return 'Comment'
        return None
