from qtpy import QtCore, QtWidgets, QtGui

def new_spacer(parent = None):
    spacer = QtWidgets.QWidget(parent)
    spacer.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
    return spacer
