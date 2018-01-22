from qtpy.QtCore import Qt
from qtpy.QtGui import QIcon
from qtpy.QtWidgets import (QApplication, QDialog, QGridLayout, QLabel,
                            QLineEdit, QPushButton, QSpinBox,
                            QWidget, QFileDialog, QListWidget,
                            QListWidgetItem, QMenu, QTableWidgetItem,
                            QTableWidget, QAbstractItemView,
                            QSizePolicy, QSpacerItem, QHBoxLayout,
                            QVBoxLayout, QCheckBox)
from functools import partial
from os import path
from collections import OrderedDict
import json


class ConfigItemText(QWidget):
    def __init__(self, parent=None, default_name=None, default_value=None):
        super(QWidget, self).__init__(parent=parent)

        self._layout = QHBoxLayout()

        self._field_name = QLineEdit('<field name>')
        self._field_name.setSizePolicy(QSizePolicy.Minimum,
                                       QSizePolicy.Minimum)
        if default_name is not None:
            self._field_name.setText(str(default_name))
        self._field_value = QLineEdit('<field value>')
        if default_value is not None:
            self._field_value.setText(str(default_value))

        self._layout.addWidget(self._field_name)
        self._layout.addWidget(self._field_value)

        self.setLayout(self._layout)

    @property
    def name(self):
        return self._field_name.text()

    @property
    def value(self):
        try:
            f = float(self._field_value.text())
            if f % 1 == 0:
                return int(f)
            else:
                return f
        except ValueError:
            return self._field_value.text()


class ConfigItemList(QWidget):
    def __init__(self, parent=None, default_name=None, default_value=None):
        super(QWidget, self).__init__(parent=parent)

        self._field_name = QLineEdit('<field name>')
        self._field_name.setSizePolicy(QSizePolicy.Minimum,
                                       QSizePolicy.Minimum)
        if default_name is not None:
            self._field_name.setText(str(default_name))
        self._field_value = QListWidget()
        self._field_value.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self._field_value.setSizePolicy(QSizePolicy.MinimumExpanding,
                                        QSizePolicy.Expanding)
        if default_value is not None and isinstance(default_value,
                                                    (list, tuple)):
            for i in default_value:
                self.add_file(i)

        self._dlayout = QVBoxLayout()
        self._dlayout.addWidget(self._field_name)
        self._dlayout.addItem(QSpacerItem(5, 5, QSizePolicy.Minimum,
                                          QSizePolicy.Expanding))

        self._add_button = QPushButton(QIcon.fromTheme('list-add'), '')
        self._add_button.setToolTip('Add values to the list')
        self._add_button.clicked.connect(self.add_value)
        self._fil_button = QPushButton(QIcon.fromTheme('document-new'), '')
        self._fil_button.setToolTip('Add file names to the list')
        self._fil_button.clicked.connect(self.open_files)
        self._rem_button = QPushButton(QIcon.fromTheme('list-remove'), '')
        self._rem_button.setToolTip('Add current row from the list')
        self._rem_button.clicked.connect(self.remove_file)
        self._cle_button = QPushButton(QIcon.fromTheme('edit-clear'), '')
        self._cle_button.setToolTip('Empty the list')
        self._cle_button.clicked.connect(self._field_value.clear)
        self._slayout = QVBoxLayout()
        self._slayout.addWidget(self._add_button)
        self._slayout.addWidget(self._fil_button)
        self._slayout.addWidget(self._rem_button)
        self._slayout.addWidget(self._cle_button)
        self._slayout.addItem(QSpacerItem(5, 5, QSizePolicy.Minimum,
                                          QSizePolicy.Expanding))

        self._layout = QHBoxLayout()
        self._layout.addLayout(self._dlayout)
        self._layout.addWidget(self._field_value)
        self._layout.addLayout(self._slayout)

        self.setLayout(self._layout)

    @property
    def name(self):
        return self._field_name.text()

    @property
    def value(self):
        return [self._field_value.item(i).text() for i in
                range(self._field_value.count())]

    def add_value(self):
        item = QListWidgetItem('<new item, click to edit>')
        item.setFlags(item.flags() | Qt.ItemIsEditable)
        self._field_value.addItem(item)

    def add_file(self, filename):
        self._field_value.addItem(path.basename(filename))

    def remove_file(self):
        items = self._field_value.selectedItems()
        for i in items:
            self._field_value.takeItem(self._field_value.row(i))

    def open_files(self):
        files = QFileDialog.getOpenFileNames(self, 'Select list file')[0]
        for i in files:
            self.add_file(i)


class ConfigItemBool(QWidget):
    def __init__(self, parent=None, default_name=None, default_value=None):
        super(QWidget, self).__init__(parent=parent)

        self._field_name = QLineEdit('<field name>')
        if default_name is not None:
            try:
                self._field_name.setText(default_name)
            except Exception as e:
                pass

        self._field_value = QCheckBox()
        if default_value is not None:
            try:
                self._field_value.setChecked(bool(default_value))
            except Exception as e:
                pass

        self._layout = QHBoxLayout()
        self._layout.addWidget(self._field_name)
        self._layout.addWidget(self._field_value)
        self._layout.addItem(QSpacerItem(5, 5, QSizePolicy.Expanding,
                                         QSizePolicy.Minimum))

        self.setLayout(self._layout)

    @property
    def name(self):
        return self._field_name.text()

    @property
    def value(self):
        return self._field_value.isChecked()


class ConfigEditorWidget(QWidget):
    def __init__(self, parent=None):
        super(QWidget, self).__init__(parent=parent)
        self.setMinimumSize(600, 300)

        self._configlist = None

        self._layout = QGridLayout(self)
        self.setLayout(self._layout)

        self._add_menu = QMenu()
        self._add_text_action = self._add_menu.addAction('Add text or number')
        self._add_text_action.triggered.connect(partial(self._add,
                                                        type='text'))
        self._add_bool_action = self._add_menu.addAction('Add boolean')
        self._add_bool_action.triggered.connect(partial(self._add,
                                                        type='bool'))
        self._add_list_action = self._add_menu.addAction('Add list field')
        self._add_list_action.triggered.connect(partial(self._add,
                                                        type='list'))

        self._add_field_button = QPushButton(QIcon.fromTheme('list-add'), '')
        self._add_field_button.setMenu(self._add_menu)
        self._remove_field_button = QPushButton(QIcon.fromTheme('list-remove'),
                                                '')
        self._remove_field_button.clicked.connect(self._remove)
        self._clear_button = QPushButton(QIcon.fromTheme('edit-clear'), '')
        self._clear_button.clicked.connect(self._clear)
        self._slayout = QVBoxLayout()
        self._slayout.addWidget(self._add_field_button)
        self._slayout.addWidget(self._remove_field_button)
        self._slayout.addWidget(self._clear_button)
        self._slayout.addItem(QSpacerItem(5, 5, QSizePolicy.Minimum,
                                          QSizePolicy.Expanding))

        self._add_config_to_list = QPushButton(QIcon.fromTheme('go-next'),
                                               'Add config to list')
        self._add_config_to_list.clicked.connect(self._finish_config)
        self._name = QLineEdit()
        self._order = QSpinBox()
        self._order.setMinimum(1)
        self._nlayout = QHBoxLayout()
        self._nlayout.addWidget(QLabel('Name'))
        self._nlayout.addWidget(self._name)
        self._nlayout.addWidget(QLabel('Order'))
        self._nlayout.addWidget(self._order)
        self._nlayout.addItem(QSpacerItem(5, 5, QSizePolicy.Expanding,
                                          QSizePolicy.Minimum))
        self._nlayout.addWidget(self._add_config_to_list)

        self._list = QListWidget()

        self._layout.addLayout(self._nlayout, 0, 0, 1, 2)
        self._layout.addWidget(self._list, 1, 1, 1, 1)
        self._layout.addLayout(self._slayout, 1, 0, 1, 1)

    def set_config_list_widget(self, widget):
        if isinstance(widget, ListWidget):
            self._configlist = widget
        else:
            raise ValueError('The given widget is not a ListWidget instace.')

    def load_config(self, order, name, config):
        self._list.clear()
        self._order.setValue(order)
        self._name.setText(name)

        for i, v in config.items():
            if isinstance(v, (list, tuple)):
                self._add('file', i, v)
            elif isinstance(v, bool):
                self._add('bool', i, v)
            else:
                self._add('text', i, v)

    def _add(self, type, default_name=None, default_value=None):
        if type == 'text':
            item = ConfigItemText(default_name=default_name,
                                  default_value=default_value)
        elif type == 'list':
            item = ConfigItemList(default_name=default_name,
                                  default_value=default_value)
        elif type == 'bool':
            item = ConfigItemBool(default_name=default_name,
                                  default_value=default_value)
        else:
            raise ValueError('Unrecognized field type.')

        i = QListWidgetItem(self._list)
        i.setSizeHint(item.sizeHint())
        self._list.addItem(i)
        self._list.setItemWidget(i, item)

        return item

    def _remove(self):
        items = self._list.selectedItems()
        for i in items:
            self._list.takeItem(self._list.row(i))

    def _clear(self):
        self._list.clear()

    def _colect_data(self):
        config = OrderedDict()
        for i in range(self._list.count()):
            wdg = self._list.itemWidget(self._list.item(i))
            config[wdg.name] = wdg.value
        return config

    def _finish_config(self):
        if self._configlist is None:
            raise ValueError('No config list widget set.')

        if self._list.count() == 0:
            return

        order = self._order.value()
        name = self._name.text()
        config = self._colect_data()

        self._configlist.set_config(order, name, config)

        self._order.setValue(order + 1)


class DisplayItemWidget(QTableWidgetItem):
    def __init__(self, name, config):
        super(QTableWidgetItem, self).__init__(name)

        self._name = name
        self._config = config

    @property
    def name(self):
        return self._name

    @property
    def config(self):
        return self._config


class ListWidget(QWidget):
    def __init__(self, parent=None):
        super(QWidget, self).__init__(parent=parent)

        self._editor_widget = None

        self._table = QTableWidget()
        self._table.setColumnCount(1)
        self._table.setRowCount(0)
        self._table.verticalHeader().setVisible(True)
        self._table.horizontalHeader().setVisible(False)
        self._table.horizontalHeader().setStretchLastSection(True)
        self._table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self.setMinimumSize(400, 300)

        self._open_button = QPushButton(QIcon.fromTheme('document-open'),
                                        'Load File')
        self._open_button.clicked.connect(self._load)
        self._save_button = QPushButton(QIcon.fromTheme('document-save'),
                                        'Save as')
        self._save_button.clicked.connect(self._save)
        self._edit_item = QPushButton(QIcon.fromTheme('go-previous'),
                                      'Edit')
        self._edit_item.clicked.connect(self._edit)
        self._clear = QPushButton(QIcon.fromTheme('edit-clear'),
                                  '')
        self._clear.setToolTip('Clear the list')
        self._clear.clicked.connect(self._table.clear)
        self._remove_item = QPushButton(QIcon.fromTheme('list-remove'),
                                        '')
        self._remove_item.setToolTip('Remove the current item from list')
        self._remove_item.clicked.connect(self._remove)
        self._move_up = QPushButton(QIcon.fromTheme('go-up'),
                                    '')
        self._move_up.setToolTip('Move the current item up')
        self._move_up.clicked.connect(partial(self._move_row, direction=-1))
        self._move_down = QPushButton(QIcon.fromTheme('go-down'),
                                      '')
        self._move_down.setToolTip('Move the current item down')
        self._move_down.clicked.connect(partial(self._move_row, direction=1))

        self._nlayout = QHBoxLayout()
        self._nlayout.addWidget(self._edit_item)
        self._nlayout.addItem(QSpacerItem(5, 5, QSizePolicy.Expanding,
                                          QSizePolicy.Minimum))
        self._nlayout.addWidget(self._open_button)
        self._nlayout.addWidget(self._save_button)

        self._slayout = QVBoxLayout()
        self._slayout.addWidget(self._remove_item)
        self._slayout.addWidget(self._clear)
        self._slayout.addWidget(self._move_up)
        self._slayout.addWidget(self._move_down)
        self._slayout.addItem(QSpacerItem(5, 5, QSizePolicy.Minimum,
                                          QSizePolicy.Expanding))

        self._layout = QGridLayout()
        self._layout.addLayout(self._nlayout, 0, 0, 1, 2)
        self._layout.addLayout(self._slayout, 1, 1, 1, 1)
        self._layout.addWidget(self._table, 1, 0, 1, 1)

        self.setLayout(self._layout)

    def get_current_config(self):
        item = self._table.item(self._table.currentRow(), 0)
        if item is None:
            raise ValueError('No current item to return.')
        return self._table.currentRow()+1, item.name, item.config

    def set_editor_widget(self, widget):
        if isinstance(widget, ConfigEditorWidget):
            self._editor_widget = widget
        else:
            raise ValueError('The given widget is not a ConfigEditorWidget'
                             ' instace.')

    def _edit(self):
        try:
            order, name, config = self.get_current_config()
            self._editor_widget.load_config(order, name, config)
        except ValueError:
            return

    def _move_row(self, direction):
        """Direction is +1 down and -1 up"""
        curr = self._table.currentRow()
        s_wdg = self._table.takeItem(curr, 0)
        d_wdg = self._table.takeItem(curr+direction, 0)

        self._table.setItem(curr, 0, d_wdg)
        self._table.setItem(curr+direction, 0, s_wdg)

        self._table.setCurrentItem(s_wdg)

    def _remove(self):
        items = self._table.selectedItems()
        for i in items:
            self._table.removeRow(self._table.row(i))

    def set_config(self, order, name, config):
        if self._table.rowCount() < order:
            self._table.setRowCount(order)
        self._table.setItem(order-1, 0, DisplayItemWidget(name, config))

    def load(self, filename):
        data = json.load(open(filename, 'r'), object_hook=OrderedDict)
        for i, v in data.items():
            if not isinstance(v, (dict, OrderedDict)):
                raise ValueError('This json file does not correspond to '
                                 'a valid config file.')

        self._table.clear()
        self._table.setRowCount(0)

        for i, v in data.items():
            last = self._table.rowCount()
            self.set_config(last+1, i, v)

    def save(self, filename):
        out = OrderedDict()
        for i in range(self._table.rowCount()):
            item = self._table.item(i, 0)
            out[item.name] = item.config

        json.dump(out, open(filename, 'w'))

    def _load(self):
        f = QFileDialog.getOpenFileName()[0]
        self.load(f)

    def _save(self):
        f = QFileDialog.getSaveFileName()[0]
        self.save(f)


class MainWindow(QDialog):
    def __init__(self, parent=None):
        super(QDialog, self).__init__(parent=parent)
        self.init_layout()

        self.list_configs = []

    def init_layout(self):
        self._layout = QGridLayout(self)

        self._config = ConfigEditorWidget()
        self._list = ListWidget()
        self._config.set_config_list_widget(self._list)
        self._list.set_editor_widget(self._config)

        self._layout.addWidget(self._config, 0, 0, 1, 1)
        self._layout.addWidget(self._list, 0, 1, 1, 1)

        self.setLayout(self._layout)


if __name__ == '__main__':
    import sys
    app = QApplication(sys.argv)
    dialog = MainWindow()
    sys.exit(dialog.exec_())
