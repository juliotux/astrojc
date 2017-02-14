try:
    from qtpy.QtCore import Qt
    from qtpy.QtWidgets import (QApplication, QComboBox, QDialog, QGridLayout,
                                QGroupBox, QHBoxLayout, QLabel, QLineEdit,
                                QPushButton, QSpinBox, QTextEdit,
                                QVBoxLayout, QWidget, QFileDialog)
except:
    try:
        from PyQt5.QtCore import Qt
        from PyQt5.QtWidgets import (QApplication, QComboBox, QDialog, QGridLayout,
                                QGroupBox, QHBoxLayout, QLabel, QLineEdit,
                                QPushButton, QSpinBox, QTextEdit,
                                QVBoxLayout, QWidget, QFileDialog)
    except:
        try:
            from PyQt4.QtCore import Qt
            from PyQt4.QtWidgets import (QApplication, QComboBox, QDialog, QGridLayout,
                                QGroupBox, QHBoxLayout, QLabel, QLineEdit,
                                QPushButton, QSpinBox, QTextEdit,
                                QVBoxLayout, QWidget, QFileDialog)
        except:
            raise ImportError('None build of Qt found.')

allowed_types = ['String', 'Integer', 'Decimal']

class BatchHeadEditMW(QDialog):
    def __init__(self, parent=None):
        super(QDialog, self).__init__(parent)

        self.file_list = []
        self.key = ''
        self.value = ''
        self.type = None
        self.comment = ''
        self.hdu = 0

        self.create_select_layout()
        self.create_key_layout()
        self.create_run_layout()

        self.log_output = QTextEdit()
        self.log_output.setFixedSize(550, 150)

        self.main_layout = QVBoxLayout()
        warning = QLabel('BE SURE TO BACKUP YOUR DATA BEFORE EDIT!')
        warning.setAlignment(Qt.AlignCenter)
        self.main_layout.addWidget(warning)
        self.main_layout.addWidget(self.select_file)
        self.main_layout.addWidget(self.key_widget)
        self.main_layout.addWidget(self.run_widget)
        self.main_layout.addWidget(self.log_output)

        self.setLayout(self.main_layout)
        self.setWindowTitle("Batch Fits Header Editor")

    def log(self, text):
        self.log_output.append(str(text))

    def create_select_layout(self):
        self.select_file = QGroupBox('Select Files')
        layout = QGridLayout()

        self.select_files_button = QPushButton('Select Fits Files')
        self.select_files_button.clicked.connect(self.select_fitses)
        self.open_list_file = QPushButton('Open .txt List File')
        self.open_list_file.clicked.connect(self.open_list)
        self.hdu_number = QSpinBox()
        self.hdu_number.setValue(0)
        self.hdu_number.setMinimum(0)
        self.hdu_number.setMaximum(100000)

        self.list_edit = QTextEdit()

        layout.addWidget(QLabel('List of FITS files to edit:'), 0, 2)
        layout.addWidget(self.select_files_button, 1, 0, 2, 2)
        layout.addWidget(self.open_list_file, 3, 0, 2, 2)
        layout.addWidget(QLabel('HDU'), 5, 0)
        layout.addWidget(self.hdu_number, 5, 1)
        layout.addWidget(self.list_edit, 1, 2, 6, 6)

        self.select_file.setLayout(layout)

    def create_key_layout(self):
        self.key_widget = QGroupBox('Header Keys')
        layout = QGridLayout()

        layout.addWidget(QLabel('Header Key'), 0, 0)
        layout.addWidget(QLabel('Key Value'), 1, 0)
        layout.addWidget(QLabel('Value Type'), 2, 0)
        layout.addWidget(QLabel('Comment'), 3, 0)

        self.header_key = QLineEdit()
        self.key_value = QLineEdit()
        self.key_type = QComboBox()
        self.key_type.addItems(allowed_types)
        self.key_comment = QLineEdit()

        layout.addWidget(self.header_key, 0, 1)
        layout.addWidget(self.key_value, 1, 1)
        layout.addWidget(self.key_type, 2, 1)
        layout.addWidget(self.key_comment, 3, 1)

        self.key_widget.setLayout(layout)

    def create_run_layout(self):
        self.run_widget = QWidget()
        layout = QHBoxLayout()

        self.run_button = QPushButton('Run')
        self.run_button.clicked.connect(self.run)
        self.clear_button = QPushButton('Clear')
        self.clear_button.clicked.connect(self.clear)

        layout.addWidget(self.clear_button)
        layout.addWidget(self.run_button)

        self.run_widget.setLayout(layout)

    def select_fitses(self):
        files =  QFileDialog.getOpenFileNames(self, 'Select Fits files')
        for path in files[0]:
            self.list_edit.append(path)
        self.log('Added %i files to list.\n---------------' % len(files))

    def open_list(self):
        file = QFileDialog.getOpenFileName(self, 'Select list file')
        file = file[0]
        try:
            self.list_edit.append(open(file, 'r').read())
        except Exception as e:
            self.log('ERROR: List not loaded.\n     - %s' % e)

    def get_values_from_fields(self):
        self.file_list = [i.strip(' ') for i in self.list_edit.toPlainText().split('\n') if i.strip(' ') != '']
        self.key = self.header_key.text()
        self.value = self.key_value.text()
        self.type = str(self.key_type.currentText())
        self.comment = self.key_comment.text()
        self.hdu = self.hdu_number.value()

    def run(self):
        self.get_values_from_fields()
        if len(self.file_list) > 0:
            if self.key.strip(' ') == '':
                self.log('ERROR: No header key name specified.\n---------------')
                return
            if self.value.strip(' ') == '':
                self.log('ERROR: No value specified.\n---------------')
                return
            if self.type not in allowed_types:
                self.log('ERROR: Invalid value type.\n---------------')
        else:
            self.log('ERROR: No files to process!\n---------------')
            return

        value = None
        if self.type == 'Integer':
            try:
                value = int(self.value)
            except Exception as e:
                self.log('ERROR: Value error:\n    - %s' % e)
        elif self.type == 'Decimal':
            try:
                value = float(self.value)
            except Exception as e:
                self.log('ERROR: Value error:\n    - %s' % e)
        else:
            value = str(self.value)

        for i in self.file_list:
            try:
                self.log('File: %s\n    - Key: %s    HDU: %i    Value: %s    Type: %s' % (i, self.hdu, self.key, value, self.type))
                fits.setval(i, self.key, value=value, comment=self.comment, ext=self.hdu)
            except Exception as e:
                self.log('ERROR: Problems in write key:\n    - %s\n---------------' % e)
        self.log('Finished\n---------------')

    def clear(self):
        self.file_list = []
        self.key = ''
        self.value = ''
        self.type = None
        self.comment = ''
        self.hdu = 0
        self.list_edit.setText('')
        self.header_key.setText('')
        self.key_value.setText('')
        self.key_comment.setText('')
        self.hdu_number.setValue(0)

if __name__ == '__main__':

    import sys

    app = QApplication(sys.argv)
    dialog = BatchHeadEditMW()
    sys.exit(dialog.exec_())
