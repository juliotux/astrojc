#!/bin/env python3

import os
import sys

import matplotlib
matplotlib.use('Qt5Agg')
from qtpy import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import colors
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy as np
import glob

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from astrojc.gui.qt_helper import new_spacer
from astrojc.mpl.ds9norm import DS9Normalize, DS9Interacter
from astrojc.mpl.zoompan import ZoomPan
from astrojc.logging import log
from astrojc.math.array import trim_array
from astrojc.reducing import polarimetry as pol

#from photutils import DAOStarFinder, CircularAperture, CircularAnnulus
#from photutils import DAOPhotPSFPhotometry, aperture_photometry
import sep

from scipy.spatial import cKDTree


def clear_ax_by_label(ax, label=None):
    '''Clear all the artists with a given label from an axes.'''
    for i in ax.get_lines():
        if i.get_label() == label or label is None:
            i.remove()
    for i in ax.patches:
        if i.get_label() == label or label is None:
            i.remove()

def find_source(data, x, y, snr, box_size=None, fwhm=5, **kwargs):
    x0, y0 = x, y
    if box_size is not None:
        data, x, y = trim_array(data, box_size, (x, y))
    #med, mean, std = sigma_clipped_stats(data)
    #star_finder = DAOStarFinder(threshold=med + snr*std, fwhm=fwhm, **kwargs)
    #sources = star_finder.find_stars(data)
    #kdt = cKDTree(list(zip(sources['xcentroid'], sources['ycentroid'])))
    #dx, indx = kdt.query((x, y), 1)
    #nx = sources['xcentroid'][indx]
    #ny = sources['ycentroid'][indx]
    data = np.array(data, '<f8')
    try:
        bkg = sep.Background(data)
    except:
        data = data.byteswap().newbyteorder()
        bkg = sep.Background(data)
    sources = sep.extract(data - bkg.back(), snr, err=bkg.globalrms)
    kdt = cKDTree(list(zip(sources['x'], sources['y'])))
    dx, indx = kdt.query((x, y), 1)
    return (x0 + (sources[indx]['x']-x), y0 + (sources[indx]['y']-y))

def aperture_photometry(data, x, y, r, ann=None, gain=1.0):
    data = np.array(data, '<f8')
    try:
        bkg = sep.Background(data)
    except:
        data = data.byteswap().newbyteorder()
        bkg = sep.Background(data)
    flux, fluxerr, flag = sep.sum_circle(data-bkg.back(), x, y, r, bkgann=ann,
                                         err=bkg.globalrms, gain=gain)
    return flux, fluxerr


class CamPolMW(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(QtWidgets.QMainWindow, self).__init__(parent)
        self.setWindowTitle('Quick Polarimetry')
        self.setMinimumSize(1000, 600)

        self._x1, self._y1 = (None, None)
        self._x2, self._y2 = (None, None)
        self.opened = None

        self._create_layout()

    def _create_layout(self):
        self.main_widget = QtWidgets.QWidget()
        self.main_grid = QtWidgets.QGridLayout(self.main_widget)

        #Main Menu
        menu = QtWidgets.QToolBar()
        menu.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                           QtWidgets.QSizePolicy.Fixed)
        menu_open = menu.addAction('Open')
        menu_open.triggered.connect(self._open_files)
        menu_clear = menu.addAction('Clear')
        menu_clear.triggered.connect(self._clear_list)
        menu.addSeparator()
        menu_auto = menu.addAction('Auto Pol.')
        menu_auto.triggered.connect(self._on_autoprocess_clicked)
        self.menu_seto = menu.addAction('Set Star')
        self.menu_seto.toggled.connect(self._on_set_star_clicked)
        self.menu_seto.setCheckable(True)
        menu.addWidget(new_spacer())
        menu_seto = menu.addAction('Config.')

        #Display image figure
        self.im_fig = Figure()
        self.im_ax = self.im_fig.add_subplot(111)
        image_widget = FigureCanvas(self.im_fig)
        self._ds9 = DS9Normalize(clip_lo = 1, clip_hi = 99)
        self._ds9interact = DS9Interacter()
        self._ds9interact.ds9norm_factory(self.im_ax, self._ds9)
        self._zp = ZoomPan()
        self._zpzoom = self._zp.zoom_factory(self.im_ax)
        self._zppan = self._zp.pan_factory(self.im_ax)
        self.im_fig.canvas.draw()

        #File List View
        self.list_view = QtWidgets.QListWidget()
        #list_view.currentItemChanged.connect(self._on_item_clicked)
        self.list_view.setSizePolicy(QtWidgets.QSizePolicy.Minimum,
                                     QtWidgets.QSizePolicy.Expanding)
        list_model = QtGui.QStandardItemModel(self.list_view)
        self.list_view.currentItemChanged.connect(self._on_item_clicked)
        list_widget = QtWidgets.QWidget()
        list_buttons = QtWidgets.QHBoxLayout(list_widget)
        add_button = QtWidgets.QPushButton()
        add_button.setText('+')
        add_button.clicked.connect(self._open_files)
        remove_button = QtWidgets.QPushButton()
        remove_button.setText('-')
        remove_button.clicked.connect(self._remove_files)
        clear_button = QtWidgets.QPushButton()
        clear_button.setText('clear')
        clear_button.clicked.connect(self._clear_list)
        list_buttons.addWidget(add_button)
        list_buttons.addWidget(remove_button)
        list_buttons.addWidget(clear_button)

        #Plot Results
        self.res_fig = Figure(figsize=(4, 3))
        self.res_ax = self.res_fig.add_subplot(111)
        self.res_ax.set_xlim(0, 360)
        self.res_ax.set_xlabel('Retarder Position')
        self.res_ax.set_ylabel('Flux Ratio')
        self.res_fig.tight_layout()
        res_widget = QtWidgets.QWidget()
        res_widget.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                 QtWidgets.QSizePolicy.Minimum)
        res_layout = QtWidgets.QVBoxLayout(res_widget)
        self.res_fig_widget = FigureCanvas(self.res_fig)
        self.res_fig_widget.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                          QtWidgets.QSizePolicy.Minimum)
        self.res_table = QtWidgets.QTableWidget()
        res_layout.addWidget(self.res_fig_widget)
        res_layout.addWidget(self.res_table)

        # Grid positioning
        self.main_grid.addWidget(menu, 0, 0, 1, 3)
        self.main_grid.addWidget(self.list_view, 1, 0, 1, 1)
        self.main_grid.addWidget(list_widget, 2, 0, 1, 1)
        self.main_grid.addWidget(image_widget, 1, 1, 2, 1)
        self.main_grid.addWidget(res_widget, 1, 2, 2, 1)
        self.setCentralWidget(self.main_widget)

    #These properties will be put in a config window
    @property
    def box_size(self):
        return 40

    @property
    def photometry(self):
        return 'aperture'

    @property
    def apertures(self):
        return np.array(range(10)) + 1

    @property
    def sky_radius(self):
        return (15, 20)

    @property
    def angle_unit(self):
        return 'degree'

    @property
    def position_key(self):
        return 'LAM-POS'

    @property
    def angle_rotation(self):
        return 22.5

    @property
    def angle_key(self):
        return None

    @property
    def find_fwhm(self):
        return 5

    @property
    def min_snr(self):
        return 3

    @property
    def hdu(self):
        return 0

    def _clear_list(self):
        '''Clear the images list'''
        self.list_view.clear()

    def _open_files(self):
        '''Open the list of files'''
        files = QtWidgets.QFileDialog.getOpenFileNames(self, 'Select the fits files')[0]
        for i in files:
            self.list_view.addItem(i)

    def _remove_files(self):
        '''Remove the selected files'''
        for i in self.list_view.selectedItems():
            self.list_view.takeItem(self.list_view.row(i))

    def _on_item_clicked(self, item, previous):
        try:
            f = item.text()
        except:
            f = ''
        self._try_plot(f)

    def _on_set_star_clicked(self, pushed):
        '''Click on set star button'''
        canvas = self.im_fig.canvas
        if pushed:
            self._x1, self._y1, self._x2, self._y2 = (None, None, None, None)
            clear_ax_by_label(self.im_ax)
            self.im_fig.canvas.draw()
            self._mpl_pair = canvas.mpl_connect('button_press_event',
                                                  self._mpl_select_points)
        else:
            try:
                canvas.mpl_disconnect(self._mpl_pair)
            except:
                raise

    def _mpl_select_points(self, event):
        '''mpl event to select stars'''
        if event.button != 1 or event.inaxes != self.im_ax:
            return
        if self.opened is None:
            self.menu_seto.setChecked(False)
            return
        if self._x1 is None or self._y1 is None:
            x, y = event.xdata, event.ydata
            self._x1, self._y1 = find_source(self.opened.data, x, y,
                                             self.min_snr,
                                             fwhm=self.find_fwhm)
            self.im_ax.plot(self._x1, self._y1, 'yo', ms=15, label='picked')
            self.im_fig.canvas.draw()
        elif self._x2 is None or self._y2 is None:
            x, y = event.xdata, event.ydata
            self._x2, self._y2 = find_source(self.opened.data, x, y,
                                             self.min_snr,
                                             fwhm=self.find_fwhm)
            self.menu_seto.setChecked(False)
            self.im_ax.plot(self._x2, self._y2, 'yo', ms=15, label='picked')
            self.im_ax.plot((self._x1, self._x2), (self._y1, self._y2),
                            'y-', lw=2.0, label='pair', picker=5)
            self._mpl_pick = self.im_fig.canvas.mpl_connect('pick_event',
                                                            self._mpl_pick_pair)
            self.im_fig.canvas.draw()

    def _mpl_pick_pair(self, event):
        mouse = event.mouseevent
        if mouse.button != 1 or mouse.inaxes != self.im_ax:
            return
        if self.opened is None:
            return
        line = event.artist
        (x1, x2), (y1, y2) = line.get_xdata(), line.get_ydata()
        self._process_star(x1, y1, x2, y2)


    def _on_autoprocess_clicked(self):
        '''Click on autoprocess button'''
        clear_ax_by_label(self.im_ax, 'picked')
        clear_ax_by_label(self.im_ax, 'pair')

    def _try_plot(self, file):
        try:
            self.opened = fits.open(file)[self.hdu]
            self.im_ax.imshow(self.opened.data, origin='lower', norm=self._ds9,
                              cmap='gray')
            self._ds9.update_clip(self.opened.data)
        except:
            self.opened = None
            self.im_ax.cla()
        self.im_fig.canvas.draw()

    def _auto_process(self):
        '''Automatic process all stars in the field.'''

    def _process_star(self, x1, y1, x2, y2):
        '''Process one star in (x1, y1) and (x2, y2)'''
        # flux, flux_error at best aperture
        ord_flux = np.zeros((self.list_view.count(), 2))
        ext_flux = np.zeros((self.list_view.count(), 2))
        lam_pos = np.zeros(self.list_view.count())
        if self.photometry == 'aperture':
            ord_tmp = np.zeros((self.list_view.count(), len(self.apertures), 2))
            ext_tmp = np.zeros((self.list_view.count(), len(self.apertures), 2))
            for i in range(self.list_view.count()):
                name = self.list_view.item(i).text()
                print('processing {}'.format(name))
                hdu = fits.open(name)[self.hdu]
                d = hdu.data
                for j in range(len(self.apertures)):
                    f, f_er = aperture_photometry(d, (x1, x2), (y1, y2),
                                                  self.apertures[j],
                                                  ann=self.sky_radius)
                    ord_tmp[i, j, 0] = f[0]
                    ord_tmp[i, j, 1] = f_er[0]
                    ext_tmp[i, j, 0] = f[1]
                    ext_tmp[i, j, 1] = f_er[1]

                if self.angle_key is not None:
                    lam_pos[i] = float(hdu.header[self.angle_key])
                else:
                    lam_pos[i] = float(hdu.header[self.position_key])*self.angle_rotation

            snr = ord_tmp[:, :, 0]/ord_tmp[:, :, 1]
            ap = np.argmax(np.sum(snr, axis=0))
            print('Selected aperture: {} px'.format(self.apertures[ap]))

            ord_flux = ord_tmp[:, ap, :]
            ext_flux = ext_tmp[:, ap, :]

        (q, u), errors = pol.calculate_polarimetry_half(ord_flux[:, 0], ext_flux[:, 0], np.deg2rad(lam_pos))

        clear_ax_by_label(self.res_ax)
        self.res_ax.plot(lam_pos, (ord_flux[:,0] - ext_flux[:,0])/(ord_flux[:,0] + ext_flux[:,0]), 'bo')
        arr = np.linspace(0, 360, 1000)
        self.res_ax.plot(arr, pol._half(np.deg2rad(arr), q, u), 'r-')
        self.res_fig.canvas.draw()


def main(argv):
    app = QtWidgets.QApplication(argv)
    mw = CamPolMW()
    mw.show()

    sys.exit(app.exec_())

if __name__ == '__main__':
    main(sys.argv)
