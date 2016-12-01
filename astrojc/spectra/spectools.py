import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from scipy.interpolate import splrep,splev

def rv_correct(wave, rv):
    '''
    Correct a wavelenght array with radial velocity, in km/s

    Parameters:
    -----------
    wave : array_like
        Array of wavelength to be corrected by radial velocity.
    rv : float
        Radial velocity of the object.

    Returns:
    --------
    array_like : The radial velocity corrected wavelength array.
    '''
    return wave*(1 - rv/3e5)

class SpecNorm:
    def __init__(self, ax, wave, flux):
        self.ax = ax
        self.flux = flux
        self.flux_norm = None
        self.wave = wave
        self.cont_x = []
        self.cont_y = []
        self.cont_flux = None
        #Draw and connect
        self.ax.plot(self.wave, self.flux, 'k-', label='spec_origin')
        self.connect()
        self.draw()

    def connect(self):
        self.click_event = self.ax.figure.canvas.mpl_connect('button_press_event',
                                                             self.onclick)
        self.pick_event = self.ax.figure.canvas.mpl_connect('pick_event',
                                                            self.onpick)
        self.type_event = self.ax.figure.canvas.mpl_connect('key_press_event',
                                                            self.ontype)

    def disconnect(self):
        self.ax.figure.canvas.mpl_disconnect(self.click_event)
        self.ax.figure.canvas.mpl_disconnect(self.pick_event)
        self.ax.figure.canvas.mpl_disconnect(self.type_event)

    def onpick(self, event):
        if event.mouseevent.button==3:
            if hasattr(event.artist,'get_label') and event.artist.get_label()=='cont_pnt':
                event.artist.remove()
        self.draw()

    def onclick(self, event):
        toolbar = self.ax.figure.canvas.toolbar
        if event.button == 1 and toolbar.mode=='':
            y = event.ydata
            self.ax.plot(event.xdata, y, 'rs', ms=10,
                         picker=5, label='cont_pnt')
            self.test.append(event)
        self.draw()

    def ontype(self, event):
        # 'enter' : Calculate the spline continuum
        # 'n' : Normalize the spectrum with the generated continuum
        # 'r' : Clear the axes and plot the original spectrum
        if event.key=='enter':
            cont_pnt_coord = []
            for artist in self.ax.get_children():
                if hasattr(artist,'get_label') and artist.get_label()=='cont_pnt':
                    cont_pnt_coord.append(artist.get_data())
                elif hasattr(artist,'get_label') and artist.get_label()=='continuum':
                    artist.remove()
            cont_pnt_coord = np.array(cont_pnt_coord)[..., 0]
            sort_array = np.argsort(cont_pnt_coord[:, 0])
            self.cont_x, self.cont_y = cont_pnt_coord[sort_array].T
            self.spline = splrep(self.cont_x, self.cont_y, k=3)
            self.cont_flux = splev(self.wave, self.spline)
            self.ax.plot(self.wave, self.cont_flux, 'r-',
                         lw=2, label='continuum')

        elif event.key=='n':
            if self.cont_flux is not None:
                self.ax.cla()
                self.flux_norm = self.flux/self.cont_flux
                self.ax.plot(self.wave, self.flux_norm,
                             'k-', label='spec_norm')

        elif event.key=='r':
            self.ax.cla()
            self.ax.plot(self.wave, self.flux, 'k-', label='spec_origin')

        self.draw()

    def draw(self):
        self.ax.figure.canvas.draw()

ll_dt = np.dtype([('lamb_0', 'f8'), ('elem', 'U2'), ('ion', 'U8'), ('mult', 'U12')])

class LineIdent:
    def __init__(self, ax, wave, flux, linelist):
        self.lamb_limit = 0.2 #limit of lambda detection
        self.wave = wave
        self.flux = flux
        self.ax = ax
        self.linelist = linelist
        self.lines = []
        self.active = []
        self.labels = []
        self.ax.plot(self.wave, self.flux, 'k-', label='spec_origin')
        self.connect()
        self.draw()

    def connect(self):
        self.click_event = self.ax.figure.canvas.mpl_connect('button_press_event',
                                                             self.onclick)
        self.pick_event = self.ax.figure.canvas.mpl_connect('pick_event',
                                                            self.onpick)
        self.type_event = self.ax.figure.canvas.mpl_connect('key_press_event',
                                                            self.ontype)

    def disconnect(self):
        self.ax.figure.canvas.mpl_disconnect(self.click_event)
        self.ax.figure.canvas.mpl_disconnect(self.pick_event)
        self.ax.figure.canvas.mpl_disconnect(self.type_event)

    def onpick(self, event):
        if event.mouseevent.button == 1:
            if hasattr(event.artist,'lineid'):
                if event.artist not in self.active:
                    event.artist.set_fontweight('bold')
                    self.active.append(event.artist)
        elif event.mouseevent.button == 3:
            if hasattr(event.artist,'lineid'):
                if event.artist in self.active:
                    event.artist.set_fontweight('normal')
                    self.active.remove(event.artist)
        self.draw()

    def onclick(self, event):
        self.draw()

    def ontype(self, event):
        # 'enter' : show the possible identifications of the line
        # 'a' : add the active lines to `lines` list
        if event.key == 'enter':
            for artist in self.labels:
                artist.remove()
            self.labels = []
            for i in self.linelist:
                if abs(i['lamb_0'] - event.xdata) < self.lamb_limit:
                    y_f = np.interp(i['lamb_0'], self.wave, self.flux)
                    l = self.ax.annotate('%.1f %s %s (m%s)' % (i['lamb_0'], i['elem'], i['ion'], i['mult']),
                                         xy = (i['lamb_0'], y_f), xytext = (i['lamb_0'], (1-0.05)*y_f),
                                         rotation=90, picker=5)
                    setattr(l, 'lineid', i)
                    self.labels.append(l)

        elif event.key == 'a':
            self.lines.append([i.lineid for i in self.active])
            self.active = []
            for i in self.labels:
                i.set_fontweight('light')

        self.draw()

    def draw(self):
        self.ax.figure.canvas.draw()

def read_illss(fname):
    f = open(fname, 'r')
    f = f.readlines()
    line = [None]*len(f)
    elem = [None]*len(f)
    ion = [None]*len(f)
    mult = [None]*len(f)
    for i in range(len(f)):
        line[i] = float(f[i][:12])
        elem[i] = f[i][12:14]
        ion[i] = f[i][14:21].replace(' ', '')
        mult[i] = f[i][29:41].replace(' ', '').replace('\n', '')
    return np.array(list(zip(line, elem, ion, mult)), dtype=ll_dt)
