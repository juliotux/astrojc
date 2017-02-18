'''
New Python implementation of IRAF imexam functionalities.

Includes a new key set to fits explorer.
'''

from astropy.modeling import models, fitting
from matplotlib import pyplot as plt
import numpy as np

from scipy.kdtree import cKDTree

from ..logging import log
from ..signal import MySignal
from ..math.array import trim_array

class ImexamConfig():
    def __init__(self):
        self.box_size = 15   #Box size for contour and surface plots and fitting

        self.n_contours = 10 #Number of contours to be draw

        self.cmap = 'viridis'#Cmap to use in plots

        self.hist_bins = 12  #Number of bins of the histogram plot


class Imexam():
    '''
    Class to mimic imexam functionalities.
    '''
    def __init__(self):
        self.commands = {'a': (self.aper_phot, 'Aperture sum, with radius region_size '),
                         'j': (self.line_fit, '1D [Gaussian1D default] line fit '),
                         'k': (self.column_fit, '1D [Gaussian1D default] column fit'),
                         'm': (self.report_stat, 'Square region stats, in [region_size],default is median'),
                         'l': (self.plot_line, 'Return line plot'),
                         'c': (self.plot_column, 'Return column plot'),
                         'r': (self.radial_profile, 'Return the radial profile plot'),
                         'h': (self.histogram, 'Return a histogram in the region around the cursor'),
                         'e': (self.contour, 'Return a contour plot in a region around the cursor'),
                         's': (self.surface, 'Display a surface plot around the cursor location')}

        self.ax = None
        self.fig = None
        self.hdu = None
        self.plot_ax= None
        self.connected = None

        self.print_text = MySignal()

        self.config = ImexamConfig()

    def _check_plot_ax(self):
        if self.plot_ax is None:
            log.error('No plot axes connected.')
        else:
            self.plot_ax.cla()

    def set_hdu(self, hdu):
        self.hdu = hdu
        log.info('HDU setted to %s' % str(self.hdu))

    def connect(self, ax, plot_ax=None):
        '''
        Turn on the imexam in an axes instance.

        Parameters:
        -----------
            ax : matplotlib.Axes
                The axes instance that will be analysed (with the image).
            plot_ax : matplotlib.Axes
                The axes where the imexam will plot the result.
        '''
        self.ax = ax
        self.fig = ax.get_figure()
        if plot_ax is not None:
            self.plot_ax = plot_ax
        self.connected = self.fig.canvas.mpl_connect('key_press_event', self.onKeyPress)

        log.info('Connecting imexam in slot %i' % self.connected)

    def disconnect(self):
        if self.connected is not None:
            log.info('Disconnecting imexam from slot %i' % self.connected)
            self.fig.canvas.mpl_disconnect(self.connected)
            self.ax = None
            self.fig = None
            self.plot_ax = None
            self.connected = None

    def _dummy_plot(self):
        self.plot_ax.plot(range(5), range(5))
        self.plot_ax.set_xlabel('dummy x')
        self.plot_ax.set_ylabel('dummy y')
        self.plot_ax.set_title('dummy title')
        self.plot_ax.set_xlim([0, 4])
        self.plot_ax.set_ylim([0, 4])
        self.draw()

    def set_labels(self, title, xlabel, ylabel):
        self.plot_ax.set_xlabel(xlabel)
        self.plot_ax.set_ylabel(ylabel)
        self.plot_ax.set_title(title)

    def draw(self):
        self.plot_ax.get_figure().tight_layout()
        self.plot_ax.get_figure().canvas.draw()

    def onKeyPress(self, event):
        log.debug('Key pressed: %s' % event.key)
        if event.inaxes != self.ax:
            return
        if self.hdu == None:
            return
        try:
            self.do_option(event.key, event.xdata, event.ydata)
        except Exception as e:
            log.error(e)

    def do_option(self, key, x, y):
        if key not in self.commands.keys():
            log.error('Key %s not valid.' % str(key))
        self.commands[key][0](x, y)

    def aper_phot(self, x, y):
        '''
        `a` key from imexam
        '''
        self.print_text.emit('x: %.1f, y:%.1f' % (x, y))

    def line_fit(self, x, y):
        '''
        `j` key from imexam
        '''
        self._check_plot_ax()
        self._dummy_plot()

    def column_fit(self, x, y):
        '''
        `k` key from imexam
        '''
        self._check_plot_ax()
        self._dummy_plot()

    def report_stat(self, x, y):
        '''
        `m` key from imexam
        '''
        self._check_plot_ax()
        self._dummy_plot()

    def plot_line(self, x, y):
        '''
        `l` key from imexam
        '''
        self._check_plot_ax()
        self._dummy_plot()

    def plot_column(self, x, y):
        '''
        `c` key from imexam
        '''
        self._check_plot_ax()
        self._dummy_plot()

    def radial_profile(self, x, y):
        '''
        `r` key from imexam
        '''
        self._check_plot_ax()
        self._dummy_plot()

    def histogram(self, x, y):
        '''
        `h` key from imexam
        '''
        self._check_plot_ax()
        data, xdata, ydata = trim_array(self.hdu.data, np.indices(self.hdu.data.shape),
                                        self.config.box_size, (x, y))
        self.plot_ax.hist(data.ravel(), self.config.hist_bins)
        self.set_labels('Histogram', 'Value', 'Number')
        self.draw()

    def contour(self, x, y):
        '''
        `e` key from imexam
        '''
        self._check_plot_ax()
        data, xdata, ydata = trim_array(self.hdu.data, np.indices(self.hdu.data.shape),
                                        self.config.box_size, (x, y))
        self.plot_ax.contour(xdata, ydata, data, self.config.n_contours)
        self.set_labels('Contour', 'Column', 'Line')
        self.draw()

    def surface(self, x, y):
        '''
        `s` key from imexam
        '''
        self._check_plot_ax()
        #data, xdata, ydata = trim_array(self.hdu.data, np.indices(self.hdu.data.shape),
        #                                self.config.box_size, (x, y))
        #self.plot_ax.plot_surface(xdata, ydata, data, rstride=1, cstride=1,
        #                          cmap=self.config.cmap, alpha=0.6)
        log.error('Surface plot not implemented yet.')

def imexamine(hdu):
    '''
    Created a Imexam instance, the proper axes in Matplotlib, and examines an
    HDUImage instance.
    '''
    def _create_plot_ax():
        fig = plt.Figure()
        return fig.add_subplot(111)

    imexam = Imexam()
    imexam.set_hdu(hdu)

    imexam.plot_ax = _create_plot_ax()
    imexam.print_text.connect(print)

    ax = _create_plot_ax()

    ax.imshow(hdu.data)

    imexam.connect(ax)

    plt.show()
