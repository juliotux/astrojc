'''
New Python implementation of IRAF imexam functionalities.

Includes a new key set to fits explorer.
'''

from astropy.modeling import models, fitting
from matplotlib import pyplot as plt
import numpy as np

from scipy.spatial import cKDTree

from ..logging import log
from ..signal import MySignal
from ..math.array import trim_array, xy2r
from ..photometry.source_find import FindAdapter

class ImexamConfig():
    def __init__(self):
        self.box_size = 15   #Box size for contour and surface plots and fitting
        self.n_contours = 10 #Number of contours to be draw
        self.cmap = 'viridis'#Cmap to use in plots
        self.hist_bins = 12  #Number of bins of the histogram plot
        self.point_plot_fmt = 'wo'
        self.point_plot_size = 6
        self.line_plot_fmt = 'w-'
        self.line_plot_width = 2

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
                         's': (self.surface, 'Display a surface plot around the cursor location'),
                         'v': (self.vector, 'Displays a vector plot between 2 selected points')}

        self.ax = None
        self.fig = None
        self.hdu = None
        self.plot_ax= None
        self.connected = None
        self.press = None

        self.plotted = []

        self.print_text = MySignal()

        self.config = ImexamConfig()

        self.finder = FindAdapter(use_photutils=True, use_sep=False)

    def _check_plot_ax(self):
        if self.plot_ax is None:
            log.error('No plot axes connected.')
        else:
            self.plot_ax.cla()
        for i in self.plotted:
            i[0].remove()
            del i
        self.plotted = []

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

    def _extract_data(self, x, y):
        return trim_array(self.hdu.data, np.indices(self.hdu.data.shape),
                          self.config.box_size, (x, y))

    def set_labels(self, title, xlabel, ylabel):
        self.plot_ax.set_xlabel(xlabel)
        self.plot_ax.set_ylabel(ylabel)
        self.plot_ax.set_title(title)

    def draw(self):
        self.plot_ax.get_figure().tight_layout()
        self.plot_ax.get_figure().canvas.draw()
        self.ax.get_figure().canvas.draw()

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
        log.error('Aperture photometry not implemented yet.')

    def line_fit(self, x, y):
        '''
        `j` key from imexam
        '''
        self._check_plot_ax()
        log.error('Line fit not implemented yet.')

    def column_fit(self, x, y):
        '''
        `k` key from imexam
        '''
        self._check_plot_ax()
        log.error('Column fit not implemented yet.')

    def report_stat(self, x, y):
        '''
        `m` key from imexam
        '''
        self._check_plot_ax()
        log.error('Statistics report not implemented yet.')

    def plot_line(self, x, y):
        '''
        `l` key from imexam
        '''
        self._check_plot_ax()
        lin = np.int(y)
        data = self.hdu.data[lin,:]

        self.plot_ax.plot(data)
        self.set_labels('Line Plot', 'Column', 'Value')
        self.draw()

    def plot_column(self, x, y):
        '''
        `c` key from imexam
        '''
        self._check_plot_ax()
        col = np.int(x)
        data = self.hdu.data[:,col]

        self.plot_ax.plot(data)
        self.set_labels('Column Plot', 'Line', 'Value')
        self.draw()

    def radial_profile(self, x, y):
        '''
        `r` key from imexam
        '''
        log.error('Radial profile plot not implemented yet.')

    def histogram(self, x, y):
        '''
        `h` key from imexam
        '''
        self._check_plot_ax()
        data, xdata, ydata = self._extract_data(x, y)
        self.plot_ax.hist(data.ravel(), self.config.hist_bins)
        self.set_labels('Histogram', 'Value', 'Number')
        self.draw()

    def contour(self, x, y):
        '''
        `e` key from imexam
        '''
        self._check_plot_ax()
        data, xdata, ydata = self._extract_data(x, y)
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

    def vector(self, x, y):
        self._check_plot_ax()
        if self.press is None:
            self.press = ('v', x, y)
            self.plotted.append(self.ax.plot(x, y, self.config.point_plot_fmt,
                                             label='imexam_vector_point',
                                             markersize = self.config.point_plot_size))
            self.plot_ax.text(0, 0, 'Press \'v\' again.', ha='center', va='center',
                              fontsize=18)
            self.plot_ax.set_xlim(-1, 1)
            self.plot_ax.set_ylim(-1, 1)
            self.draw()
        elif self.press[0] == 'v':
            key, x0, y0 = self.press
            self.press = None

            n = int(np.hypot(x-x0, y-y0))
            xdata, ydata = np.linspace(x0, x, n), np.linspace(y0, y, n)
            data = self.hdu.data[ydata.astype('int16'), xdata.astype('int16')]
            dist, data = xy2r(xdata, ydata, data, x0, y0)

            self.plot_ax.plot(dist, data, 'k-')
            self.plotted.append(self.ax.annotate("", xy=(x0, y0), xytext=(x1, y1),
                                                 arrowprops={arrowstyle="<-",
                                                             connectionstyle="arc3",
                                                             color=self.config.line_plot_fmt[0],
                                                             linewidth=self.config.line_plot_width})
            self.plotted.append(self.ax.plot([x0, x], [y0, y], self.config.point_plot_fmt,
                                             label='imexam_vector_point',
                                             markersize = self.config.point_plot_size))
            self.set_labels('Vector Plot', 'Distance', 'Value')
            self.draw()
        else:
            log.error('You changed the pressing key. Clearing!')

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
