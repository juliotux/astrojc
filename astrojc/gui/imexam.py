'''
New Python implementation of IRAF imexam functionalities.

Includes a new key set to fits explorer.
'''

from astropy.modeling import models, fitting
from matplotlib import pyplot as plt

from ..logging import log
from .my_signal import MySignal

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
                         '2': (self.new_plot_axes, 'Make the next plot in a new axes'),}

        self.ax = None
        self.fig = None
        self.hdu = None
        self.plot_ax= None

        self.create_plot_axes = MySignal()
        self.print_text = MySignal()

    def _check_plot_ax(self):
        if not isinstance(self.plot_ax, plt.Axes):
            self.plot_ax = self.create_plot_axes.emit()

    def set_hdu(self, hdu):
        self.hdu = hdu

    def imexam_factory(self, ax, plot_ax=None):
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
        self.plot_ax = plot_ax
        self.connected = self.fig.canvas.mpl_connect('key_press_event',
                                                     self.onKeyPress)

    def onKeyPress(self, event):
        if event.inaxes != self.ax:
            return
        if self.hdu == None:
            return
        self.do_option(event.key, event.xdata, event.ydata)

    def do_option(self, key, x, y):
        if key not in self.commands.keys():
            log.error('Key %s not valid.' % str(key))
        self.commands[key][0](x, y)

    def aper_phot(self, x, y):
        '''
        `a` key from imexam
        '''
        self.print_text.emit('Not implemented yet')

    def line_fit(self, x, y):
        '''
        `j` key from imexam
        '''
        self._check_plot_ax()

    def column_fit(self, x, y):
        '''
        `k` key from imexam
        '''
        self._check_plot_ax()

    def report_stat(self, x, y):
        '''
        `m` key from imexam
        '''
        self._check_plot_ax()

    def plot_line(self, x, y):
        '''
        `l` key from imexam
        '''
        self._check_plot_ax()

    def plot_column(self, x, y):
        '''
        `c` key from imexam
        '''

    def radial_profile(self, x, y):
        '''
        `r` key from imexam
        '''
        self._check_plot_ax()

    def histogram(self, x, y):
        '''
        `h` key from imexam
        '''
        self._check_plot_ax()

    def contour(self, x, y):
        '''
        `e` key from imexam
        '''
        self._check_plot_ax()

    def surface(self, x, y):
        '''
        `s` key from imexam
        '''
        self._check_plot_ax()

    def new_plot_axes(self, x, y):
        '''
        `2` key from imexam
        '''
        self.plot_ax = self.new_ax.emit()

def imexamine(hdu):
    '''
    Created a Imexam instance, the proper axes in Matplotlib, and examines an
    HDUImage instance.
    '''
    def _create_plot_ax():
        fig = plt.gcf()
        return fig.gca()

    imexam = Imexam()
    imexam.set_hdu(hdu)

    imexam.create_plot_axes.connect(_create_plot_ax)
    imexam.print_text.connect(print)

    ax = _create_plot_ax()

    ax.imshow(hdu.data)

    imexam.imexam_factory(ax)

    plt.show()
