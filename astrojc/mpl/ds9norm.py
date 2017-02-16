'''
Module to mimic the normalization of images like DS9, with interacting with
matplotlib.

Similar to glue-viz ds9norm package, but with interacting.
'''

import numpy as np
from matplotlib.colors import Normalize

def get_clip(data, percentile_lo, percentile_hi):
    '''
    Get the values of a low and high percentiles of an array.
    '''
    d = np.asarray(data)
    if ~np.isfinite(d).any():
        return (0.0, 1.0)

    d = d[np.isfinite(data)]
    lo, hi = np.percentile(d, [percentile_lo, percentile_hi])
    return lo, hi

def norm(x, lo, hi):
    '''
    Clip the data with lo and hi limits.
    '''
    result = np.clip(x, lo, hi)
    return result

def cscale(x, bias, contrast):
    '''
    Apply bias and contrast scaling. Overwrite input.

    Parameters
    ----------
        x : array
            Values between 0 and 1
        bias : float
        contrast : float

    Returns
    -------
    The input x, scaled inplace
    '''
    x = np.subtract(x, bias, out=x)
    x = np.multiply(x, contrast, out=x)
    x = np.add(x, 0.5, out=x)
    x = np.clip(x, 0, 1, out=x)
    return x

def linear_warp(x, vmin, vmax, bias, contrast):
    return cscale(clip(x, vmin, vmax), bias, contrast)

def log_warp(x, vmin, vmax, bias, contrast, exp=1000.0):
    black = x < vmin
    x = norm(x, vmin, vmax)
    x = np.multiply(exp, x, out=x)
    # sidestep numpy bug that masks log(1)
    # when out is provided
    x = np.add(x, 1.001, out=x)
    x = np.log(x, out=x)
    x = np.divide(x, np.log(exp + 1.0), out=x)
    x = cscale(x, bias, contrast)
    return x

def pow_warp(x, vmin, vmax, bias, contrast, exp=1000.0):
    x = norm(x, vmin, vmax)
    x = np.power(exp, x, out=x)
    x = np.subtract(x, 1, out=x)
    x = np.divide(x, exp - 1)
    x = cscale(x, bias, contrast)
    return x

def sqrt_warp(x, vmin, vmax, bias, contrast):
    x = norm(x, vmin, vmax)
    x = np.sqrt(x, out=x)
    x = cscale(x, bias, contrast)
    return x

def squared_warp(x, vmin, vmax, bias, contrast):
    x = norm(x, vmin, vmax)
    x = np.power(x, 2, out=x)
    x = cscale(x, bias, contrast)
    return x

def asinh_warp(x, vmin, vmax, bias, contrast):
    x = norm(x, vmin, vmax)
    x = np.divide(np.arcsinh(np.multiply(x, 10, out=x), out=x), 3, out=x)
    x = cscale(x, bias, contrast)
    return x

warpers = dict(linear=linear_warp,
               log=log_warp,
               sqrt=sqrt_warp,
               power=pow_warp,
               squared=squared_warp,
               arcsinh=asinh_warp)

class DS9Norm(Normalize, object):
    def __init__(self, stretch='linear',
                 bias=0.5, contrast=1.0,
                 clip_lo=0., clip_hi=100.):
        self.ax = None
        self.colorbar = None
        self._stretch = stretch
        self.bias = bias
        self.contrast = contrast
        self._clip_lo = clip_lo
        self._clip_hi = clip_hi

        Normalize.__init__(self, None, None, None)

        self.press = None

    @property
    def clip_lo(self):
        return self._clip_lo

    @clip_lo.setter
    def clip_lo(self, value):
        self._clip_lo = value
        self.vmin = self.vmax = None

    @property
    def clip_hi(self):
        return self._clip_hi

    @clip_hi.setter
    def clip_hi(self, value):
        self._clip_hi = value
        self.vmin = self.vmax = None

    @property
    def stretch(self):
        return self._stretch

    @stretch.setter
    def stretch(self, value):
        if value not in warpers:
            raise ValueError("Invalid stretch: %s\n Valid options are: %s" %
                             (value, warpers.keys()))
        self._stretch = value

    def get_bias_contrast(self, x_pos, y_pos):
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()

        #bias is defined by the x axis
        bias = float(x_pos - xlim[0])/(xlim[1]-xlim[0])
        bias = np.abs(bias*(self._clip_hi - self._clip_lo)*self._clip_hi)

        #contrast is defined by the y axis
        contrast = float(y_pos - ylim[0])/(ylim[1]-ylim[0])

        return bias, contrast

    def ds9norm_factory(self, ax, colorbar=None):
        self.ax = ax
        self.colorbar = colorbar
        self.fig = self.ax.get_figure()

        self.fig.canvas.mpl_connect('button_press_event', self.onPress)
        self.fig.canvas.mpl_connect('button_release_event', self.onRelease)
        self.fig.canvas.mpl_connect('motion_notify_event', self.onMotion)

    def onPress(self, event):
        if event.inaxes != self.ax or event.button != 3:
            return
        self.press = event.xdata, event.ydata
        self.bias, self.contrast = self.get_bias_contrast(event.xdata, event.ydata)

    def onMotion(self, event):
        if press == None:
            return
        if event.inaxes != self.ax:
            return
        self.bias, self.contrast = self.get_bias_contrast(event.xdata, event.ydata)

    def onRelease(self, event):
        self.press = None
        self.ax.figure.canvas.draw()

    def autoscale(self, A):
        self.update_clip(A)

    def autoscale_None(self, A):
        if self.vmin is None or self.vmax is None:
            self.update_clip(A)

    def update_clip(self, image):
        self.vmin, self.vmax = get_clip(image, self._clip_lo, self._clip_hi)

    def __call__(self, value, clip=False):
        self.autoscale_None(value)

        inverted = self.vmax <= self.vmin

        hi, lo = max(self.vmin, self.vmax), min(self.vmin, self.vmax)

        warp = warpers[self.stretch]
        result = warp(value, lo, hi, self.bias, self.contrast)

        if inverted:
            result = np.subtract(1, result, out=result)

        result = np.ma.MaskedArray(result, copy=False)

        return value
