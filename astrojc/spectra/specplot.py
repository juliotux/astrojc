'''
Generic module to plot spectra in matplotlib.
'''

import numpy as np
from .spectrum import Spectrum

__all__ = ['to_vel', 'plot_spec', 'plot_spec_vel', 'ticks_ids', 'plot_line_id']

def to_vel(wcenter, wave):
    '''
    Returns the dispersion in units of velocity to a reference wavelength.
    '''
    return ((wave - wcenter)/wcenter)*3E5

def _join_strings(arr1, arr2):
    '''
    Joins two special strings of multiplets
    '''
    outarray = []
    for i in xrange(len(arr1)):
        if arr2[i] == '-':
            outarray.append(str(arr1[i]))
        else:
            outarray.append(str(arr1[i]) + ' m' + str(arr2[i]))
    return outarray

def plot_spec(ax, spec, offset=0, fmt='k-', label=None, pltkwargs={}, **kwargs):
    '''
    Plots just one spectrum with a offset in and axis.
    '''
    ax.plot(spec.wave, spec.flux + offset, fmt, label=label, **pltkwargs,
            **kwargs)

def plot_spec_vel(ax, spec, wcenter, fmt='k-', offset=0, label=None,
                  pltkwargs={}, **kwargs):
    '''
    Plots just one spectrum with a offset in and axis.
    '''
    ax.plot(to_vel(wcenter, spec.wave), spec.flux + offset, fmt, label=label,
            **pltkwargs, **kwargs)

def ticks_ids(ax, wave, text):
    axn = ax.twiny()
    axn.set_xticks(wave)
    axn.set_xticklabels(text, rotation=90)
    axn.set_xlim(ax.get_xlim())
    axn.grid()

def plot_line_id(ax, wave, y, text, box_size, spec=None, flux=None,
                 wavetext=None, pltkwargs={}):
    if wavetext is None:
        wavetext = wave #colocar iterador para definir novo wavetext
    if spec is not None:
        ax.annotate(text, xy=(wave, np.interp(spec.wave, spec.flux, wave)),
                    xytext=(wavetext, y), arrowprops=dict(arrowstyle='-'),
                    ha='center', va='bottom', **pltkwargs)
    elif flux is not None:
        ax.annotate(text, xy=(wave, flux), xytext=(wavetext, y),
                    arrowprops=dict(arrowstyle='-'), ha='center', va='bottom',
                    **pltkwargs)
