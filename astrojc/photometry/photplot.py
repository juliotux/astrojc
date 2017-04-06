import numpy as np

from ..time import year2jd

def plot_twiny(ax, begin, end, step):
    year_labels = np.arange(begin,end,step)
    year_ticks = year2jd(year_labels)
    axyear = ax.twiny()
    axyear.set_xticks(year_ticks)
    axyear.set_xticklabels(year_labels)
    axyear.set_xlim(ax.get_xlim())
    axyear.set_xlabel('Year')

def plot_lightcurve(ax, jd, flux, jd_error=None, flux_error=None, offset=0, fmt='ko',
                    label=None, pltkwargs={}):
    if flux_error is not None or jd_error is not None:
        ax.errorbar(jd, flux + offset, xerr=jd_error, yerr=flux_error, fmt=fmt, label=label
                    **pltkwargs)
    else:
        ax.plot(jd, flux + offset, fmt, label=label, **pltkwargs)
