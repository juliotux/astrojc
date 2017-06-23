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
        ax.errorbar(jd, flux + offset, xerr=jd_error, yerr=flux_error, fmt=fmt, label=label,
                    **pltkwargs)
    else:
        ax.plot(jd, flux + offset, fmt, label=label, **pltkwargs)

def resample_lightcurve(jd, flux, n_box, method = 'median', errors = True):
    '''Resample a light curve in n_box points in interval.'''
    if method == 'median':
        func = np.nanmedian
    else:
        func = np.nanmean

    lo = np.nanmin(jd)
    hi = np.nanmax(jd)
    interval = (hi - lo)/n_box

    n_jd = np.array([np.nan]*n_box)
    n_flux = np.array([np.nan]*n_box)
    n_flux_error = np.array([np.nan]*n_box)

    for i in range(n_box):
        lim = (lo+(i*interval), lo+((i+1)*interval))
        filt = np.where((jd >= lim[0]) & (jd < lim[1]))
        if len(filt) > 0:
            n_jd[i] = func(jd[filt])
            n_flux[i] = func(flux[filt])
            n_flux_error[i] = np.nanstd(flux[filt])

    if errors:
        return n_jd, n_flux, n_flux_error
    else:
        return n_jd, n_flux
