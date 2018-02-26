from astropy.time import Time
import numpy as np

from .math.opd_utils import opd2jd

__all__ = ['year2jd', 'to_datetime', 'opd2jd']


def to_datetime(jd_times, sbfmt=None):
    time = [None]*len(jd_times)
    for i in range(len(jd_times)):
        j = Time(jd_times[i], format='jd')
        if sbfmt is not None:
            j.out_subfmt = sbfmt
        time[i] = np.datetime64(j.iso)
    return time


def year2jd(year):
    '''
    Convert a year to a Julian Date.
    '''
    return Time(year, format='decimalyear').jd
