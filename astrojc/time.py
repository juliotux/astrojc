from astropy.time import Time
import numpy as np

__all__ = ['year2jd', 'opd2jd', 'to_datetime']

def to_datetime(jd_times, sbfmt=None):
    time = [None]*len(jd_times)
    for i in xrange(len(jd_times)):
        j = Time(jd_times[i], format='jd')
        if sbfmt is not None:
            j.out_subfmt = sbfmt
        time[i] = np.datetime64(j.iso)
    return time

dict_meses = {v: k+1 for k,v in enumerate(['jan','fev','mar','abr','mai','jun',
                                           'jul','ago','set','out','nov','dez'])}

def opd2jd(opd_dates):
    '''
    Convert the OPD std date designation '99ago01' to Julian Date.
    '''
    res_arr = []
    opd_dates = np.array(opd_dates)
    for i in opd_dates:
        yy = int(i[0:2])
        mmm = dict_meses[i[2:5]]
        dd = int(i[5:])

        #Workaround for only 2 digits year
        if yy < 50:
            yy = 2000 + yy
        else:
            yy = 1900 + yy

        res_arr.append("%i-%i-%i" % (yy, mmm, dd))
    return Time(res_arr).jd

def year2jd(year):
    '''
    Convert a year to a Julian Date.
    '''
    return Time(year, format='decimalyear').jd
