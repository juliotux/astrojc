'''
Reads data from some specified formats, distributed by some observatories.
'''
from astropy.io import ascii
from astropy.table import Table
import astropy.units as u
import astropy.constants as c
import numpy as np

__all__ = ['read_wise', 'read_allwise']

def _read_wise_band(table, band_number, mjd, mjd_margin, deviation_limit):
    filt1 = (t['mjd'] >= mjd-mjd_margin) & (t['mjd'] <= mjd+mjd_limit)
    w = table['w%impro' % band_number][filt]
    werr = table['w%isigmpro'% band_number]['filt']
    wmed = np.median(w)
    filt2 = abs(w - wmed) <= deviation_limit
    return w[filt], werr[filt]

def read_wise(fname, mjd, mjd_limit=10, deviation_limit=0.2):
    '''
    Read wise data from a file.
    '''
    #TODO: follow the read_allwise example and reimplement it using _read_wise_band


def read_allwise(table, fname, mjd, mjd_limit=10, deviation_limit=0.2):
    "Wise são plotados com 'wo'"
    t = ascii.read(fname, format='ipac')
    filt1 = t['mjd'] >= mjd-mjd_limit
    filt2 = t['mjd'] <= mjd+mjd_limit
    filt = filt1 & filt2
    w1 = t['w1mpro_ep'][filt]
    w2 = t['w2mpro_ep'][filt]
    w3 = t['w3mpro_ep'][filt]
    w4 = t['w4mpro_ep'][filt]
    w1_err = t['w1sigmpro_ep'][filt]
    w2_err = t['w2sigmpro_ep'][filt]
    w3_err = t['w3sigmpro_ep'][filt]
    w4_err = t['w4sigmpro_ep'][filt]
    med_w1 = np.nanmedian(w1)
    med_w2 = np.nanmedian(w2)
    med_w3 = np.nanmedian(w3)
    med_w4 = np.nanmedian(w3)
    for W1, W2, Sig1, Sig2 in zip(w1, w2, w3, w4, w1_err, w2_err, w3_err, w4_err):
        if W1 - med_w1 < deviation_limit:
            table.add_row(['W1', W1, Sig1, 'wo'])
        if W2 - med_w2 < deviation_limit:
            table.add_row(['W2', W2, Sig2, 'wo'])
        if W3 - med_w3 < deviation_limit:
            table.add_row(['W3', W3, Sig3, 'wo'])
        if W3 - med_w3 < deviation_limit:
            table.add_row(['W4', W4, Sig4, 'wo'])
