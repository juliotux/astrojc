'''
Reads data from some specified formats, distributed by some observatories.
'''
from astropy.io import ascii
from astropy.table import Table
from astropy.io.votable import parse as vo
import astropy.units as u
import astropy.constants as c
import numpy as np

from .band_converter import uvby2UBV
from ..logging import log

__all__ = ['read_wise', 'read_allwise', 'read_asas', 'read_aavso', 'read_vg',
           'read_dasch', 'read_ltpv', 'read_generic_table', 'LTPV_data']

#TODO: better documentate these functions.

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
    raise NotImplementedError

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

class LTPV_data():
    def __init__(self, data_dir, star_name):
        self.data_type = np.dtype([('HJD',np.float),
                                   ('u',np.float),
                                   ('b',np.float),
                                   ('v',np.float),
                                   ('y',np.float),
                                   ('m1',np.float),
                                   ('c1',np.float),
                                   ('U',np.float),
                                   ('B',np.float),
                                   ('V',np.float)])
        self.star_name = star_name
        self.data_dir = data_dir
        if(self.data_dir[-1] != '/'):
            self.data_dir = self.data_dir + '/'
        self.star_id = self.get_id()
        self.data = self.load_data()

    def get_id(self):
        star = ''
        star_name = self.star_name.replace(' ','')
        for i in ['1','2','3','4']:
            mycat = ascii.read(self.data_dir+'stars'+i,readme=self.data_dir+'ReadMe'+i)
            for j in mycat:
                if( str(j['HD']).replace(' ','') == star_name or str(j['HR']).replace(' ','') == star_name or
                    str(j['DM']).replace(' ','') == star_name or str(j['Name2']).replace(' ','') == star_name ):
                    star = j['Name']
        if star == '':
            print('Star '+self.star_name+' not found. Check the name.')
            return ''
        else:
            return star

    def load_data(self):
        tempdata = []
        for i in ['1','2','3','4']:
            mycat = ascii.read(self.data_dir+'ltpv'+i,
                               readme=self.data_dir+'ReadMe'+i)
            for j in mycat[mycat['Name'] == self.star_id]:
                if i == '1' or i == '2' or i == '3':
                    hjd = j['HJD']
                    u = j['umag']
                    b = j['bmag']
                    v = j['vmag']
                    y = j['ymag']
                    m1 = (v-b)-(b-y)
                    c1 = (u-v)-(u-b)
                    (U,B,V) = uvby2UBV(u, b, v, y)
                    tempdata.append((hjd,u,b,v,y,m1,c1,U,B,V))
                elif i == '4':
                    hjd = j['HJD']
                    y = j['ymag']
                    b = y + j['b-y']
                    m1 = j['m1']
                    c1 = j['c1']
                    v = m1 + 2*b - y
                    u = c1 + 2*v - b
                    (U,B,V) = uvby2UBV(u, b, v, y)
                    tempdata.append((hjd,u,b,v,y,m1,c1,U,B,V))
        return np.array(tempdata, dtype=self.data_type)

def read_asas(fname, filter_grade=('A', 'B'), filter_mag=29.999):
    '''
    Read All Sky Automated Survey (ASAS) data, filtering the `GRADE` (data
    quality) and invalid magnitude.
    '''
    dtype = np.dtype([('HJD',np.float),('MAG_0',np.float),('MAG_1',np.float),
                      ('MAG_2',np.float),('MAG_3',np.float),('MAG_4',np.float),
                      ('MER_0',np.float),('MER_1',np.float),('MER_2',np.float),
                      ('MER_3',np.float),('MER_4',np.float),('GRADE','S8'),
                      ('FRAME',np.int)])
    asas = ascii.read(fname, header_start=53)
    filt = np.zeros(len(asas), dtype=bool)
    for i in filter_grade:
        filt = filt | (asas['GRADE'] == i)
    for i in range(5):
        filt = filt & (asas['MAG_%i' % i] < filter_mag)

    log.debug('Size of table before filtering: %i' % len(asas))
    asas = asas[filt]
    log.debug('Size of table after filtering: %i' % len(asas))

    fmag = np.nanmean(np.array([asas['MAG_0'],asas['MAG_1'],asas['MAG_2'],
                                asas['MAG_3'],asas['MAG_4']]), axis=0)
    stdmag = np.nanstd(np.array([asas['MAG_0'],asas['MAG_1'],asas['MAG_2'],
                                 asas['MAG_3'],asas['MAG_4']]), axis=0)

    return np.array(list(zip(asas['HJD']+2450000, fmag, stdmag)),
                    dtype=np.dtype([('JD', 'f8'), ('MAG', 'f8'),
                                    ('MAG_ERR', 'f8')]))

def read_aavso(fname, band, validation_flags=None):
    '''
    Reads one band AAVSO data from a file. Validation flags can be filtered.

    For Vis. filter, it's recomended to use `validation_flags` = ('V', 'Z').
    '''
    t = ascii.read(fname)
    filt = t['Band'] == band
    if validation_flags is not None:
        filt1 = np.zeros(len(t), dtype=bool)
        vf = np.array(validation_flags)
        for i in vf:
            filt1 = filt1 | t['Validation Flag'] == i
        filt = filt & filt1

    log.debug('Size of table before filtering: %i' % len(t))
    t = t[filt]
    log.debug('Size of table after filtering: %i' % len(t))
    t['MAG'] = t['Magnitude']
    t['MAG_ERR'] = t['Uncertainty']

    return np.array(t['JD', 'MAG', 'MAG_ERR'])

def read_dasch(fname, max_error):
    '''
    Read DASCH VO table '.xml.gz'.
    '''
    dasch = vo(fname).resources[0].tables[0]
    filt = dasch['magcal_local_rms'] <= max_error

    log.debug('Size of table before filtering: %i' % len(dasch))
    dasch = dasch[filt]
    log.debug('Size of table after filtering: %i' % len(dasch))

    dasch['JD'] = dasch['ExposureDate']
    dasch['MAG'] = dasch['magcal_magdep']
    dasch['MAG_ERR'] = dasch['magcal_local_rms']

    return np.array(dasch['JD', 'MAG', 'MAG_ERR'])

def read_vg(fname, band):
    '''
    Read data from Van Genderen works, converted to a ascii table with the
    following columns:
    JD-2440000, V, V-B, B-U, U-W, B-L, Vj, B-Vj.

    You may specify the wanted filter.
    '''
    vg = ascii.read(fname)
    vg['JD'] = vg['JD-2440000'] + 2440000
    if band in ['V', 'Vj']:
        vg['MAG'] = vg[band]
    elif band == 'B':
        vg['MAG'] = vg['V'] - vg['V-B']
    elif band == 'U':
        vg['MAG'] = vg['V'] - vg['V-B'] - vg['B-U']
    elif band == 'W':
        vg['MAG'] = vg['V'] - vg['V-B'] - vg['B-U'] - vg['U-W']
    elif band == 'L':
        vg['MAG'] = vg['V'] - vg['V-B'] - vg['B-L']
    elif band == 'Bj':
        vg['MAG'] = vg['Vj'] + vg['B-Vj']
    vg['MAG_ERR'] = [0]*len(vg)

    return np.array(vg['JD', 'MAG', 'MAG_ERR'])

def read_ltpv(ltpv_location, star_hd, band):
    '''
    Read LTPV data.

    filter : ('u', 'v', 'b', 'y', 'U', 'V', 'B', 'm1', 'c1')
    '''
    if band not in ('u', 'v', 'b', 'y', 'U', 'V', 'B', 'm1', 'c1'):
        raise ValueError('`band` must be: (u, v, b, y, U, V, B, m1, c1)')
    ltpv = LTPV_data(ltpv_location, star_hd).data
    return np.array(list(zip(ltpv['HJD'], ltpv[band], [0]*len(ltpv))),
                    dtype=np.dtype([('JD', 'f8'), ('MAG', 'f8'),
                                    ('MAG_ERR', 'f8')]))

def read_generic_table(table, jd_key, mag_key, mag_err_key=None, jd_correct=0):
    '''
    Read a generic astropy table containing the data and correct the Julian Date by a
    value if it's needed (like the data in MJD).
    '''
    t = Table()
    if not isinstance(table, Table):
        try:
            ascii.read(table)
        except:
            raise ValueError('You must pass a astropy table instance or a file'
            'name that can be parsed by astropy.io.ascii.read')
    t['JD'] = table[jd_key] + jd_correct
    t['MAG'] = table[mag_key]
    if mag_err_key is not None:
        t['MAG_ERR'] = table[mag_err_key]
    else:
        t['MAG_ERR'] = [0]*len(table)

    return np.array(t)
