from astropy.table import Table
from astroquery.vizier import Vizier
from ftplib import FTP
import gzip
import io
import numpy as np

hip_file = None


class ftpstore():
    def __init__(self):
        self.data = []

    def string(self):
        return b"".join(self.data)

    def store(self, data):
        self.data.append(data)


def get_tyc_file(tyc):
    f = FTP('cdsarc.u-strasbg.fr')
    f.login()
    f.cwd('pub/cats/I/239/epophot')
    fil = ftpstore()
    f.cwd("tep"+tyc[0][:2])
    f.retrbinary('RETR {}.gz'.format(tyc[0]), fil.store)
    f.close()
    return fil.string()


def downlad_tyc_hip(star_name):
    """Auto download tycho and hipparcos epoch photometry when available."""
    # get tycho and hipparcos names
    global hip_file

    tyc_query = Vizier.query_object(star_name, catalog='I/239/tyc_main')
    if len(tyc_query) > 0:
        tyc = tyc_query[0][0]['TYC'].split()
    else:
        tyc = None

    hip_query = Vizier.query_object(star_name, catalog='I/239/hip_main')
    if len(hip_query) > 0:
        hip = hip_query[0][0]['HIP']
    else:
        hip = None

    t = None
    if tyc is not None:
        phot = gzip.decompress(get_tyc_file(tyc))
        phot = phot.decode("utf-8").split('\n')
        process = False
        head = None
        transits = []
        for line in phot:
            if not process:
                line = [i.strip() for i in line.split('|')]
                if line[0] == tyc[0] and line[1] == tyc[1] and line[2] == tyc[2]:
                    head = line
                    process = True
            else:
                if line.strip() == '':
                    break
                line = [i.strip() for i in line.split('|')]
                line = [float(i) if i != '' else 0 for i in line]
                transits.append(tuple(line))
        if head is not None:
            a = np.array(transits,
                         dtype=np.dtype([('JD-2440000', 'f8'),
                                         ('BTmag', 'f4'),
                                         ('e_BTmag', 'f4'),
                                         ('BTbg', 'f4'),
                                         ('VTmag', 'f4'),
                                         ('e_VTmag', 'f4'),
                                         ('VTbg', 'f4'),
                                         ('z', 'f4'),
                                         ('PAz', 'f4'),
                                         ('Du', 'i4'),
                                         ('e_Di', 'i4'),
                                         ('F2', 'f4'),
                                         ('Tflg', 'i4')]))
            t = Table(a)
            for i,j in zip(head, ['TYC1', 'TYC2', 'TYC3', 'Ntr',
                                  'Nb', 'BTmag', 'e_BTmag',
                                  'BTmax', 'BTmin', 'Nv',
                                  'VTmag', 'e_VTmag', 'VTmax',
                                  'VTmin', 'Flag1', 'Flag2']):
                try:
                    fi = float(i)
                    ii = int(i)
                    if np.abs(fi - ii) < 0.001:
                        i = ii
                    else:
                        i = fi
                except:
                    pass
                t.meta[j] = i

    h = None
    if hip is not None:
        if hip_file is None:
            f = FTP('cdsarc.u-strasbg.fr')
            f.login()
            f.cwd('pub/cats/I/239/epophot')
            fil = ftpstore()
            f.retrbinary('RETR hep.gz', fil.store)
            f.close()
            hip_file = gzip.decompress(fil.string()).decode('utf-8').split('\n')
            del fil

        hip_str = str(hip).rjust(6)
        head = None
        transits = []
        process = False
        for line in hip_file:
            if not process and line[:6] == hip_str:
                head = line.split('|')
                process = True
            elif line == '' and process:
                process = False
                break
            elif process:
                line = line.split('|')
                line = [float(i) if i.strip() != '' else 0 for i in line]
                transits.append(tuple(line))

        if head is not None:
            a = np.array(transits, dtype=np.dtype([('JD-2440000', 'f8'),
                                                   ('HPmag', 'f4'),
                                                   ('e_HPmag', 'f4'),
                                                   ('Tflg', 'i4')]))
            h = Table(a)
            for i,j in zip(head, ['HIP', 'm_HIP', 'V-I', 'Ntr', 'Nh',
                                  'HPmag', 'e_HPmag', 'HPmax', 'HPmin',
                                  'Per', 'T0', 'VType', 'Flg1', 'Flg2']):
                try:
                    fi = float(i)
                    ii = int(i)
                    if np.abs(fi - ii) < 0.001:
                        i = ii
                    else:
                        i = fi
                except:
                    pass
                h.meta[j] = i

    return t, h
