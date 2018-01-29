"""Read PCCDPACK result files. Generate scripts to execute in iraf for
pccdpack."""

import os
import re
from collections import OrderedDict
import numpy as np

from astropy.io import ascii as asci
from astropy.io import fits
from tempfile import NamedTemporaryFile


def read_log(log_file):
    """Read the .log file from pccdpack and return a dict with the results"""
    f = open(log_file, 'r')
    f = f.readlines()

    result = OrderedDict()
    # {star_n:{aperture:{'params':{},z:[]},'positions':[]}
    star = None

    wave_n_pos = 16
    wave_pos = None

    read_star = False
    read_apertures = False
    _aperture = 0

    params_index = None
    z_index = None

    i = 0
    while i < len(f):
        lin = f[i].strip('\n').strip(' ')
        if re.match('No. of waveplate positions :\s+\d+', lin):
            wave_n_pos = int(re.findall('\d+', lin)[0])
        if re.match('Waveplate pos. : [\s+\d+]+ = \d+', lin):
            wave_pos = [int(v) for v in re.findall('\d+', lin)[:-1]]
        if re.match('No. of apertures observed:\s+\d+', lin):
            n_apertures = int(re.findall('\d+', lin)[0])

        if re.match('STAR # \s+\d+\s+ .*', lin):
            star = int(re.findall('\d+', lin)[0])
            result[star] = OrderedDict()
            result[star]['positions'] = wave_pos
            read_star = True

        if re.match('APERTURE = \s+\d+.\d+.', lin):
            aperture = float(re.findall('\d+.\d+', lin)[0])
            result[star][aperture] = {'z': [], 'params': OrderedDict()}
            _aperture += 1
            read_apertures = True
            params_index = i+2
            z_index = i+5

        if read_star and read_apertures:
            plin = f[params_index].strip('\n').strip(' ')
            params = [float(k) for k in plin.split()]
            plin = f[params_index-1].strip('\n').strip(' ')
            pnames = plin.split()
            del plin
            for n, v in zip(pnames, params):
                result[star][aperture]['params'][n] = v

            z = []
            while len(z) < wave_n_pos:
                zlin = f[z_index].strip('\n').strip(' ')
                for k in zlin.split():
                    z.append(float(k))
                z_index += 1
                del zlin
            result[star][aperture]['z'] = z

            if _aperture == n_apertures:
                read_star = False
            read_apertures = False
            i = z_index
        i += 1
    return result


def read_out(file_out, file_ord=None):
    """Read the out file, with optionally the ord file for x,y coords."""
    fout = asci.read(file_out)
    if file_ord is not None:
        ford = asci.read(file_ord)
        fout['X0'] = [ford['XCENTER'][2*i] for i in range(len(fout))]
        fout['Y0'] = [ford['YCENTER'][2*i] for i in range(len(fout))]
        fout['X1'] = [ford['XCENTER'][2*i+1] for i in range(len(fout))]
        fout['Y1'] = [ford['YCENTER'][2*i+1] for i in range(len(fout))]
    return fout


def create_script(result_dir, image_list, star_name,
                  filter, plot=True, apertures=np.arange(1, 21, 1),
                  r_ann=60.0, r_dann=10.0, n_retarder_positions=16,
                  gain_key='GAIN', readnoise_key='RDNOISE',
                  auto_pol=False, normalize=True, retarder='half'):
    """Creates a script to easy execute pccdpack for a set of images."""
    wd = result_dir
    try:
        os.makedirs(wd)
    except Exception as e:
        print(e)

    # FIXME: for a bug in pccdpack, the images have to be in the data_dir
    # FIXME: when not, the pccdpack cannot read the generated .mag.1 files
    imlist = [os.path.basename(i) for i in image_list]
    d = os.path.dirname(image_list[0])
    f = NamedTemporaryFile(prefix='pccdpack', suffix='.txt')
    with open(f.name, 'w') as tmp:
        [tmp.write(i + '\n') for i in imlist]

    kwargs = {'object': star_name,
              'lista': '@' + f.name,
              'imagem': imlist[0],
              'pospars': np.min([len(image_list), n_retarder_positions]),
              'nume_la': len(image_list),
              'nab': len(apertures),
              'abertur': ','.join([str(i) for i in apertures]),
              'anel': r_ann,
              'danel': r_dann,
              'autoabe': 'no',
              'desloca': 'no',
              'confirm': 'no'}

    ref = fits.open(image_list[0])[0]
    kwargs['gan'] = float(ref.header[gain_key])
    kwargs['readnoi'] = float(ref.header[readnoise_key])
    kwargs['tamanho'] = int(np.max(ref.data.shape))
    if retarder == 'quarter':
        kwargs['circula'] == 'yes'

    if auto_pol:
        comm = 'auto_pol'
    else:
        comm = 'padrao_pol'

    command = comm + ' '
    command += ' '.join(['{}={}'.format(i, j) for (i, j) in kwargs.items()])

    print("cd {}".format(d))
    print(command)
    print("mv *.pdf *.ord *.mag.1 *.coo *.log *.shift *.out *.old *.eps *.par"
          " list_mag *.dat *.dao dat.* *.nfo {0}\n".format(wd))
    input('Press enter when finished.')
