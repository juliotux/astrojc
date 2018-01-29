'''Wrapper from ccdproc package.'''
import six
import ccdproc
from ccdproc.combiner import combine
from ccdproc import CCDData
import numpy as np
from astropy import units as u
from astropy.io import fits

mem_limit = 1e8


def read_fits(fname):
    try:
        return ccdproc.CCDData.read(fname)
    except ValueError:
        # assume if not gain corrected, the unit is adu
        return ccdproc.CCDData.read(fname, unit='adu')


def save_fits(ccd, fname):
    ccd.write(fname)


_ccd_procces_keys = ['oscan', 'trim', 'error', 'master_bias', 'dark_frame',
                     'master_flat', 'bad_pixel_mask', 'gain', 'readnoise',
                     'oscan_median', 'oscan_model', 'min_value',
                     'dark_exposure', 'data_exposure', 'exposure_key',
                     'exposure_unit', 'dark_scale', 'gain_corrected',
                     'add_keyword']


def check_ccddata(image):
    """Check if a image is a CCDData. If not, try to convert it to CCDData."""
    if not isinstance(image, ccdproc.CCDData):
        if isinstance(image, six.string_types):
            return read_fits(image)
        elif isinstance(image, (fits.HDUList)):
            return CCDData(image[0].data, meta=image[0].header)
        elif isinstance(image, (fits.ImageHDU, fits.PrimaryHDU,
                                fits.ComImageHDU)):
            return CCDData(image.data, meta=image.header)
        else:
            raise ValueError('image type {} not supported'.format(type(image)))

    return image


def process_image(ccd, gain_key=None, readnoise_key=None,
                  *args, **kwargs):
    nkwargs = {}
    for i in kwargs.keys():
        if i in _ccd_procces_keys:
            nkwargs[i] = kwargs.get(i)

    for i in ['master_bias', 'master_flat', 'dark_frame']:
        if i in nkwargs.keys() and \
           not isinstance(nkwargs[i], CCDData) and \
           nkwargs[i] is not None:
            nkwargs[i] = read_fits(nkwargs[i])

    if gain_key is not None:
        if isinstance(ccd, ccdproc.CCDData):
            nkwargs['gain'] = float(ccd.header[gain_key])*u.electron/u.adu
        else:
            nkwargs['gain'] = float(fits.getval(gain_key))*u.electron/u.adu

    if readnoise_key is not None:
        if isinstance(ccd, ccdproc.CCDData):
            nkwargs['readnoise'] = float(ccd.header[readnoise_key])*u.electron
        else:
            nkwargs['readnoise'] = float(fits.getval(gain_key))*u.electron

    if kwargs.get('master_bias', None):
        kwargs['master_bias'] = check_ccddata(kwargs['master_bias'])
    if kwargs.get('master_flat', None):
        kwargs['master_flat'] = check_ccddata(kwargs['master_flat'])
    if kwargs.get('dark_frame', None):
        kwargs['dark_frame'] = check_ccddata(kwargs['dark_frame'])

    return ccdproc.ccd_process(ccd, *args, **nkwargs)


def combine_bias(ccd_list, combine_method='median', mem_limit=mem_limit,
                 **reductwargs):
    for i in range(len(ccd_list)):
        ccd_list[i] = process_image(ccd_list[i], **reductwargs)
    return combine(ccd_list, method=combine_method, mem_limit=mem_limit)


def combine_flat(ccd_list, master_bias, combine_method='median',
                 mem_limit=mem_limit, **reductwargs):
    for i in range(len(ccd_list)):
        ccd_list[i] = process_image(ccd_list[i], master_bias=master_bias,
                                    **reductwargs)
        ccd_list[i].data = ccd_list[i].data/np.median(ccd_list[i].data)
    return combine(ccd_list, method=combine_method, mem_limit=mem_limit)
