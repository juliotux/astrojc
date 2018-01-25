import json
from collections import OrderedDict
import copy

from astropy.io import fits, ascii
from astropy import units as u
from astropy.logger import log as logger
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.time import Time
import numpy as np

from ..logging import log as logger
from . import photutils_wrapper as phot
from . import sep_wrapper as sep
from . import ccdproc_wrapper as ccdproc

# Important: this is the default keys of the configuration files!
'''
#General instrument config
ra_key                  # Right Ascension header key
dec_key                 # Declination header key
gain_key                # Instrument gain header key (e-/adu)
rdnoise_key             # Read noise header key
lamina_key              # Polarimetry retarder position key
exposure_key            # Key of the exposure time of the image

jd_key                  # Julian day header key
time_key                # Time header key will be used
time_key_type           # Type of time stored in time_key: 'jd', 'isot', 'mjd'

#Per product config
pipeline                # pipeline name: 'photometry', 'lightcurve',
                        # 'polarimetry', 'calib'
calib_type              # 'flat', 'bias', 'dark', 'science'
result_file             # file to store the results
sources                 # source files to process
prebin                  # number of pixels to bin
filter                  # filter of the image
master_bias             # bias file to correct the image
master_flat             # flat file to correct the image
badpixmask              # bad pixel mas file. pipeline=calib will created it
plate_scale             # the plate scale of the final image
gain                    # manual set of the gain of the image,
                        # if gain_key is not set

#pipeline config
photometry_type         # aperture or psf or both
psf_model               # model name of the psf: 'gaussian' or 'moffat'
r                       # aperture radius or list of apertures (the best snr
                        # will be used)
r_in                    # internal sky subtract radius
r_out                   # external sky subtract radius
box_size                # box size to extract and fit stars with psf
detect_sigma            # sigma of daofind detect
detect_snr              # minimum signal to noise ratio to detect stars
remove_cosmics          # bool: remove cosmic rays with astroscrappy
plot                    # plot results
montecarlo_iters        # number of iterations of montecarlo magnitude
montecarlo_percentage   # percentage of star in each montecarlo magnitude
identify_catalog_file   # json file containing the catalog definitions
identify_catalog_name   # name of the definition dataset in json file
identify_limit_angle    # maximum distance to identify stars in arcseconds
shift                   # if calculate fft_shift and apply the corrections
science_catalog         # file catalog with science stars
combine_method          # method for combine images. Can be 'median', 'average'
                        # or sum
combine_sigma           # sigma of sigma_clip when combining images
mem_limit               # maximum memory limit
'''


_mem_limit = 1e7


def _multi_process(_func, iterator, *args, **kwargs):
    """Run a function func for all i in a iterator list."""
    return [_func(i, *args, **kwargs) for i in iterator]


def create_calib(sources, result_file, calib_type=None, master_bias=None,
                 master_flat=None, dark_frame=None, badpixmask=None,
                 prebin=None, gain_key=None, rdnoise_key=None, gain=None,
                 combine_method='median', combine_sigma=None,
                 exposure_key=None, mem_limit=_mem_limit):
    """Create calibration frames."""
    s = _multi_process(ccdproc.read_fits, sources)

    if prebin is not None:
        s = _multi_process(ccdproc.ccdproc.block_reduce, s, prebin,
                           func=np.sum)
        if rdnoise_key is not None:
            def set_rdnoise(f, key, binning):
                f.header[key] = float(f.header[key])*binning
                return f
            s = _multi_process(set_rdnoise, s, rdnoise_key, prebin)

    s = _multi_process(ccdproc.process_image, s, master_bias=master_bias,
                       master_flat=master_flat, gain=gain, gain_key=gain_key,
                       readnoise_key=rdnoise_key, exposure_key=exposure_key)

    if combine_sigma is not None:
        sig_clip = True
    else:
        sig_clip = False

    if calib_type == 'flat':
        def scaling_func(arr):
            return 1/np.ma.average(arr)
    else:
        scaling_func = None

    res = ccdproc.combine(s, output_file=result_file,
                          method=combine_method, sigma_clip=sig_clip,
                          sigma_clip_low_thresh=combine_sigma,
                          sigma_clip_high_thresh=combine_sigma,
                          mem_limit=mem_limit, scale=scaling_func)

    if badpixmask is not None and calib_type == 'flat':
        badpix = np.logical_or(res.data <= 0.2, res.data >= 10).astype('int16')
        badpix = ccdproc.CCDData(badpix, unit='')
        badpix.write(badpixmask)

    return res


def calib_science(sources,
                  result_file,
                  master_bias=None,
                  master_flat=None,
                  badpixmask=None,
                  prebin=1,
                  gain_key='GAIN',
                  rdnoise_key='RDNOISE'):
    """Calib science images with gain, flat, bias and binning."""


def process_photometry(image,
                       detect_sigma,
                       detect_snr,
                       box_size,
                       photometry_type,
                       r=np.arange(2, 20, 1),
                       r_in=50,
                       r_out=60,
                       r_find_best=True,
                       psf_model='gaussian'):
    """Process standart photometry in one image, without calibrations."""


def process_calib_photometry(image,
                             identify_catalog,
                             identify_catalog_name=None,
                             identify_limit_angle='2 arcsec',
                             science_catalog=None,
                             montecarlo_iters=100,
                             montecarlo_percentage=0.5,
                             **kwargs):
    """Process photometry with magnitude calibration using catalogs."""


def process_light_curve(image_set,
                        jd_key='JD',
                        shift=True,
                        **kwargs):
    """Process photometry in different files and make lightcurve with them."""


class ReduceScript():
    """Simple class to process pipeline scripts qith configuration files."""
    def __init__(self, config=None):
        """Dummy entry for the class"""
        self._default_kwargss = OrderedDict()

        if config is not None:
            if isinstance(config, (dict, OrderedDict)):
                self._default_kwargss.update(config)
            elif isinstance(config, str):
                self.load_default_file(config)
            else:
                raise ValueError('Config not supported.')

    def load_default_file(self, filename):
        """Loads a json file with default variables for the reduction.
        If the configuration shares a variable with detaful values, it will be
        overrided."""
        # TODO: make it safer
        j = json.load(open(filename, 'r'))
        self._default_kwargss.update(j)

    def process_product(self, filename, dataset='all'):
        """Process the products of a json file.
        If the configuration shares a variable with detaful values, it will be
        overrided."""
        prod_file = json.load(open(filename, 'r'))
        for i, v in prod_file.items():
            prod = copy.copy(self._default_kwargss)
            prod.update(v)
            try:
                self.run(**prod)
            except Exception as e:
                logger.warn('Problem in the process of {} product from'
                            ' {} file. Passing it.'.format(i, filename))

    def run(self, **config):
        """Run a single product. Config is the dictionary of needed
        parameters."""
        raise NotImplementedError('This pipeline is not a valid implemented'
                                  ' pipeline!')
