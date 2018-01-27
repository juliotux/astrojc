import json
from collections import OrderedDict
import copy
import six
import re

from astropy.io import fits  # , ascii
# from astropy import units as u
from astropy.table import Table
# from astropy.time import Time
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.wcs import WCS
import numpy as np

from ..logging import log as logger

try:
    from . import photutils_wrapper as phot
    use_phot = True
except ModuleNotFoundError:
    use_phot = False
    logger.warn('Photutils not found, ignoring it.')

try:
    from . import sep_wrapper as sep
    use_sep = True
except ModuleNotFoundError:
    use_sep = False
    logger.warn('SEP not found, ignoring it')

from . import ccdproc_wrapper as ccdproc
from .image_shifts import ccddata_shift_images
from .catalog_wrapper import (Catalog, solve_photometry_montecarlo,
                              solve_photometry_median,
                              solve_photometry_average)
from .astrometry_wrapper import solve_astrometry_xy, wcs_radec2xy, wcs_xy2radec

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
save_calib_path         # path to save calibrated images

#pipeline config
photometry_type         # aperture or psf or both
psf_model               # model name of the psf: 'gaussian' or 'moffat'
r                       # aperture radius or list of apertures (the best snr
                        # will be used)
r_in                    # internal sky subtract radius
r_out                   # external sky subtract radius
psf_niters              # number of iterations of daofind
box_size                # box size to extract and fit stars with psf
detect_fwhm             # fwhm of daofind detect
detect_snr              # minimum signal to noise ratio to detect stars
remove_cosmics          # bool: remove cosmic rays with astroscrappy
plot                    # plot results
solve_photometry_type   # 'montecarlo', 'median', 'average'
montecarlo_iters        # number of iterations of montecarlo magnitude
montecarlo_percentage   # percentage of star in each montecarlo magnitude
identify_catalog_file   # json file containing the catalog definitions
identify_catalog_name   # name of the definition dataset in json file
identify_limit_angle    # maximum distance to identify stars in arcseconds
lightcurve_shift        # if calculate and apply shift in lightcurve
science_catalog         # file catalog with science stars
science_id_key          # column name of the id of the star
science_ra_key          # column name of the RA of the star
science_dec_key         # column name of the DEC of the star
combine_method          # method for combine images. Can be 'median', 'average'
                        # or sum
combine_sigma           # sigma of sigma_clip when combining images
combine_align_method    # if calculate and apply shift in combine: 'fft', 'wcs'
mem_limit               # maximum memory limit
'''


_mem_limit = 1e7


def _multi_process(_func, iterator, *args, **kwargs):
    """Run a function func for all i in a iterator list."""
    return [_func(i, *args, **kwargs) for i in iterator]


def _check_ccddata(image):
    """Check if a image is a CCDData. If not, try to convert it to CCDData."""
    if not isinstance(image, ccdproc.CCDData):
        if isinstance(image, six.string_types):
            return ccdproc.read_fits(image)
        elif isinstance(image, (fits.HDUList)):
            return ccdproc.CCDData(image[0].data, meta=image[0].header)
        elif isinstance(image, (fits.ImageHDU, fits.PrimaryHDU,
                                fits.ComImageHDU)):
            return ccdproc.CCDData(image.data, meta=image.header)
        else:
            raise ValueError('image type {} not supported'.format(type(image)))

    return image


def _check_iterable(value):
    """Check if a value is iterable (list), but not a string."""
    try:
        iter(value)
        if not isinstance(value, six.string_types):
            return True
        else:
            return False
    except Exception as e:
        pass

    return False


def create_calib(sources, result_file, calib_type=None, master_bias=None,
                 master_flat=None, dark_frame=None, badpixmask=None,
                 prebin=None, gain_key=None, rdnoise_key=None, gain=None,
                 combine_method='median', combine_sigma=None,
                 exposure_key=None, mem_limit=_mem_limit):
    """Create calibration frames."""
    s = _multi_process(_check_ccddata, sources)

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


def _calib_image(image, output_file=None, master_bias=None, master_flat=None,
                 dark_frame=None, badpixmask=None, prebin=None, gain=None,
                 gain_key='GAIN', rdnoise_key='RDNOISE', exposure_key=None,
                 save_calib_path=None):
    """Calib one single science image with calibration frames."""
    # TODO: save calibrated file. Problem: get the file name if image!=str
    image = _check_ccddata(image)

    if prebin is not None:
        image = ccdproc.ccdproc.block_reduce(image, prebin)
        if rdnoise_key is not None:
            image.header[rdnoise_key] = float(image.header[rdnoise_key])*prebin

    image = ccdproc.process_image(image, master_bias=master_bias,
                                  master_flat=master_flat, gain=gain,
                                  gain_key=gain_key, readnoise_key=rdnoise_key,
                                  exposure_key=exposure_key,
                                  badpixmask=badpixmask)

    return image


def calib_science(sources, master_bias=None, master_flat=None, dark_frame=None,
                  badpixmask=None, prebin=None, gain_key='GAIN', gain=None,
                  rdnoise_key='RDNOISE', combine_method=None,
                  combine_sigma=None, exposure_key=None, mem_limit=_mem_limit,
                  save_calib_path=None, combine_align_method=None):
    """Calib science images with gain, flat, bias and binning."""
    s = _multi_process(_calib_image, sources, master_bias=master_bias,
                       master_flat=master_flat, dark_frame=dark_frame,
                       badpixmask=badpixmask, prebin=prebin, gain=gain,
                       rdnoise_key=rdnoise_key, exposure_key=exposure_key)

    if combine_method is not None:
        if combine_align_method is not None:
            s = ccddata_shift_images(s, combine_align_method)
        if combine_sigma is not None:
            sig_clip = True
        s = ccdproc.combine(s, method=combine_method, sigma_clip=sig_clip,
                            sigma_clip_low_thresh=combine_sigma,
                            sigma_clip_high_thresh=combine_sigma,
                            mem_limit=mem_limit)

    return s


def process_photometry(image, photometry_type, detect_fwhm=None,
                       detect_snr=None, box_size=None,
                       r=np.arange(2, 20, 1), r_in=50,
                       r_out=60, r_find_best=True, psf_model='gaussian',
                       psf_niters=1):
    """Process standart photometry in one image, without calibrations."""
    image = _check_ccddata(image)
    data = image.data
    result = {'aperture': None,
              'psf': None}

    # aperture photometry
    if photometry_type == 'aperture' or photometry_type == 'both':
        if use_sep:
            detect = sep.find_sources
            detect_kwargs = {}
            background = sep.calculate_background
            aperture = sep.aperture_photometry
        elif use_phot:
            detect = phot.find_sources
            detect_kwargs = {'method': 'daofind', 'fwhm': detect_fwhm}
            background = phot.calculate_background
            aperture = phot.aperture_photometry
        else:
            raise ValueError('Sep and Photutils aren\'t installed. You must'
                             ' have at last one of them.')
        sky, rms = background(data)
        s = detect(data, bkg=sky, rms=rms, snr=detect_snr, **detect_kwargs)
        if _check_iterable(r):
            apertures = Table()
            for i in r:
                ap = aperture(data, s['x'], s['y'], i, r_in, r_out,
                              err=rms)
                apertures['flux_r{}'.format(i)] = ap['flux']
                apertures['flux_r{}_error'.format(i)] = ap['flux_error']

            if r_find_best:
                snr = np.zeros(len(r))
                for i in range(len(r)):
                    flux = 'flux_r{}'.format(r[i])
                    error = 'flux_r{}_error'.format(r[i])
                    snr[i] = np.median(apertures[flux]/apertures[error])
                arg = np.argmax(snr)
                logger.debug('Best aperture matched: {}'.format(r[arg]))
                apertures['flux'] = apertures['flux_r{}'.format(r[arg])]
                apertures['flux_error'] = apertures['flux_r{}'.format(r[arg]) +
                                                    '_error']
                apertures = apertures[sorted(apertures.colnames)]
        else:
            ap = aperture(data, s['x'], s['y'], i, r_in, r_out, err=rms)
            apertures = Table([ap['flux'], ap['flux_error']],
                              names=('flux', 'flux_error'))
        res_ap = Table()
        res_ap['x'] = s['x']
        res_ap['y'] = s['y']
        for i in apertures.colnames:
            res_ap[i] = apertures[i]

        result['aperture'] = res_ap

    if photometry_type == 'psf' or photometry_type == 'both':
        if not use_phot:
            raise ValueError('You must have Photutils installed for psf.')

        if photometry_type == 'both':
            x = s['x']
            y = s['y']
        else:
            x = None
            y = None

        sigma = detect_fwhm/gaussian_sigma_to_fwhm
        ph = phot.psf_photometry(data, x, y, sigma_psf=sigma, snr=detect_snr,
                                 box_size=box_size, model=psf_model,
                                 niters=psf_niters)

        result['psf'] = Table(ph)

    return result


def process_calib_photometry(image, identify_catalog,
                             identify_catalog_name=None,
                             identify_limit_angle='2 arcsec',
                             science_catalog=None, science_ra_key=None,
                             science_dec_key=None, science_id_key=None,
                             montecarlo_iters=100,
                             montecarlo_percentage=0.5, filter=None,
                             solve_photometry_type=None, **kwargs):
    """Process photometry with magnitude calibration using catalogs."""
    image = _check_ccddata

    result = {'aperture': None, 'psf': None}

    photkwargs = {}
    for i in kwargs.keys():
        if i in ['photometry_type', 'detect_fwhm', 'detect_snr', 'box_size',
                 'r', 'r_out', 'r_find_best', 'psf_model', 'psf_niters']:
            photkwargs[i] = kwargs.get(i)

    ph = process_photometry(image, **photkwargs)

    header = image.header
    wcs = WCS(header, relax=True)
    if not wcs.wcs.ctype[0]:
        imw, imh = image.data.shape
        if ph['aperture'] is not None:
            x, y = ph['aperture']['x'], ph['aperture']['y']
            flux = ph['aperture']['flux']
        elif ph['psf'] is not None:
            x, y = ph['psf']['x'], ph['psf']['y']
            flux = ph['psf']['flux']
        else:
            raise ValueError('No photometry was processed!')
        wcs = solve_astrometry_xy(x, y, flux, header, imw, imh)

    cat = Catalog.load_from_json(identify_catalog, identify_catalog_name)

    if solve_photometry_type == 'montecarlo':
        solver = solve_photometry_montecarlo
        solver_kwargs = {'n_iter': montecarlo_iters,
                         'n_star': montecarlo_percentage}
    elif solve_photometry_type == 'median':
        solver = solve_photometry_median
        solver_kwargs = {}
    elif solve_photometry_type == 'average':
        solver = solve_photometry_average
        solver_kwargs = {}
    else:
        raise ValueError('solve_photometry_type {} not'
                         ' supported'.format(solve_photometry_type))

    if ph['aperture'] is not None:
        ap = ph['aperture']
        x, y = ap['x'], ap['y']
        ra, dec = wcs_xy2radec(x, y, wcs)
        name, mag, mag_err = cat.query_id_mag(ra, dec, filter,
                                              limit_angle=identify_limit_angle)

        mags = Table()

        if 'flux' in ap.colnames:
            mags['mag'], mags['mag_err'] = solver(ap['flux'], ap['flux_error'],
                                                  mag, **solver_kwargs)

        for i in ap.colnames:
            m = re.match('^flux_r\d+[\.\d]+$', i)
            if m:
                r = re.findall('\d+[\.\d]+', i)[0]
                ma, err = solver(ap['flux_r{}'.format(r)],
                                 ap['flux+r{}_error'.format(r)],
                                 mag, **solver_kwargs)
                mags['mag_r{}'.format(r)] = ma
                mags['mag_r{}_err'.format(r)] = err

        res = Table()

        if science_catalog is not None:
            sci = Catalog.load_from_ascii(science_catalog,
                                          id_key=science_id_key,
                                          ra_key=science_ra_key,
                                          dec_key=science_dec_key,
                                          flux_key=None,
                                          flux_error_key=None,
                                          flux_unit=None)
            res['sci_id'], _, _ = sci.query(ra, dec, None,
                                            limit_angle=identify_limit_angle)

        res['cat_id'] = name
        res['cat_mag'] = mag
        res['cat_mag_err'] = mag_err

        for i in ap.colnames:
            res[i] = ap[i]
        for i in mags.colnames:
            res[i] = mags[i]

        result['aperture'] = res

    if ph['psf'] is not None:
        psf = ph['psf']
        x, y = psf['x'], psf['y']
        ra, dec = wcs_xy2radec(x, y, wcs)
        name, mag, mag_err = cat.query_id_mag(ra, dec, filter,
                                              limit_angle=identify_limit_angle)
        psf['mag'], psf['mag_err'] = solver(psf['flux'], psf['flux_error'],
                                            mag, **solver_kwargs)

        res = Table()

        if science_catalog is not None:
            sci = Catalog.load_from_ascii(science_catalog,
                                          id_key=science_id_key,
                                          ra_key=science_ra_key,
                                          dec_key=science_dec_key,
                                          flux_key=None,
                                          flux_error_key=None,
                                          flux_unit=None)
            res['sci_id'], _, _ = sci.query(ra, dec, None,
                                            limit_angle=identify_limit_angle)

            res['cat_id'] = name
            res['cat_mag'] = mag
            res['cat_mag_err'] = mag_err

            for i in psf.colnames:
                res[i] = psf[i]

        result['psf'] = res

    return result


def process_light_curve(image_set,
                        jd_key='JD',
                        shift=True,
                        **kwargs):
    """Process photometry in different files and make lightcurve with them."""


def process_polarimetry(image_set, **kwargs):
    """Process the photometry and polarimetry of a set of images."""


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
