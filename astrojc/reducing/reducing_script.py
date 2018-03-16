import json
from collections import OrderedDict
import copy
import six
import os
import glob

from astropy.io import fits
from astropy.table import Table
import numpy as np

from ..logging import log as logger

from .photometry import process_calib_photometry
from .polarimetry import run_pccdpack, process_polarimetry
from .image_shifts import ccddata_shift_images
from .image_processing import process_image, check_hdu, combine
from ..io.mkdir import mkdir_p
from ..py_utils import process_list, batch_key_replace, check_iterable


# Important: this is the default keys of the configuration files!
'''
#General instrument config
ra_key                  # Right Ascension header key
dec_key                 # Declination header key
gain_key                # Instrument gain header key (e-/adu)
rdnoise_key             # Read noise header key
exposure_key            # Key of the exposure time of the image
jd_key                  # Julian day header key
time_key                # Time header key will be used
time_key_type           # Type of time stored in time_key: 'jd', 'isot', 'mjd'
retarder_key            # Polarimetry retarder position key
retarder_type           # type of the retarder: 'half', 'quarter'
retarder_rotation       # rotation of each position of the retarder, in deg
                        # retarder, like 16
retarder_direction      # direction of the rotation of the retarder.
                        # can be 'cw', 'ccw' or -1, +1

#Per product config
pipeline                # pipeline name: 'photometry', 'lightcurve',
                        # 'polarimetry', 'calib'
calib_type              # 'flat', 'bias', 'dark', 'science'
result_file             # file to store the results
save_calib_path         # alternative folder to save calibrated images
sources                 # source files to process
source_ls_pattern       # filename pattern to load sources from ls
prebin                  # number of pixels to bin
filter                  # filter of the image
master_bias             # bias file to correct the image
master_flat             # flat file to correct the image
badpixmask              # bad pixel mas file. pipeline=calib will created it
plate_scale             # the plate scale of the final image
gain                    # manual set of the gain of the image,
                        # if gain_key is not set

#pipeline config
raw_dir                 # directory containing raw data
product_dir             # directory to store product data
calib_dir               # directory to store calib data
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
align_images            # if calculate and apply shift in lightcurve and
                        # polarimetry
solve_photometry_type   # 'montecarlo', 'median', 'average'
match_pairs_tolerance   # maximum tolerance for matching pairs for polarimetry,
                        # in pixels
montecarlo_iters        # number of iterations of montecarlo magnitude
montecarlo_percentage   # percentage of star in each montecarlo magnitude
identify_catalog_file   # json file containing the catalog definitions
identify_catalog_name   # name of the definition dataset in json file
identify_limit_angle    # maximum distance to identify stars in arcseconds
science_catalog         # file catalog with science stars
science_id_key          # column name of the id of the star
science_ra_key          # column name of the RA of the star
science_dec_key         # column name of the DEC of the star
combine_method          # method for combine images. Can be 'median', 'average'
                        # or 'sum'
combine_sigma           # sigma of sigma_clip when combining images
combine_align_method    # if calculate and apply shift in combine: 'fft', 'wcs'
mem_limit               # maximum memory limit

# In the values, {key_name} will be formated to the key_name value
'''


_mem_limit = 1e9


def _calib_image(image, product_dir=None, master_bias=None, master_flat=None,
                 dark_frame=None, badpixmask=None, prebin=None, gain=None,
                 gain_key='GAIN', rdnoise_key='RDNOISE', exposure_key=None,
                 calib_dir=None, remove_cosmics=True,
                 bias_check_keys=[], flat_check_keys=[],
                 dark_check_keys=[]):
    """Calib one single science image with calibration frames."""
    image = check_hdu(image)

    kwargs = {'rebin_size': prebin,
              'gain_key': gain_key,
              'gain': gain,
              'readnoise_key': rdnoise_key,
              'exposure_key': exposure_key,
              'lacosmic': remove_cosmics,
              'bias_check_keys': bias_check_keys,
              'flat_check_keys': flat_check_keys,
              'dark_check_keys': dark_check_keys}

    if master_bias:
        master_bias = os.path.join(calib_dir, master_bias)
        kwargs['master_bias'] = master_bias
    if master_flat:
        master_flat = os.path.join(calib_dir, master_flat)
        kwargs['master_flat'] = master_flat
    if dark_frame:
        dark_frame = os.path.join(calib_dir, dark_frame)
        kwargs['dark_frame'] = dark_frame
    if badpixmask:
        badpixmask = os.path.join(calib_dir, badpixmask)

    image = process_image(image, **kwargs)

    return image


def create_calib(sources, result_file=None, calib_type=None, master_bias=None,
                 master_flat=None, dark_frame=None, badpixmask=None,
                 prebin=None, gain_key=None, rdnoise_key=None, gain=None,
                 combine_method='median', combine_sigma=None,
                 combine_align_method=None, remove_cosmics=False,
                 exposure_key=None, mem_limit=_mem_limit, calib_dir=None,
                 bias_check_keys=[], flat_check_keys=[],
                 dark_check_keys=[]):
    """Create calibration frames."""
    s = process_list(_calib_image, sources, master_bias=master_bias,
                     master_flat=master_flat, dark_frame=dark_frame,
                     badpixmask=badpixmask, prebin=prebin, gain_key=gain_key,
                     rdnoise_key=rdnoise_key, gain=gain,
                     exposure_key=exposure_key, remove_cosmics=remove_cosmics,
                     calib_dir=calib_dir,
                     bias_check_keys=bias_check_keys,
                     flat_check_keys=flat_check_keys,
                     dark_check_keys=dark_check_keys)

    if combine_sigma is not None:
        reject = 'sigmaclip'
    else:
        reject = []

    if calib_type == 'flat':
        def scaling_func(arr):
            return 1/np.ma.average(arr)
    else:
        scaling_func = None

    if result_file is not None and calib_dir is not None:
        result_file = os.path.join(calib_dir, result_file)

    if combine_align_method in ['fft', 'wcs']:
        s = ccddata_shift_images(s, combine_align_method)

    res = combine(s, output_file=result_file, method=combine_method,
                  reject=reject, sigma_clip_low=combine_sigma,
                  sigma_clip_high=combine_sigma, mem_limit=mem_limit,
                  scale=scaling_func)

    if badpixmask is not None and calib_type == 'flat':
        badpix = np.logical_or(res.data <= 0.2, res.data >= 10).astype('uint8')
        badpix = check_hdu(badpix)
        if calib_dir is not None:
            badpixmask = os.path.join(calib_dir, badpixmask)
        badpix.writeto(badpixmask)

    return res


def calib_science(sources, master_bias=None, master_flat=None, dark_frame=None,
                  badpixmask=None, prebin=None, gain_key='GAIN', gain=None,
                  rdnoise_key='RDNOISE', combine_method=None,
                  combine_sigma=None, exposure_key=None, mem_limit=_mem_limit,
                  product_dir=None, combine_align_method=None, calib_dir=None,
                  save_calib_path=None, remove_cosmics=True,
                  return_filename=False, bias_check_keys=[],
                  flat_check_keys=[], dark_check_keys=[]):
    """Calib science images with gain, flat, bias and binning."""
    s = process_list(_calib_image, sources, master_bias=master_bias,
                     master_flat=master_flat, dark_frame=dark_frame,
                     badpixmask=badpixmask, prebin=prebin, gain=gain,
                     rdnoise_key=rdnoise_key, exposure_key=exposure_key,
                     calib_dir=calib_dir, remove_cosmics=remove_cosmics,
                     bias_check_keys=bias_check_keys,
                     flat_check_keys=flat_check_keys,
                     dark_check_keys=dark_check_keys)

    s = process_list(check_hdu, s)

    if combine_align_method in ['fft', 'wcs']:
        s = ccddata_shift_images(s, method=combine_align_method)

    if save_calib_path is not None and isinstance(sources[0],
                                                  six.string_types):
        mkdir_p(save_calib_path)
        names = process_list(lambda x: os.path.join(save_calib_path,
                                                    os.path.basename(x)),
                             sources)
        process_list(lambda x: fits.writeto(*x, overwrite=True),
                     zip(names, [i.data for i in s], [i.header for i in s]))
        files_saved = True
    else:
        files_saved = False

    if combine_method is not None:
        if combine_align_method is not None:
            s = ccddata_shift_images(s, combine_align_method)
        if combine_sigma is not None:
            reject = ['sigmaclip']
        else:
            reject = []
        s = combine(s, method=combine_method, reject=reject,
                    sigma_clip_low=combine_sigma,
                    sigma_clip_high=combine_sigma,
                    mem_limit=mem_limit)

    if return_filename and files_saved:
        return names

    return s


class ReduceScript():
    """Simple class to process pipeline scripts with configuration files."""
    # TODO: introduce parameters checking
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

    def clear(self):
        self._default_kwargss = OrderedDict()

    def process_product(self, filename, dataset='all', raise_error=False,
                        **kwargs):
        """Process the products of a json file.
        If the configuration shares a variable with detaful values, it will be
        overrided."""
        prod_file = json.load(open(filename, 'r'))
        default = copy.copy(self._default_kwargss)
        if '__preload__' in prod_file.keys():
            preload = copy.copy(self._default_kwargss)
            preload.update(prod_file.pop('__preload__'))
            batch_key_replace(preload)
            if 'preload_config_files' in preload.keys():
                files = preload.pop('preload_config_files')
                other = {}
                if check_iterable(files):
                    for f in files:
                        other.update(json.load(open(f, 'r')))
                else:
                    other.update(json.load(open(files, 'r')))
                other.update(preload)
                preload = other
            default.update(preload)

        if dataset not in ['all', None]:
            valid = dataset
        else:
            valid = prod_file.keys()
        for i, v in prod_file.items():
            if i in valid:
                prod = copy.copy(default)
                prod.update(v)
                prod.update(kwargs)
                batch_key_replace(prod)
                logger.info('Reducing {} from {}'.format(i, filename))
                if raise_error:
                    self.run(i, **prod)
                else:
                    try:
                        self.run(i, **prod)
                    except Exception as e:
                        logger.error('Problem in the process of {} product'
                                     ' from {} file. Passing it.'
                                     .format(i, filename) +
                                     '\nError: {}'.format(e))

    def run(self, name, **config):
        """Run a single product. Config is the dictionary of needed
        parameters."""
        raise NotImplementedError('This pipeline is not a valid implemented'
                                  ' pipeline!')

    def __call__(self, name, **config):
        return self.run(name, **config)


class CalibScript(ReduceScript):
    def __init__(self, config=None):
        super(CalibScript, self).__init__(config=config)

    def run(self, name, **config):
        for k in config.keys():
            batch_key_replace(config, k)
        s = [os.path.join(config['raw_dir'], i) for i in config['sources']]

        calib_kwargs = {}
        for i in ['calib_type', 'master_bias', 'master_flat', 'dark_frame',
                  'badpixmask', 'prebin', 'gain_key', 'rdnoise_key', 'gain',
                  'combine_method', 'combine_sigma', 'exposure_key',
                  'mem_limit', 'calib_dir',
                  'bias_check_keys', 'flat_check_keys', 'dark_check_keys']:
            if i in config.keys():
                calib_kwargs[i] = config[i]

        if 'result_file' not in config.keys():
            outname = name
        else:
            outname = config['result_file']

        mkdir_p(config['calib_dir'])
        return create_calib(s, outname, **calib_kwargs)


class PhotometryScript(ReduceScript):
    def __init__(self, config=None):
        super(PhotometryScript, self).__init__(config=config)

    def run(self, name, **config):
        """Run this pipeline script"""
        product_dir = config['product_dir']
        night = config['night']
        phot_prod = os.path.join(product_dir,
                                 "{}_photometry_{}".format(night, name))
        s = [os.path.join(config['raw_dir'], i) for i in config['sources']]

        if config.get('astrojc_cal', True):
            calib_kwargs = {}
            for i in ('master_bias', 'master_flat', 'dark_frame', 'badpixmask',
                      'prebin', 'gain_key', 'gain', 'rdnoise_key',
                      'combine_method', 'combine_sigma', 'exposure_key',
                      'mem_limit', 'save_calib_path', 'combine_align_method',
                      'calib_dir', 'product_dir', 'remove_cosmics',
                      'bias_check_keys', 'flat_check_keys', 'dark_check_keys'):
                if i in config.keys():
                    calib_kwargs[i] = config[i]
            ccd = calib_science(s, **calib_kwargs)
        elif 'save_calib_path' in config.keys():
            ccd = process_list(os.path.basename, s)
            ccd = [os.path.join(config['save_calib_path'], i) for i in ccd]
            ccd = combine(ccd, method=config['combine_method'])
        else:
            ccd = combine(s, method=config['combine_method'])

        photkwargs = {}
        for i in ['ra_key', 'dec_key', 'gain_key', 'rdnoise_key',
                  'filter', 'plate_scale', 'photometry_type',
                  'psf_model', 'r', 'r_in', 'r_out', 'psf_niters',
                  'box_size', 'detect_fwhm', 'detect_snr', 'remove_cosmics',
                  'align_images', 'solve_photometry_type',
                  'montecarlo_iters', 'montecarlo_percentage',
                  'identify_catalog_file', 'identify_catalog_name',
                  'identify_limit_angle', 'science_catalog', 'science_id_key',
                  'science_ra_key', 'science_dec_key']:
            if i in config.keys():
                photkwargs[i] = config[i]

        t = process_calib_photometry(ccd, **photkwargs)

        hdus = []
        for i in [i for i in t.keys() if t[i] is not None]:
            header_keys = ['solve_photometry_type', 'plate_scale', 'filter']
            if i == 'aperture':
                header_keys += ['r', 'r_in', 'r_out', 'detect_fwhm',
                                'detect_snr']
            elif i == 'psf':
                header_keys += ['psf_model', 'box_size', 'psf_niters']

            if config.get('solve_photometry_type', None) == 'montecarlo':
                header_keys += ['montecarlo_iters', 'montecarlo_percentage']

            if config.get('identify_catalog_name', None) is not None:
                header_keys += ['identify_catalog_name',
                                'identify_limit_angle']

            hdu = fits.BinTableHDU(t[i], name="{}_photometry".format(i))
            for k in header_keys:
                if k in config.keys():
                    v = config[k]
                    key = 'hierarch astrojc {}'.format(k)
                    if check_iterable(v):
                        hdu.header[key] = ','.join([str(m) for m in v])
                    else:
                        hdu.header[key] = v
            hdus.append(hdu)

            best = Table(dtype=t[i].dtype)
            for group in t[i].group_by('star_index').groups:
                b = np.argmax(group['flux']/group['flux_error'])
                best.add_row(group[b])
            best['snr'] = best['flux']/best['flux_error']

            hdu = fits.BinTableHDU(best, name="{}_best_snr".format(i))
            for k in header_keys:
                if k in config.keys():
                    v = config[k]
                    key = 'hierarch astrojc {}'.format(k)
                    if check_iterable(v):
                        hdu.header[key] = ','.join([str(m) for m in v])
                    else:
                        hdu.header[key] = v
            hdus.append(hdu)

        mkdir_p(product_dir)
        hdulist = fits.HDUList([ccd, *hdus])
        hdulist.writeto(phot_prod, overwrite=True)


class PolarimetryScript(ReduceScript):
    def __init__(self, config=None):
        super(PolarimetryScript, self).__init__(config=config)

    def run(self, name, **config):
        """Run this pipeline script"""
        product_dir = config['product_dir']
        night = config['night']
        s = [os.path.join(config['raw_dir'], i) for i in config['sources']]

        check_exist = config.get('check_exist', False)

        astrojc_prod = os.path.join(product_dir, "{}_polarimetry_astrojc_{}"
                                    .format(night, name))
        pccd_prod = os.path.join(product_dir, "{}_polarimetry_pccdpack_{}"
                                 .format(night, name))
        if check_exist:
            process_astrojc = (config.get('astrojc_pol', True) and not
                               os.path.isfile(astrojc_prod))
            process_pccd = (config.get('pccdpack', False) and not
                            os.path.isfile(pccd_prod))
        else:
            process_pccd = config.get('pccdpack', False)
            process_astrojc = config.get('astrojc_pol', True)

        if not process_pccd and not process_astrojc:
            return

        calib_kwargs = {}
        for i in ('master_bias', 'master_flat', 'dark_frame', 'badpixmask',
                  'prebin', 'gain_key', 'gain', 'rdnoise_key',
                  'combine_method', 'combine_sigma', 'exposure_key',
                  'mem_limit', 'combine_align_method',
                  'calib_dir', 'product_dir', 'remove_cosmics',
                  'save_calib_path'):
            if i in config.keys():
                calib_kwargs[i] = config[i]

        if config.get('astrojc_cal', True) and (process_astrojc or
                                                process_pccd):
            ccds = calib_science(s, **calib_kwargs)
        else:
            ccds = process_list(check_hdu, s)

        polkwargs = {}
        for i in ['ra_key', 'dec_key', 'gain_key', 'rdnoise_key',
                  'retarder_key', 'retarder_type', 'retarder_direction',
                  'filter', 'plate_scale', 'photometry_type',
                  'psf_model', 'r', 'r_in', 'r_out', 'psf_niters',
                  'box_size', 'detect_fwhm', 'detect_snr', 'remove_cosmics',
                  'align_images', 'solve_photometry_type',
                  'match_pairs_tolerance', 'montecarlo_iters',
                  'montecarlo_percentage', 'identify_catalog_file',
                  'identify_catalog_name', 'identify_limit_angle',
                  'science_catalog', 'science_id_key', 'science_ra_key',
                  'science_dec_key', 'astrometry_calib',
                  'delta_x', 'delta_y', 'brightest_star_dec',
                  'brightest_star_ra', 'image_flip', 'image_north_direction']:
            if i in config.keys():
                polkwargs[i] = config[i]

        if process_astrojc:
            if not config.get('astrojc_cal', True):
                if 'save_calib_path' in config.keys():
                    ccds = [os.path.join(config['save_calib_path'],
                                         os.path.basename(i))
                            for i in s]
                    ccds = process_list(check_hdu, ccds)
            logger.info('Processing polarimetry with astrojc.')
            logger.debug('Processing {} images'.format(len(ccds)))
            t, wcs, ret = process_polarimetry(ccds, **polkwargs)
            config['retarder_positions'] = ret
        else:
            wcs = None
            t = {}

        mkdir_p(product_dir)

        image = combine(ccds, method='sum', mem_limit=config.get('mem_limit',
                                                                 _mem_limit))

        hdus = []
        for i in [i for i in t.keys() if t[i] is not None]:
            header_keys = ['retarder_type', 'retarder_rotation',
                           'retarder_direction', 'retarder_positions',
                           'align_images', 'solve_photometry_type',
                           'plate_scale', 'filter', 'night']
            if i == 'aperture':
                header_keys += ['r', 'r_in', 'r_out', 'detect_fwhm',
                                'detect_snr']
            elif i == 'psf':
                header_keys += ['psf_model', 'box_size', 'psf_niters']

            if config.get('solve_photometry_type', None) == 'montecarlo':
                header_keys += ['montecarlo_iters', 'montecarlo_percentage']

            if config.get('identify_catalog_name', None) is not None:
                header_keys += ['identify_catalog_name',
                                'identify_limit_angle']

            hdu = fits.BinTableHDU(t[i], name="{}_log".format(i))
            for k in header_keys:
                if k in config.keys():
                    v = config[k]
                    key = 'hierarch astrojc {}'.format(k)
                    if check_iterable(v):
                        hdu.header[key] = ','.join([str(m) for m in v])
                    else:
                        hdu.header[key] = v
            hdus.append(hdu)

            out = Table(dtype=t[i].dtype)
            for group in t[i].group_by('star_index').groups:
                m = np.argmax(group['p']/group['p_error'])
                out.add_row(group[m])
            out['snr'] = out['p']/out['p_error']

            hdu = fits.BinTableHDU(out, name="{}_out".format(i))
            for k in header_keys:
                if k in config.keys():
                    v = config[k]
                    key = 'hierarch astrojc {}'.format(k)
                    if check_iterable(v):
                        hdu.header[key] = ','.join([str(m) for m in v])
                    else:
                        hdu.header[key] = v
            hdus.append(hdu)

        if process_astrojc:
            if wcs is not None:
                image.header.update(wcs.to_header(relax=True))
            hdulist = fits.HDUList([image, *hdus])
            hdulist.writeto(astrojc_prod, overwrite=True)

        if process_pccd:
            if not config.get('astrojc_cal', True):
                if 'save_calib_path' in config.keys():
                    ccds = [os.path.join(config['save_calib_path'],
                                         os.path.basename(i))
                            for i in s]
                    ccds = process_list(check_hdu, ccds)
            logger.info('Processing polarimetry with pccdpack.')
            pccd = run_pccdpack(ccds, wcs=wcs, **polkwargs)
            wcs = pccd[3]
            hdus = []
            hdus.append(fits.BinTableHDU(pccd[0], name='out_table'))
            hdus.append(fits.BinTableHDU(pccd[1], name='dat_table'))
            hdus.append(fits.BinTableHDU(pccd[2], name='log_table'))
            header_keys = ['retarder_type', 'retarder_rotation',
                           'retarder_direction', 'retarder_positions',
                           'align_images', 'solve_photometry_type',
                           'plate_scale', 'filter', 'night',
                           'r', 'r_in', 'r_out']
            if config.get('identify_catalog_name', None) is not None:
                header_keys += ['identify_catalog_name',
                                'identify_limit_angle']
            for i in hdus:
                for k in header_keys:
                    if k in config.keys():
                        v = config[k]
                        key = 'hierarch astrojc {}'.format(k)
                        if check_iterable(v):
                            hdu.header[key] = ','.join([str(m) for m in v])
                        else:
                            hdu.header[key] = v
            hdulist = fits.HDUList([image, *hdus])
            hdulist.writeto(pccd_prod, overwrite=True)


class MasterReduceScript(ReduceScript):
    def __init__(self, config=None):
        super(MasterReduceScript, self).__init__(config=config)

    def run(self, name, **config):
        if 'sources' in config.keys() and 'source_ls_pattern' in config.keys():
            logger.warn('sources and sources_ls_pattern given. Using sources.')
        elif 'sources_ls_pattern' in config.keys():
            ls_pat = os.path.join(config['raw_dir'],
                                  config['sources_ls_pattern'])
            fs = glob.glob(ls_pat)
            config['sources'] = [os.path.basename(i) for i in sorted(fs)]
            if len(config['sources']) == 0:
                raise FileNotFoundError("Could not determine sources."
                                        " glob pattern: {}".format(ls_pat))

        config['sources'] = [i for i in config['sources']
                             if i not in config['exclude_images']]

        logger.debug("Product {} config:{}".format(name, str(config)))
        if 'pipeline' not in config.keys():
            raise ValueError('The config must specify what pipeline will be'
                             ' used!')
        if config['pipeline'] == 'photometry':
            p = PhotometryScript()
        elif config['pipeline'] == 'polarimetry':
            p = PolarimetryScript()
        elif config['pipeline'] == 'calib':
            p = CalibScript()
        elif config['pipeline'] == 'lightcurve':
            logger.error('lightcurve pipeline not implemented yet')
        else:
            raise ValueError('Pipeline {} not'
                             ' supported.'.format(config['pipeline']))

        p.run(name, **config)
