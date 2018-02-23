import json
from collections import OrderedDict
import copy
import six
import re
import os
import glob
import shutil
import time

from astropy.io import fits  # , ascii
# from astropy import units as u
from astropy.table import Table, Column
# from astropy.time import Time
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.wcs import WCS
import numpy as np
import pandas as pd
from tempfile import mkdtemp, NamedTemporaryFile

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
from .ccdproc_wrapper import check_ccddata as _check_ccddata
from .image_shifts import ccddata_shift_images
from .catalog_wrapper import (Catalog, solve_photometry_montecarlo,
                              solve_photometry_median,
                              solve_photometry_average)
from .astrometry_wrapper import solve_astrometry_xy, wcs_xy2radec
from .polarimetry import estimate_dxdy, match_pairs, calculate_polarimetry
from . import pccdpack_wrapper as pccd
from ..io.mkdir import mkdir_p
from ..py_utils import process_list, batch_key_replace, check_iterable


DEBUG = True

if DEBUG:
    from matplotlib import pyplot as plt
    failed = 'fail.log'

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


def create_calib(sources, result_file=None, calib_type=None, master_bias=None,
                 master_flat=None, dark_frame=None, badpixmask=None,
                 prebin=None, gain_key=None, rdnoise_key=None, gain=None,
                 combine_method='median', combine_sigma=None,
                 combine_align_method=None,
                 exposure_key=None, mem_limit=_mem_limit, calib_dir=None):
    """Create calibration frames."""
    s = process_list(_check_ccddata, sources)

    if prebin is not None:
        s = process_list(ccdproc.ccdproc.block_reduce, s, prebin,
                         func=np.sum)
        if rdnoise_key is not None:
            def set_rdnoise(f, key, binning):
                f.header[key] = float(f.header[key])*binning
                return f
            s = process_list(set_rdnoise, s, rdnoise_key, prebin)

    if master_bias is not None:
            master_bias = os.path.join(calib_dir, master_bias)
    if master_flat is not None:
            master_flat = os.path.join(calib_dir, master_flat)

    s = process_list(ccdproc.process_image, s, master_bias=master_bias,
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

    if result_file is not None and calib_dir is not None:
        result_file = os.path.join(calib_dir, result_file)

    if combine_align_method in ['fft', 'wcs']:
        s = ccddata_shift_images(s, combine_align_method)

    res = ccdproc.combine(s, output_file=result_file,
                          method=combine_method, sigma_clip=sig_clip,
                          sigma_clip_low_thresh=combine_sigma,
                          sigma_clip_high_thresh=combine_sigma,
                          mem_limit=mem_limit, scale=scaling_func)

    if badpixmask is not None and calib_type == 'flat':
        badpix = np.logical_or(res.data <= 0.2, res.data >= 10).astype('uint8')
        badpix = ccdproc.CCDData(badpix, unit='')
        if calib_dir is not None:
            badpixmask = os.path.join(calib_dir, badpixmask)
        badpix.write(badpixmask)

    return res


def _calib_image(image, product_dir=None, master_bias=None, master_flat=None,
                 dark_frame=None, badpixmask=None, prebin=None, gain=None,
                 gain_key='GAIN', rdnoise_key='RDNOISE', exposure_key=None,
                 calib_dir=None):
    """Calib one single science image with calibration frames."""
    image = _check_ccddata(image)

    if prebin is not None:
        image = ccdproc.ccdproc.block_reduce(image, prebin)
        if rdnoise_key is not None:
            image.header[rdnoise_key] = float(image.header[rdnoise_key])*prebin

    if master_bias:
        master_bias = os.path.join(calib_dir, master_bias)
    if master_flat:
        master_flat = os.path.join(calib_dir, master_flat)
    if dark_frame:
        dark_frame = os.path.join(calib_dir, dark_frame)
    if badpixmask:
        badpixmask = os.path.join(calib_dir, badpixmask)

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
                  product_dir=None, combine_align_method=None, calib_dir=None,
                  save_calib_path=None, remove_cosmics=True,
                  return_filename=False):
    """Calib science images with gain, flat, bias and binning."""
    s = process_list(_calib_image, sources, master_bias=master_bias,
                     master_flat=master_flat, dark_frame=dark_frame,
                     badpixmask=badpixmask, prebin=prebin, gain=gain,
                     rdnoise_key=rdnoise_key, exposure_key=exposure_key,
                     calib_dir=calib_dir)

    s = process_list(_check_ccddata, s)

    if combine_align_method in ['fft', 'wcs']:
        s = ccddata_shift_images(s, method=combine_align_method)

    if remove_cosmics:
        s = process_list(ccdproc.ccdproc.cosmicray_lacosmic, s)

    if save_calib_path is not None and isinstance(sources[0],
                                                  six.string_types):
        mkdir_p(save_calib_path)
        names = process_list(lambda x: os.path.join(save_calib_path,
                                                    os.path.basename(x)),
                             sources)
        process_list(lambda x: ccdproc.save_fits(*x), zip(s, names))
        files_saved = True
    else:
        files_saved = False

    if combine_method is not None:
        if combine_align_method is not None:
            s = ccddata_shift_images(s, combine_align_method)
        if combine_sigma is not None:
            sig_clip = True
        else:
            sig_clip = False
        s = ccdproc.combine(s, method=combine_method, sigma_clip=sig_clip,
                            sigma_clip_low_thresh=combine_sigma,
                            sigma_clip_high_thresh=combine_sigma,
                            mem_limit=mem_limit)

    if return_filename and files_saved:
        return names

    return s


def _do_aperture(data, detect_fwhm=None, detect_snr=None, x=None, y=None,
                 r=np.arange(2, 20, 1), r_in=50, r_out=60, r_find_best=True):
    """Perform aperture photometry in a image"""
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
    if x is not None and y is not None:
        s = np.array(list(zip(x, y)), dtype=[('x', 'f8'), ('y', 'f8')])
    else:
        s = detect(data, bkg=sky, rms=rms, snr=detect_snr, **detect_kwargs)
    if check_iterable(r):
        apertures = Table()
        for i in r:
            # logger.debug('Aperture photometry in R={}'.format(i))
            ap = aperture(data, s['x'], s['y'], i, r_in, r_out,
                          err=rms)
            apertures['flux_r{}'.format(i)] = ap['flux']
            apertures['flux_r{}_error'.format(i)] = ap['flux_error']

        # FIXME: Allways return the smallest aperture
        if r_find_best:
            snr = np.zeros(len(r))
            for i in range(len(r)):
                flux = 'flux_r{}'.format(r[i])
                error = 'flux_r{}_error'.format(r[i])
                snr[i] = np.median(apertures[flux]/apertures[error])
            arg = np.argmax(snr)
            # logger.debug([(i, j) for (i, j) in zip(r, snr)])
            # logger.debug('Best aperture matched: {}'.format(r[arg]))
            apertures['flux'] = apertures['flux_r{}'.format(r[arg])]
            apertures['flux_error'] = apertures['flux_r{}'.format(r[arg]) +
                                                '_error']
            apertures = apertures[sorted(apertures.colnames)]
    else:
        ap = aperture(data, s['x'], s['y'], r, r_in, r_out, err=rms)
        apertures = Table([ap['flux'], ap['flux_error']],
                          names=('flux', 'flux_error'))
    res_ap = Table()
    res_ap['x'] = s['x']
    res_ap['y'] = s['y']
    for i in apertures.colnames:
        res_ap[i] = apertures[i]

    return res_ap


def process_photometry(image, photometry_type, detect_fwhm=None,
                       detect_snr=None, box_size=None,
                       r=np.arange(2, 20, 1), r_in=50,
                       r_out=60, r_find_best=True, psf_model='gaussian',
                       psf_niters=1, x=None, y=None):
    """Process standart photometry in one image, without calibrations."""
    image = _check_ccddata(image)
    data = image.data
    result = {'aperture': None,
              'psf': None}

    # aperture photometry
    if photometry_type == 'aperture' or photometry_type == 'both':
        result['aperture'] = _do_aperture(data, detect_fwhm=detect_fwhm,
                                          detect_snr=detect_snr, r=r,
                                          r_in=r_in, r_out=r_out,
                                          r_find_best=r_find_best,
                                          x=x, y=y)
        x = result['aperture']['x']
        y = result['aperture']['y']

    if photometry_type == 'psf' or photometry_type == 'both':
        if not use_phot:
            raise ValueError('You must have Photutils installed for psf.')

        sigma = detect_fwhm/gaussian_sigma_to_fwhm
        ph = phot.psf_photometry(data, x, y, sigma_psf=sigma, snr=detect_snr,
                                 box_size=box_size, model=psf_model,
                                 niters=psf_niters)

        result['psf'] = Table(ph)

    return result


def _identify_star(table, wcs, filter, identify_catalog_file,
                   identify_catalog_name=None, identify_limit_angle='2 arcsec',
                   science_catalog=None, science_id_key=None,
                   science_ra_key=None, science_dec_key=None):
    cat = Catalog.load_from_json(identify_catalog_file, identify_catalog_name)
    x, y = table['x'], table['y']
    ra, dec = wcs_xy2radec(x, y, wcs)

    name, mag, mag_err = cat.query_id_mag(ra, dec, filter,
                                          limit_angle=identify_limit_angle)

    res = Table()
    if science_catalog is not None:
        sci = Catalog.load_from_ascii(science_catalog,
                                      id_key=science_id_key,
                                      ra_key=science_ra_key,
                                      dec_key=science_dec_key,
                                      flux_key=None,
                                      flux_error_key=None,
                                      flux_unit=None,
                                      filters=None,
                                      prepend_id_key=False)
        limit_angle = identify_limit_angle
        res['sci_id'], _, _ = sci.query_id_mag(ra, dec, None,
                                               limit_angle=limit_angle)
        res['ra'] = ra
        res['dec'] = dec

    res['cat_id'] = name
    res['cat_mag'] = mag
    res['cat_mag_err'] = mag_err

    return res


def _solve_photometry(table, wcs=None, cat_mag=None,
                      identify_catalog_file=None, identify_catalog_name=None,
                      identify_limit_angle='2 arcsec', science_catalog=None,
                      science_id_key=None, science_ra_key=None,
                      science_dec_key=None, montecarlo_iters=100,
                      montecarlo_percentage=0.5, filter=None,
                      solve_photometry_type=None):
    """Solve the absolute photometry of a field using a catalog."""

    if DEBUG:
        start = time.time()

    if solve_photometry_type == 'montecarlo':
        solver = solve_photometry_montecarlo
        solver_kwargs = {'n_iter': montecarlo_iters,
                         'n_stars': montecarlo_percentage}
    elif solve_photometry_type == 'median':
        solver = solve_photometry_median
        solver_kwargs = {}
    elif solve_photometry_type == 'average':
        solver = solve_photometry_average
        solver_kwargs = {}
    else:
        raise ValueError('solve_photometry_type {} not'
                         ' supported'.format(solve_photometry_type))

    if cat_mag is None:
        id_table = _identify_star(table=table, wcs=wcs, filter=filter,
                                  identify_catalog_file=identify_catalog_file,
                                  identify_catalog_name=identify_catalog_name,
                                  science_catalog=science_catalog,
                                  science_id_key=science_id_key,
                                  science_ra_key=science_ra_key,
                                  science_dec_key=science_dec_key)
        cat_mag = id_table['cat_mag']

    mags = Table()

    if 'flux' in table.colnames:
        mags['mag'], mags['mag_err'] = solver(table['flux'],
                                              table['flux_error'],
                                              cat_mag, **solver_kwargs)

    for i in table.colnames:
        m = re.match('^flux_r\d+(\.\d+)?$', i)
        if m:
            r = i.strip('flux_r')
            ma, err = solver(table['flux_r{}'.format(r)],
                             table['flux_r{}_error'.format(r)],
                             cat_mag, **solver_kwargs)
            mags['mag_r{}'.format(r)] = ma
            mags['mag_r{}_err'.format(r)] = err

    if DEBUG:
        end = time.time()
        logger.debug("solving photometry took {} seconds".format(end-start))

    return mags


def _solve_astrometry(header, table, shape, ra_key=None, dec_key=None,
                      plate_scale=None):
    """Solves the astrometry of a field and return a valid wcs."""
    if DEBUG:
        start = time.time()
    wcs = WCS(header, relax=True)
    if not wcs.wcs.ctype[0]:
        im_params = {}
        if ra_key is not None and dec_key is not None:
            im_params['ra_key'] = ra_key
            im_params['dec_key'] = dec_key
        if plate_scale is not None:
            im_params['pltscl'] = plate_scale
            im_params['radius'] = 5*plate_scale*np.max(shape)/3600
        imw, imh = shape
        x, y = table['x'], table['y']
        # TODO: when note r_find_best, this key do not exists
        flux = table['flux']
        wcs = solve_astrometry_xy(x, y, flux, header, imw, imh,
                                  image_params=im_params, return_wcs=True)
    if DEBUG:
        end = time.time()
        logger.debug("solving astrometry took {} seconds".format(end-start))
    return wcs


def process_calib_photometry(image, identify_catalog_file=None,
                             identify_catalog_name=None,
                             identify_limit_angle='2 arcsec',
                             science_catalog=None, science_ra_key=None,
                             science_dec_key=None, science_id_key=None,
                             montecarlo_iters=100,
                             montecarlo_percentage=0.2, filter=None,
                             solve_photometry_type=None, **kwargs):
    """Process photometry with magnitude calibration using catalogs."""
    image = _check_ccddata(image)

    result = {'aperture': None, 'psf': None}

    photkwargs = {}
    for i in kwargs.keys():
        if i in ['photometry_type', 'detect_fwhm', 'detect_snr', 'box_size',
                 'r', 'r_out', 'r_find_best', 'psf_model', 'psf_niters']:
            photkwargs[i] = kwargs.get(i)

    ph = process_photometry(image, **photkwargs)

    processed = [i for i in ['aperture', 'psf'] if ph[i] is not None]
    if len(processed) == 0:
        raise ValueError('No photometry was processed!')
    else:
        wcs = _solve_astrometry(image.header, ph[processed[0]],
                                image.data.shape,
                                ra_key=kwargs['ra_key'],
                                dec_key=kwargs['dec_key'],
                                plate_scale=kwargs['plate_scale'])
        ids = _identify_star(ph[processed[0]], wcs, filter=filter,
                             identify_catalog_file=identify_catalog_file,
                             identify_catalog_name=identify_catalog_name,
                             identify_limit_angle=identify_limit_angle,
                             science_catalog=science_catalog,
                             science_id_key=science_id_key,
                             science_ra_key=science_ra_key,
                             science_dec_key=science_dec_key)

    for i in processed:
        res = _solve_photometry(ph[i], wcs, ids['cat_mag'],
                                montecarlo_iters=montecarlo_iters,
                                montecarlo_percentage=montecarlo_percentage,
                                solve_photometry_type=solve_photometry_type)
        result[i] = copy.copy(ids)
        for j in ph[i].colnames:
            result[i][j] = ph[i][j]
        for j in res.colnames:
            result[i][j] = res[j]

    return result


def process_light_curve(image_set, jd_key='JD', align_images=True, **kwargs):
    """Process photometry in different files and make lightcurve with them."""


def _find_pairs(x, y, match_pairs_tolerance):
    dx, dy = estimate_dxdy(x, y)
    pairs = match_pairs(x, y, dx, dy, tolerance=match_pairs_tolerance)

    o_idx = pairs['o']
    e_idx = pairs['e']

    tmp = Table()
    tmp['xo'] = x[o_idx]
    tmp['yo'] = y[o_idx]
    tmp['xe'] = x[e_idx]
    tmp['ye'] = y[e_idx]

    return tmp, pairs


def _do_polarimetry(phot_table, retarder_positions, retarder_type,
                    retarder_direction, pairs,
                    retarder_rotation=22.5):
    """Calculate the polarimetry of a given photometry table.

    phot_tables is a list of tables containing ['flux', 'flux_error']
    keys.
    """
    ph = phot_table

    if 'flux' not in ph[0].colnames or 'flux_error' not in ph[0].colnames:
        raise ValueError('Table for polarimetry must contain "flux" and'
                         ' "flux_error" keys.')

    if retarder_direction == 'cw':
        retarder_direction = -1
    if retarder_direction == 'ccw':
        retarder_direction = 1
    pos = np.array(retarder_positions)*retarder_rotation*retarder_direction

    tmp = Table()

    def _process(idx):
        f = 'flux'
        fe = 'flux_error'
        o = np.array([ph[j][f][pairs[idx]['o']] for j in range(len(ph))])
        e = np.array([ph[j][f][pairs[idx]['e']] for j in range(len(ph))])
        oe = np.array([ph[j][fe][pairs[idx]['o']] for j in range(len(ph))])
        ee = np.array([ph[j][fe][pairs[idx]['e']] for j in range(len(ph))])
        # res, err, therr = calculate_polarimetry(o, e, pos,
        #                                         retarder=retarder_type,
        #                                         o_err=oe, e_err=ee)
        for k in ['flux']:  # + list(res.keys()) + ['sigma_theor']:
            if k not in tmp.colnames:
                tmp.add_column(Column(name=k, dtype='f8',
                                      length=len(pairs)))
                if k != 'sigma_theor':
                    tmp.add_column(Column(name='{}_error'.format(k),
                                          dtype='f8',
                                          length=len(pairs)))
            if k == 'sigma_theor':
                # tmp[k][idx] = therr
                0 == 0
            elif k == 'flux':
                tmp[k][idx] = np.sum(o)+np.sum(e)
                tmp['{}_error'.format(k)][idx] = np.sum(oe)+np.sum(ee)
            else:
                # tmp[k][idx] = res[k]
                # tmp['{}_error'.format(k)][idx] = err[k]
                0 == 0

    for i in range(len(pairs)):
        _process(idx=i)

    return tmp


def process_polarimetry(image_set, align_images=True, retarder_type=None,
                        retarder_key=None, match_pairs_tolerance=1.0,
                        retarder_rotation=22.5, retarder_direction=None,
                        **kwargs):
    """Process the photometry and polarimetry of a set of images.

    kwargs are the arguments for the following functions:
    process_photometry, _solve_photometry
    """
    s = process_list(_check_ccddata, image_set)
    result = {'aperture': None, 'psf': None}

    sources = _do_aperture(np.sum([i.data for i in s], axis=0),
                           detect_fwhm=kwargs['detect_fwhm'],
                           detect_snr=kwargs['detect_snr'],
                           r=5)

    res_tmp, pairs = _find_pairs(sources['x'], sources['y'],
                                 match_pairs_tolerance=match_pairs_tolerance)
    try:
        wcs = _solve_astrometry(s[0].header, sources[pairs['o']],
                                s[0].data.shape,
                                ra_key=kwargs['ra_key'],
                                dec_key=kwargs['dec_key'],
                                plate_scale=kwargs['plate_scale'])
        idkwargs = {}
        for i in ['identify_catalog_file', 'identify_catalog_name', 'filter',
                  'identify_limit_angle', 'science_catalog',
                  'science_id_key', 'science_ra_key', 'science_dec_key']:
            if i in kwargs.keys():
                idkwargs[i] = kwargs[i]
        ids = _identify_star(Table([res_tmp['xo'], res_tmp['yo']],
                                   names=('x', 'y')), wcs,
                             **idkwargs)
        solve_phot = True
    except Exception as e:
        solve_phot = False
        ids = Table()
        logger.warn('Astrometry not solved. Ignoring photometry solving. {}'
                    .format(e))

    apkwargs = {}
    for i in ['photometry_type', 'detect_fwhm', 'detect_snr', 'box_size',
              'r_in', 'r_out', 'r_find_best', 'psf_model', 'psf_niters']:
        if i in kwargs.keys():
            apkwargs[i] = kwargs.get(i)

    ph = {}
    if check_iterable(kwargs['r']):
        for i in kwargs['r']:
            p = process_list(process_photometry, s, x=sources['x'],
                             y=sources['y'],
                             r=i, **apkwargs)
            ph[i] = {}
            ph[i]['aperture'] = [Table([j['aperture']['flux'],
                                        j['aperture']['flux_error']])
                                 if j['aperture'] is not None else None
                                 for j in p]
            ph[i]['psf'] = [Table([j['psf']['flux'], j['psf']['flux_error']])
                            if j['psf'] is not None else None
                            for j in p]
    else:
        p = process_list(process_photometry, s, x=sources['x'],
                         y=sources['y'],
                         r=kwargs['r'], **apkwargs)
        ph[kwargs['r']] = {}
        ph[kwargs['r']]['aperture'] = [Table([j['aperture']['flux'],
                                              j['aperture']['flux_error']])
                                       if j['aperture'] is not None else None
                                       for j in p]
        ph[kwargs['r']]['psf'] = [Table([j['psf']['flux'],
                                         j['psf']['flux_error']])
                                  if j['psf'] is not None else None
                                  for j in p]

    if solve_phot:
        solvekwargs = {}
        for i in ['montecarlo_iters', 'montecarlo_percentage',
                  'solve_photometry_type']:
            if i in kwargs.keys():
                solvekwargs[i] = kwargs.get(i)

    ret = [int(i.header[retarder_key], 16) for i in s]

    for ri in ph:
        for i in ['aperture', 'psf']:
            if result[i] is None and ph[ri][i][0] is not None:
                result[i] = Table()
                for m in ids.colnames:
                    result[i][m] = ids[m]

            if ph[ri][i][0] is not None:
                ap = _do_polarimetry([k for k in ph[ri][i]], ret,
                                     retarder_type=retarder_type,
                                     retarder_direction=retarder_direction,
                                     pairs=pairs,
                                     retarder_rotation=retarder_rotation)
                if solve_phot:
                    tmp = _solve_photometry(ap, cat_mag=ids['cat_mag'],
                                            **solvekwargs)
                    ap['mag'] = tmp['mag']
                    ap['mag_err'] = tmp['mag_err']

                for a in ap.colnames:
                    result[i]['r={}_{}'.format(ri, a)] = ap[a]

    return result


def run_pccdpack(image_set, retarder_type=None, retarder_key=None,
                 retarder_rotation=22.5, retarder_direction=None,
                 save_calib_path=None, r=np.arange(1, 21, 1),
                 r_in=60, r_out=70, gain_key=None, rdnoise_key=None,
                 **kwargs):
    files = []
    dtmp = mkdtemp(prefix='pccdpack')
    for i in range(len(image_set)):
        if isinstance(image_set[i], ccdproc.CCDData):
            name = os.path.join(dtmp, "image{:02d}.fits".format(i))
            image_set[i].to_hdu()[0].writeto(name)
            logger.debug("image {} saved to {}".format(i, name))
            files.append(name)
        elif isinstance(image_set[i], six.string_types):
            files.append(os.path.join(save_calib_path,
                                      os.path.basename(image_set[i])))
            try:
                shutil.copy(image_set[i], files[-1])
            except Exception:
                pass

    script = pccd.create_script(result_dir=dtmp, image_list=files,
                                star_name='object', apertures=r, r_ann=r_in,
                                r_dann=r_out-r_in,
                                gain_key=gain_key,
                                readnoise_key=rdnoise_key,
                                retarder=retarder_type,
                                auto_pol=True)
    print('\n\nExecute the following script:\n-----------------------------\n')
    print(script)
    time.sleep(0.5)
    print('------------------------------------\n')
    input('Press Enter when finished!')

    out_table = pccd.read_out(os.path.join(dtmp, 'object.out'),
                              os.path.join(dtmp, 'object.ord'))
    dat_table = Table.read(os.path.join(dtmp, 'dat.001'),
                           format='ascii.no_header')
    log_table = pccd.read_log(os.path.join(dtmp, 'object.log'))

    x, y = out_table['X0'], out_table['Y0']
    data = ccdproc.check_ccddata(files[0])
    ft = _do_aperture(data.data, x=x, y=y, r=5)
    astkwargs = {}
    for i in ['ra_key', 'dec_key', 'plate_scale']:
        if i in kwargs.keys():
            astkwargs[i] = kwargs[i]
    wcs = _solve_astrometry(data.header, ft, data.data.shape, **astkwargs)

    idkwargs = {}
    for i in ['identify_catalog_file', 'identify_catalog_name', 'filter',
              'identify_limit_angle', 'science_catalog',
              'science_id_key', 'science_ra_key', 'science_dec_key']:
        if i in kwargs.keys():
            idkwargs[i] = kwargs[i]
    ids = _identify_star(Table([x, y], names=('x', 'y')), wcs, **idkwargs)
    out_table['cat_id'] = ids['cat_id']
    out_table['sci_id'] = ids['sci_id']

    shutil.rmtree(dtmp)

    return out_table, dat_table, log_table


class ReduceScript():
    """Simple class to process pipeline scripts qith configuration files."""
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

    def process_product(self, filename, dataset='all'):
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
                batch_key_replace(prod)
                try:
                    self.run(i, **prod)
                except Exception as e:
                    logger.warn('Problem in the process of {} product from'
                                ' {} file. Passing it.\n'.format(i, filename) +
                                'Error: {}'.format(e))
                    if DEBUG:
                        os.system("echo \"{:<32} {:<32} {}\" >> {}"
                                  .format(filename, i, e, failed))

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
                  'mem_limit', 'calib_dir']:
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
        s = [os.path.join(config['raw_dir'], i) for i in config['sources']]

        calib_kwargs = {}
        for i in ('master_bias', 'master_flat', 'dark_frame', 'badpixmask',
                  'prebin', 'gain_key', 'gain', 'rdnoise_key',
                  'combine_method', 'combine_sigma', 'exposure_key',
                  'mem_limit', 'save_calib_path', 'combine_align_method',
                  'calib_dir', 'product_dir', 'remove_cosmics'):
            if i in config.keys():
                calib_kwargs[i] = config[i]
        ccd = calib_science(s, **calib_kwargs)

        # photkwargs = {}
        # for i in ['ra_key', 'dec_key', 'gain_key', 'rdnoise_key',
        #           'filter', 'plate_scale', 'photometry_type',
        #           'psf_model', 'r', 'r_in', 'r_out', 'psf_niters',
        #           'box_size', 'detect_fwhm', 'detect_snr', 'remove_cosmics',
        #           'align_images', 'solve_photometry_type',
        #           'montecarlo_iters', 'montecarlo_percentage',
        #           'identify_catalog_file', 'identify_catalog_name',
        #           'identify_limit_angle', 'science_catalog', 'science_id_key',
        #           'science_ra_key', 'science_dec_key']:
        #     if i in config.keys():
        #         photkwargs[i] = config[i]
        #
        # t = process_calib_photometry(ccd, **photkwargs)
        #
        # hdus = []
        # for i in t.keys():
        #     header_keys = ['solve_photometry_type', 'plate_scale']
        #     if i == 'aperture':
        #         header_keys += ['r', 'r_in', 'r_out', 'detect_fwhm',
        #                         'detect_snr']
        #     elif i == 'psf':
        #         header_keys += ['psf_model', 'box_size', 'psf_niters']
        #
        #     if config.get('solve_photometry_type', None) == 'montecarlo':
        #         header_keys += ['montecarlo_iters', 'montecarlo_percentage']
        #
        #     if config.get('identify_catalog_name', None) is not None:
        #         header_keys += ['identify_catalog_name',
        #                         'identify_limit_angle']
        #
        #     hdu = fits.BinTableHDU(t[i], name="{}_photometry".format(i))
        #     for k in header_keys:
        #         if k in config.keys():
        #             v = config[k]
        #             key = 'hierarch astrojc {}'.format(k)
        #             if check_iterable(v):
        #                 hdu.header[key] = ','.join([str(m) for m in v])
        #             else:
        #                 hdu.header[key] = v
        #     hdus.append(hdu)
        #
        # hdulist = fits.HDUList([*ccd.to_hdu(wcs_relax=True), *hdus])
        # hdulist.writeto(os.path.join(product_dir,
        #                              "photometry_{}".format(name)))


class PolarimetryScript(ReduceScript):
    def __init__(self, config=None):
        super(PolarimetryScript, self).__init__(config=config)

    def run(self, name, **config):
        """Run this pipeline script"""
        product_dir = config['product_dir']
        s = [os.path.join(config['raw_dir'], i) for i in config['sources']]

        calib_kwargs = {}
        for i in ('master_bias', 'master_flat', 'dark_frame', 'badpixmask',
                  'prebin', 'gain_key', 'gain', 'rdnoise_key',
                  'combine_method', 'combine_sigma', 'exposure_key',
                  'mem_limit', 'save_calib_path', 'combine_align_method',
                  'calib_dir', 'product_dir', 'remove_cosmics'):
            if i in config.keys():
                calib_kwargs[i] = config[i]
        ccds = calib_science(s, **calib_kwargs)

        # polkwargs = {}
        # for i in ['ra_key', 'dec_key', 'gain_key', 'rdnoise_key',
        #           'retarder_key', 'retarder_type', 'retarder_direction',
        #           'filter', 'plate_scale', 'photometry_type',
        #           'psf_model', 'r', 'r_in', 'r_out', 'psf_niters',
        #           'box_size', 'detect_fwhm', 'detect_snr', 'remove_cosmics',
        #           'align_images', 'solve_photometry_type',
        #           'match_pairs_tolerance', 'montecarlo_iters',
        #           'montecarlo_percentage', 'identify_catalog_file',
        #           'identify_catalog_name', 'identify_limit_angle',
        #           'science_catalog', 'science_id_key', 'science_ra_key',
        #           'science_dec_key']:
        #     if i in config.keys():
        #         polkwargs[i] = config[i]

        # t = process_polarimetry(ccds, **polkwargs)
        #
        # hdus = []
        # for i in t.keys():
        #     header_keys = ['retarder_type', 'retarder_rotation',
        #                    'retarder_direction', 'align_images',
        #                    'solve_photometry_type', 'plate_scale']
        #     if i == 'aperture':
        #         header_keys += ['r', 'r_in', 'r_out', 'detect_fwhm',
        #                         'detect_snr']
        #     elif i == 'psf':
        #         header_keys += ['psf_model', 'box_size', 'psf_niters']
        #
        #     if config.get('solve_photometry_type', None) == 'montecarlo':
        #         header_keys += ['montecarlo_iters', 'montecarlo_percentage']
        #
        #     if config.get('identify_catalog_name', None) is not None:
        #         header_keys += ['identify_catalog_name',
        #                         'identify_limit_angle']
        #
        #     hdu = fits.BinTableHDU(t[i], name="{}_polarimetry".format(i))
        #     for k in header_keys:
        #         if k in config.keys():
        #             v = config[k]
        #             key = 'hierarch astrojc {}'.format(k)
        #             if check_iterable(v):
        #                 hdu.header[key] = ','.join([str(m) for m in v])
        #             else:
        #                 hdu.header[key] = v
        #     hdus.append(hdu)
        #
        # image = ccdproc.combine(ccds, method='sum',
        #                         mem_limit=config.get('mem_limit', _mem_limit))
        #
        # hdulist = fits.HDUList([*image.to_hdu(wcs_relax=True), *hdus])
        # hdulist.writeto(os.path.join(product_dir,
        #                              "polarimetry_astrojc_{}".format(name)))

        # pccd = run_pccdpack(ccds, **polkwargs)
        # hdus = []
        # hdus.append(fits.BinTableHDU(pccd[0], name='out_table'))
        # hdus.append(fits.BinTableHDU(pccd[1], name='dat_table'))
        # hdulist = fits.HDUList([fits.PrimaryHDU(header=image.header),
        #                         *hdus])
        # hdulist.writeto(os.path.join(product_dir,
        #                              "polarimetry_pccdpack_{}".format(name)))


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
