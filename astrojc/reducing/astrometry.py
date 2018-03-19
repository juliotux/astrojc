import six
import numpy as np
from astropy.wcs import WCS
from astropy.table import Table
from .astrometry_net_wrapper import (solve_astrometry_xy, wcs_xy2radec,
                                     guess_coordinates)
from .catalog_wrapper import Catalog
from ..py_utils import process_list, string_fix


__all__ = ['solve_astrometry', 'identify_stars', 'wcs_from_coords']


_angles = {
    'left': 180,
    'right': 0,
    'top': 90,
    'bottom': 270
}


def wcs_from_coords(x, y, ra, dec, plate_scale, north, flip=None):
    """Giving coordinates and plate scale, creates a WCS.
    x, y: float, pixel coordinates in image
    ra, dec: float, sky coordinates

    plate_scale in arcsec/pix

    north direction can be angles ('cw', from top axis) or
    ['right', 'top', 'left', 'bottom'] to angles [270, 0, 90, 180].

    flip means if some direction is mirrored in image.
    Can be 'ra', 'dec', 'all' or None
    The standard coordinates are ra and dec in the following order, ccw:
    W - N - E - S

    Problems:
    This algorith is not good for coordinates far from crpix. But is useful
    when you cannot solve with other algorithms. (Like just one star in field).
    """
    sk = guess_coordinates(ra, dec)
    ra, dec = sk.ra.degree, sk.dec.degree

    if isinstance(north, six.string_types):
        if north in _angles.keys():
            north = _angles[north]
        else:
            raise ValueError('invalid value for north: {}'.format(north))

    # convert arcsec/pix to degree/pix
    plate_scale /= 3600

    # following the standard astrometry.net, all the delta informations
    # will be located inside the rotation/scale matrix pc
    deltra = -plate_scale if flip in ['ra', 'all'] else plate_scale
    deltde = -plate_scale if flip in ['dec', 'all'] else plate_scale

    rot = [[np.sin(np.radians(north)), -np.cos(np.radians(north))],
           [np.cos(np.radians(north)), np.sin(np.radians(north))]]

    pc = np.multiply(rot, [[deltra]*2, [deltde]*2])

    # and finally, we can construct the wcs
    w = WCS(naxis=2)
    w.wcs.crpix = [x, y]
    w.wcs.crval = [ra, dec]
    w.wcs.cdelt = [1, 1]
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.cunit = ['deg', 'deg']
    w.wcs.pc = pc
    # FIXME: Check why astrometry.net put pole in the crpix
    w.wcs.latpole = dec
    w.wcs.lonpole = 180

    return w


def solve_astrometry(table, header, shape, ra_key=None, dec_key=None,
                     plate_scale=None):
    """Solves the astrometry of a field and return a valid wcs."""
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
        flux = table['flux']
        wcs = solve_astrometry_xy(x, y, flux, header, imw, imh,
                                  image_params=im_params, return_wcs=True)
    return wcs


def identify_stars(table, wcs, filter, identify_catalog_file,
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
        sci_names, _, _ = sci.query_id_mag(ra, dec, None,
                                           limit_angle=limit_angle)
        res['sci_id'] = process_list(string_fix, sci_names)

    res['cat_id'] = process_list(string_fix, name)
    res['ra'] = ra
    res['dec'] = dec
    res['cat_mag'] = mag
    res['cat_mag_err'] = mag_err

    return res
