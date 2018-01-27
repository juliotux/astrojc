from astropy.coordinates.angles import Angle
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.io import ascii as asci
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy import units as u
import functools
import bottleneck as bn
import numpy as np
from scipy.spatial.ckdtree import cKDTree
import json
import os

from .astrometry_wrapper import wcs_xy2radec as xy2radec
from ..logging import log
from ..config import get_config_file


__all__ = ['Catalog', 'from_simbad', 'from_vizier', 'from_table',
           'solve_photometry_montecarlo', 'solve_photometry_median',
           'solve_photometry_average']


def from_simbad(center, radius, filter):
    s = Simbad()
    s.add_votable_fields('fluxdata({filter})'.format(filter=filter))

    query = s.query_region(center, radius)

    id = query['MAIN_ID'].data
    coords = SkyCoord(query['RA'], query['DEC'], unit=(u.hourangle, u.degree))
    try:
        flux = query['FLUX_{filter}'.format(filter=filter)].data
        flux_error = query['FLUX_ERROR_{filter}'.format(filter=filter)].data
    except NameError:
        flux = np.array([np.nan]*len(id))
        flux_error = np.array([np.nan]*len(id))
    ra, dec = coords.ra.degree, coords.dec.degree

    return np.array(list(zip(id, ra, dec, flux, flux_error)),
                    np.dtype(list(zip(['id', 'ra', 'dec', 'flux',
                                       'flux_error'],
                                  ['S32']+['f8']*4))))


def from_vizier(table, center, radius, id_key='ID', ra_key='RAJ2000',
                dec_key='DEJ2000', flux_key='FLUX',
                flux_error_key=None, prepend_id_key=True):
    v = Vizier()
    v.ROW_LIMIT = -1

    query = v.query_region(center, radius=Angle(radius), catalog=table)[0]

    id = query[id_key].data.astype('S32')
    if prepend_id_key:
        for i in range(len(id)):
            id[i] = "{id_key} {id}".format(id_key=id_key, id=id[i])

    # FIXME: we assume, at now, that the units are hexa
    coords = SkyCoord(query[ra_key], query[dec_key],
                      unit=(u.hourangle, u.degree))

    try:
        flux = query[flux_key].data
    except NameError:
        flux = np.array([np.nan]*len(id))

    if flux_error_key is not None:
        try:
            flux_error = query[flux_error_key].data
        except NameError:
            flux_error = np.array([np.nan]*len(id))

    ra, dec = coords.ra.degree, coords.dec.degree

    return np.array(list(zip(id, ra, dec, flux, flux_error)),
                    np.dtype(list(zip(['id', 'ra', 'dec', 'flux',
                                       'flux_error'],
                                  ['S32']+['f8']*4))))


def from_ascii(filename, id_key=None, ra_key=None, dec_key=None, flux_key=None,
               flux_error_key=None, prepend_id_key=False, **readkwargs):
    """Load from a ascii local catalog, with named columns."""
    t = asci.read(filename, **asciikwargs)
    id = t[id_key]
    ra = t[ra_key]
    dec = t[dec_key]

    try:
        flux = t[flux_key]
    except KeyError:
        flux = [np.nan]*len(t)

    try:
        flux_error = t[flux_error_key]
    except KeyError:
        flux_error = [np.nan]*len(t)

    if prepend_id_key:
        id = ['{id_key} {id}'.format(id_key=id_key, id=id)]

    return np.array(list(zip(id, ra, dec, flux, flux_error)),
                    np.dtype(list(zip(['id', 'ra', 'dec', 'flux',
                                       'flux_error'],
                                  ['S32']+['f8']*4))))


class Catalog():
    _mandatory = ['filters', 'prepend_id_key', 'id_key', 'ra_key',
                  'dec_key', 'flux_key', 'flux_unit']
    _vizier_mandatory = ['vizier_table']
    _simbad_mandatory = []
    _local_mandatory = ['file_name']
    _optional_keys = ['flux_error_key', 'comment', 'ads_ref_string']
    """Container for store the definitions of a catalog and related routines"""
    def __init__(self, definitions):
        """
        definitions : dict_like
            A dictionary containing the needed informations for the catalog.
            The mandatory keys are:
            - filters : list of available filters for the catalog
                like ["B", "V", "g", "r", "i"]
            - source : the catalog source
                can be "local", "vizier", "simbad"
            - vizier_table : the table id in vizier_table (if source=="vizier")
                like "I/322A"
            - file_name : the local file name (if source=="local")
                like "/home/user/mytable.dat"
            - prepend_id_key : True if you want to prepend the id_key to
            the name of the star in the output id in query. Default is False.
            - id_key : key for the column containing the ID or name of the
            star.
                like "UCAC4"
            - ra_key : key for the column containing the RA of the star.
                like "RAJ2000"
            - dec_key : key for the column containing the DEC of the star.
                like "DEJ2000"
            - flux_key : key for the column containing the flux of the star.
            The string {filter} will be replaced by the filter name passed to
            query functions.
                like "{filter}mag",
            - flux_error_key : key for the column containing the flux error
            of the star. The string {filter} will be replaced by the filter
            name passed to query functions. Optional.
                like "e_{filter}mag",
            - flux_unit : key for the column containing the flux unit of the
            star.
                can be ["mag", "flux"]
            - ads_ref_string : reference string for ADS. Optional.
                like "2013AJ....145...44Z"
            - comment : Aditional comment. Optional.

        The catalog definitions can be loaded from a json file using
        Catalog.load_from_json(filename) function.
        """
        def _check_key(key, list_of_keys):
            if key not in list_of_keys:
                raise ValueError("{} mandatory key not found in"
                                 " definitions.".format(key))

        self._table = None

        _check_key('source', definitions.keys())

        source = definitions['source']
        if source != 'simbad':
            # simbad catalog do not need all other mandatory keys
            for i in self._mandatory:
                _check_key(i, definitions.keys())

        if source == 'vizier':
            for i in self._vizier_mandatory:
                _check_key(i, definitions.keys())
        elif source == 'simbad':
            for i in self._simbad_mandatory:
                definitions['flux_unit'] = 'mag'
                _check_key(i, definitions.keys())
        elif source == 'local':
            for i in self._local_mandatory:
                _check_key(i, definitions.keys())
        else:
            raise ValueError("source {} not supported.".format(source))

        for i, v in definitions.items():
            if i in list(self._mandatory + self._simbad_mandatory +
                         self._vizier_mandatory + self._local_mandatory +
                         self._optional_keys):
                self.__setattr__(i, v)

        if self.souce == 'local':
            self._load_table()

    @staticmethod
    def load_from_json(filename, key=None):
        j = json.load(open(filename, 'r'))
        if key is not None:
            j = j[key]
        return Catalog(j)

    @staticmethod
    def load_from_ascii(filename, id_key, ra_key, dec_key, flux_key,
                        flux_error_key, flux_unit, **readkwargs):
        new = Catalog({'souce': 'local', 'file_name': filename,
                       'ra_key': ra_key, 'dec_key': dec_key,
                       'flux_key': flux_key, 'flux_error_key': flux_error_key,
                       'flux_unit': flux_unit})

    def _load_table(self):
        self._table = from_ascii(filename, id_key=self.id_key,
                                 ra_key=self.ra_key, dec_key=self.dec_key,
                                 flux_key=self.flux_key,
                                 flux_error_key=self.flux_error_key,
                                 **readkwargs)

    def _query(self, ra, dec, filter, radius):
        if self.source == 'simbad':
            return from_simbad(center=SkyCoord(ra, dec,
                                               unit=(u.degree, u.degree)),
                               filter=filter, radius=radius)
        else:
            if filter not in self.filters:
                raise ValueError('Filter {} is not available.'.format(filter))

        if self.source == 'vizier':
            kwargs = dict()
            for i in self._mandatory - ['filters', 'flux_unit']:
                kwargs[i] = self.__dict__[i]
            return from_vizier(center=SkyCoord(ra, dec,
                                               unit=(u.degree, u.degree)),
                               id_key=self.id_key, ra_key=self.ra_key,
                               dec_key=self.dec_key, flux_key=self.flux_key,
                               flux_error_key=self.flux_error_key,
                               filter=filter, radius=radius, **kwargs)
        elif self.source == 'local':
            if self._table is None:
                raise ValueError('Table not loaded!')
            return self._table
        else:
            raise ValueError('Source not available. This should not happen!')

    def query_id_mag(self, ra, dec, filter, limit_angle='2 arcsec'):
        """Query the photometry parameters of a list of RA and DEC of objects
        in this catalog.

        Return: id, mag, mag_err (or fluxes)
        """
        center_ra = (np.max(ra) + np.min(ra))/2
        center_dec = (np.max(dec) + np.min(dec))/2
        radius = np.max(np.max(ra) - np.min(ra),
                        np.max(dec) - np.min(dec))*u.degree

        cat = self._query(center_ra, center_dec, filter, radius)

        i, d, _ = match_coordinates_sky(SkyCoord(ra, dec, unit=('degree',
                                                                'degree'),
                                                 frame='icrs'),
                                        SkyCoord(cat['ra'], cat['dec'],
                                                 unit=('degree', 'degree'),
                                                 frame='icrs'))

        id = np.zeros(len(ra), dtype='S32')
        mag = np.zeros(len(ra), dtype='f8')
        mag_err = np.zeros(len(ra), dtype='f8')
        mag.fill(np.nan)
        mag_err.fill(np.nan)

        lim = Angle(limit_angle)
        for k in range(len(ra)):
            if d[k] <= lim:
                id[k] = cat['id'][i[k]]
                mag[k] = cat['flux'][i[k]]
                mag_err[k] = cat['flux_error'][i[k]]

        return id, mag, mag_err


class PhotometrySolver():
    "Solve the photometry of a field by median or montecarlo comparisions."
    # TODO: Unfinished. Think better on it.
    def __init__(self, filter, catalog, catalog_config_file=None):
        """catalog can be a string, from predefineds catalog in
        ~/config/.astrojc/photometry_catalogs.json, or a Catalog object."""
        if isinstance(catalog, Catalog):
            self._catalog = catalog
        elif isinstance(catalog, string):
            try:
                if catalog_config_file is not None:
                    f = get_config_file('photometry_catalogs.json')
                else:
                    f = catalog_config_file
                self._catalog = Catalog.load_from_json(f)
            except FileNotFoundError:
                raise('File {} not found. Please check its existence to '
                      'use predefined catalogs.'.format(f))
            except Exception as e:
                raise('This file is not a valid config catalog.'
                      ' Due to: {}'.format(e))

        self._filter = filter

        self._operator = None
        self._inverse_operator = None

    def _to_cat_scale(self, data, data_error, data_scale):
        """Transform the user data to the same unit of catalog data."""
        target = self._catalog.flux_unit
        if target == 'mag':
            self._operator = lambda x, y: x - y
            self._inverse_operator = lambda x, y: x + y
            if data_unit in ['linear']:
                if data_error is None:
                    err = None
                else:
                    err = 10.086*np.divide(flux_error, flux)
                return -2.5*np.log10(data), err
            elif data_unit in ['log']:
                if data_error is None:
                    err = None
                else:
                    err = np.multiply(2.5, data_error)
                return np.multiply(-2.5, data), err
            elif data_unit in ['mag']:
                return data, data_error
        elif target == 'log':
            self._operator = lambda x, y: x - y
            self._inverse_operator = lambda x, y: x + y
            if data_unit in ['linear', 'adu', 'count', 'flux']:
                if data_error is None:
                    err = None
                else:
                    err = 0.4342*np.divide(flux_error, flux)
                return np.log10(data), err
            elif data_unit in ['log10', 'log']:
                return data, data_error
            elif data_unit in ['mag']:
                if data_error is None:
                    err = None
                else:
                    err = np.divide(flux_error, 2.5)
                return np.divide(data, -2.5), err
        elif target == 'linear':
            self._operator = lambda x, y: x / y
            self._inverse_operator = lambda x, y: x * y
            if data_unit in ['linear', 'adu', 'count', 'flux']:
                return data, data_error
            elif data_unit in ['log10', 'log']:
                if data_error is None:
                    err = None
                else:
                    err = 0.4342 * np.power(10, data_error)
                return 10**data, err
            elif data_unit in ['mag']:
                if data_error is None:
                    err = None
                else:
                    err = 0.921*np.exp(-0.921*data_error)
                return np.power(10, -0.4*data)
        else:
            raise NotImplementedError

    def solve_median(self, ra, dec, flux, flux_error=None, mag_limits=(5, 18),
                     flux_unit='linear'):
        """Solve the photometry of a dataset by median comparision."""

    def solve_montecarlo(fluxes, flux_error, references, limits=(5, 18),
                         n_iter=100, n_stars=0.2):
        """Solve the photometry of a dataset by montecarlo comparision.""""


def solve_photometry_median(fluxes, flux_error, references, limits=(5, 18)):
    mags = -2.5*np.log10(fluxes)

    a, b = limits
    a, b = a, b if a < b else b, a
    args = np.where(np.logical_and(references>=a, references<=b))

    diff = references - mags
    dif = np.nanmedian(diff[args])
    err = np.nanstd(diff[args])

    error = 1.086*((flux_error + np.sqrt(fluxes))/fluxes) + err
    return mags + dif, error


def solve_photometry_average(fluxes, flux_error, references, limits=(5, 18)):
    mags = -2.5*np.log10(fluxes)

    a, b = limits
    a, b = a, b if a < b else b, a
    args = np.where(np.logical_and(references>=a, references<=b))

    diff = references - mags
    dif = np.nanaverage(diff[args], weights=np.divide(1, flux_error[args]))
    err = np.nanstd(diff[args])

    error = 1.086*((flux_error + np.sqrt(fluxes))/fluxes) + err
    return mags + dif, error


def solve_photometry_montecarlo(fluxes, flux_error, references, limits=(5, 18),
                                n_iter=100, n_stars=0.2):
    mags = -2.5*np.log10(fluxes)

    if float(n_stars).is_integer():
        n_stars = n_stars
    else:
        n_stars = int(n_stars*len(fluxes))

    iter_mags = np.zeros((len(fluxes), n_iter), dtype='f8')
    for i in range(n_iter):
        for j in range(len(fluxes)):
            choices = np.random.choice(len(fluxes), n_stars)
            iter_mags[j, i] = mags[j] + bn.nanmedian(references[choices] -
                                                     mags[choices])

    result = np.nanmedian(iter_mags, axis=1)
    errors = np.nanstd(iter_mags, axis=1)
    return result, errors
