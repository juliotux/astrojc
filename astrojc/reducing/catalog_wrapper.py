from astropy.coordinates.angles import Angle
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.io import ascii as asci
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy import units as u
import numpy as np
try:
    import bottleneck as bn
except ModuleNotFoundError:
    bn = np
import json
import six

from ..config import get_config_file
from ..py_utils import batch_key_replace


__all__ = ['Catalog', 'from_simbad', 'from_vizier', 'from_ascii',
           'solve_photometry_montecarlo', 'solve_photometry_median',
           'solve_photometry_average']


def _filter_replace(flux_key, flux_error_key, filter):
    d = {'flux_key': flux_key, 'flux_error_key': flux_error_key,
         'filter': filter}
    batch_key_replace(d)
    return d['flux_key'], d['flux_error_key']


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


def from_vizier(table, center, radius, filter, id_key='ID', ra_key='RAJ2000',
                dec_key='DEJ2000', flux_key='FLUX',
                flux_error_key=None, prepend_id_key=True, **ignoredargs):
    # TODO: implement the case for multiple tables
    v = Vizier()
    v.ROW_LIMIT = -1

    flux_key, flux_error_key = _filter_replace(flux_key, flux_error_key,
                                               filter)

    query = v.query_region(center, radius=Angle(radius), catalog=table)
    if len(query) == 0:
        raise RuntimeError("No online catalog results were found.")
    query = query[0]

    id = query[id_key].data
    if prepend_id_key:
        id = ["{id_key} {id}".format(id_key=id_key, id=i) for i in id]
        id = np.array(id)

    # FIXME: we assume, at now, that the units are hexa
    coords = SkyCoord(query[ra_key], query[dec_key],
                      unit=(u.hourangle, u.degree))

    try:
        flux = query[flux_key].data
    except Exception as e:
        flux = np.array([np.nan]*len(id))

    if flux_error_key is not None:
        try:
            flux_error = query[flux_error_key].data
        except Exception as e:
            flux_error = np.array([np.nan]*len(id))

    ra, dec = coords.ra.degree, coords.dec.degree

    return np.array(list(zip(id, ra, dec, flux, flux_error)),
                    np.dtype(list(zip(['id', 'ra', 'dec', 'flux',
                                       'flux_error'],
                                  ['S32']+['f8']*4))))


def from_ascii(filename, filter=None, id_key=None, ra_key=None, dec_key=None,
               flux_key=None, flux_error_key=None, prepend_id_key=False,
               **readkwargs):
    """Load from a ascii local catalog, with named columns."""
    t = asci.read(filename, **readkwargs)
    id = t[id_key].data
    ra = t[ra_key].data
    dec = t[dec_key].data

    try:
        ra = np.array([float(i) for i in ra])
        dec = np.array([float(i) for i in dec])
    except ValueError:
        coords = SkyCoord(ra, dec, unit=(u.hourangle, u.degree))
        ra = coords.ra.degree
        dec = coords.dec.degree

    flux_key, flux_error_key = _filter_replace(flux_key, flux_error_key,
                                               filter)

    try:
        flux = t[flux_key]
    except Exception as e:
        flux = [np.nan]*len(t)

    try:
        flux_error = t[flux_error_key]
    except Exception as e:
        flux_error = [np.nan]*len(t)

    if prepend_id_key:
        id = ['{id_key} {id}'.format(id_key=id_key, id=id)]

    return np.array(list(zip(id, ra, dec, flux, flux_error)),
                    np.dtype(list(zip(['id', 'ra', 'dec', 'flux',
                                       'flux_error'],
                                  ['S32']+['f8']*4))))


class Catalog():
    _mandatory = ['filters', 'prepend_id_key', 'id_key', 'ra_key',
                  'dec_key', 'flux_key', 'flux_unit', 'source']
    _vizier_mandatory = ['vizier_table']
    _simbad_mandatory = []
    _local_mandatory = ['file_name']
    _optional_keys = ['flux_error_key', 'comment', 'ads_ref_string',
                      '_readkwargs']
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
        self._defdict = {}

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
            self._defdict[i] = v

        batch_key_replace(self._defdict)

        if self['source'] == 'local':
            self._load_table()

    def __getitem__(self, key):
        if key not in self._defdict.keys():
            raise KeyError(key)
        return self._defdict[key]

    def __setitem__(self, key, value):
        self._defdict[key] = value

    @staticmethod
    def load_from_json(filename, key=None):
        j = json.load(open(filename, 'r'))
        if key is not None:
            j = j[key]
        return Catalog(j)

    @staticmethod
    def load_from_ascii(filename, id_key, ra_key, dec_key, flux_key,
                        flux_error_key, flux_unit, filters, prepend_id_key,
                        **readkwargs):
        new = Catalog({'source': 'local', 'file_name': filename,
                       'id_key': id_key, 'ra_key': ra_key, 'dec_key': dec_key,
                       'flux_key': flux_key, 'flux_error_key': flux_error_key,
                       'flux_unit': flux_unit, 'filters': filters,
                       'prepend_id_key': prepend_id_key,
                       '_readkwargs': readkwargs})
        return new

    def _load_table(self):
        self._table = from_ascii(self['file_name'], id_key=self['id_key'],
                                 ra_key=self['ra_key'],
                                 dec_key=self['dec_key'],
                                 flux_key=self['flux_key'],
                                 flux_error_key=self['flux_error_key'],
                                 **self['_readkwargs'])

    def _query(self, ra, dec, filter, radius):
        if self['source'] == 'simbad':
            return from_simbad(center=SkyCoord(ra, dec,
                                               unit=(u.degree, u.degree)),
                               filter=filter, radius=radius)
        else:
            if self['filters'] is not None:
                if filter not in self['filters']:
                    raise ValueError('Filter {} is not available.'
                                     .format(filter))

        if self['source'] == 'vizier':
            return from_vizier(center=SkyCoord(ra, dec,
                                               unit=(u.degree, u.degree)),
                               id_key=self['id_key'], ra_key=self['ra_key'],
                               dec_key=self['dec_key'],
                               flux_key=self['flux_key'],
                               flux_error_key=self['flux_error_key'],
                               filter=filter, radius=radius,
                               prepend_id_key=self['prepend_id_key'],
                               table=self['vizier_table'])
        elif self['source'] == 'local':
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
        radius = np.max([np.max(ra) - np.min(ra),
                         np.max(dec) - np.min(dec)])*u.degree

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


def solve_photometry_median(fluxes, flux_error, references, limits=(5, 18)):
    mags = -2.5*np.log10(fluxes)

    a, b = limits
    a, b = a, b if a < b else b, a
    args = np.where(np.logical_and(references >= a, references <= b))

    diff = references - mags
    dif = np.nanmedian(diff[args])
    err = np.nanstd(diff[args])

    error = 1.086*((flux_error + np.sqrt(fluxes))/fluxes) + err
    return mags + dif, error


def solve_photometry_average(fluxes, flux_error, references, limits=(5, 18)):
    mags = -2.5*np.log10(fluxes)

    a, b = limits
    a, b = a, b if a < b else b, a
    args = np.where(np.logical_and(references >= a, references <= b))

    diff = references - mags
    dif = np.nanaverage(diff[args], weights=np.divide(1, flux_error[args]))
    err = np.nanstd(diff[args])

    error = 1.086*((flux_error + np.sqrt(fluxes))/fluxes) + err
    return mags + dif, error


def _montecarlo_loop(args):
    mags = args[0]
    references = args[1]
    n_stars = args[2]

    iter_mags = np.zeros(len(mags))
    iter_mags[:] = np.nan

    choices = np.random.choice(len(mags), n_stars)
    iter_mags = mags + bn.nanmedian(references[choices] -
                                    mags[choices])
    return iter_mags


def solve_photometry_montecarlo(fluxes, flux_error, ref_mags, limits=(5, 18),
                                n_iter=100, n_stars=0.2):
    mags = -2.5*np.log10(fluxes)

    if float(n_stars).is_integer():
        n_stars = n_stars
    else:
        n_stars = max(1, int(n_stars*len(fluxes)))

    nrefs = np.array(ref_mags)
    lim = sorted(limits)
    filt = np.where(np.logical_or(nrefs < lim[0], nrefs > lim[1]))
    nrefs[filt] = np.nan

    args = (mags, nrefs, n_stars)
    iter_mags = [_montecarlo_loop(args) for i in range(n_iter)]

    result = np.nanmedian(iter_mags, axis=0)
    errors = np.nanstd(iter_mags, axis=0)
    return result, errors
