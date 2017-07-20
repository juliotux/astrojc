from astropy.coordinates.angles import Angle
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy import units as u
import functools
import bottleneck as bn
import numpy as np

def from_simbad(center, radius, filter):
    s = Simbad()
    s.add_votable_fields('fluxdata(%s)' % filter)

    query = s.query_region(center, radius)

    id = query['MAIN_ID'].data
    coords = SkyCoord(query['RA'], query['DEC'], unit=(u.hourangle, u.degree))
    try:
        flux = query['FLUX_%s' % filter].data
        flux_error = query['FLUX_ERROR_%s' % filter].data
    except:
        flux = np.array([np.nan]*len(id))
        flux_error = np.array([np.nan]*len(id))
    ra, dec = coords.ra.degree, coords.dec.degree

    return np.array(list(zip(id, ra, dec, flux, flux_error)),
                    np.dtype(list(zip(['id', 'ra', 'dec', 'flux', 'flux_error'],
                                 ['S32']+['f8']*4))))

def from_vizier(table, center, radius, idk='ID', rak='RAJ2000', deck='DEJ2000', fluxk='FLUX', fluxek='FLUX_ERROR',
                prepend_idk=True):
    v = Vizier()
    v.ROW_LIMIT = -1

    try:
        query = v.query_region(center, radius=Angle(radius), catalog=table)[0]
    except:
        raise AssertError('Catalog unable.')

    id = query[idk].data.astype('S32')
    if idk not in ['ID', 'MAIN_ID']:
        for i in range(len(id)):
            id[i] = "%s %s" % (idk, id[i])
    coords = SkyCoord(query[rak], query[deck], unit=(u.hourangle, u.degree))
    try:
        flux = query[fluxk].data
        flux_error = query[fluxek].data
    except:
        flux = np.array([np.nan]*len(id))
        flux_error = np.array([np.nan]*len(id))
    ra, dec = coords.ra.degree, coords.dec.degree

    return np.array(list(zip(id, ra, dec, flux, flux_error)),
                    np.dtype(list(zip(['id', 'ra', 'dec', 'flux', 'flux_error'],
                                 ['S32']+['f8']*4))))

def from_bon10(center, radius, **ks):
    f1 = functools.partial(from_vizier, 'J/AJ/140/416/table3', idk='Name', prepend_idk=False)
    f2 = functools.partial(from_vizier, 'J/AJ/138/1003/table3', idk='Name', prepend_idk=False)
    try:
        return np.append(f1(center, radius, **ks), f2(center, radius, **ks))
    except:
        try:
            return f1(center, radius, **ks)
        except:
            return f2(center, radius, **ks)

from_ucac4 = functools.partial(from_vizier, 'I/322A', idk='UCAC4')
from_denis = functools.partial(from_vizier, 'B/denis', idk='DENIS')
from_2mass = functools.partial(from_vizier, 'II/246', idk='_2MASS')

def match_catalog(sources, catalog):
    '''
    everything in degrees:
    sources = (ra, dec)
    catalog = (ra, dec)
    '''

    s_cat = SkyCoord(ra=sources[0], dec=sources[1], frame='icrs', unit=('degree', 'degree'))
    r_cat = SkyCoord(ra=catalog[0], dec=catalog[1], frame='icrs', unit=('degree', 'degree'))

    idx, sep2d, sep3d = match_coordinates_sky(s_cat, r_cat, 1)
    return idx, sep2d

def median_comparision(fluxes, flux_error, references, limits=(5,18)):
    mags = -2.5*np.log10(fluxes)
    for i in range(len(references)):
        if not limits[0] <= references[i] <= limits[1]:
            references[i] = np.nan
    diff = references - mags
    dif = bn.nanmedian(diff)
    err = bn.nanstd(diff)

    error = 1.086*((flux_error + np.sqrt(fluxes))/fluxes) + err
    return mags + dif, error

def montecarlo_comparision(fluxes, flux_error, references, limits=(5, 18), n_iter=100, n_stars=0.5):
    mags = -2.5*np.log10(fluxes)

    if float(n_stars).is_integer():
        n_stars = n_stars
    else:
        n_stars = int(n_stars*len(fluxes))

    iter_mags = np.zeros((len(fluxes), n_iter), dtype='f8')
    for i in range(n_iter):
        for j in range(len(fluxes)):
            choices = np.random.choice(len(fluxes), n_stars)
            iter_mags[j,i] = mags[j] + bn.nanmedian(references[choices] - mags[choices])

    result = bn.nanmedian(iter_mags, axis=1)
    errors = bn.nanstd(iter_mags, axis=1)
    return result, errors
