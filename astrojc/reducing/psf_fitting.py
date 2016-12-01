'''Standalone psf fitting functions. May be problematic.'''

# Fitting Functions
import numpy as np
from scipy.optimize import curve_fit
from scipy import integrate
from astropy.nddata.utils import extract_array
from astropy.stats import sigma_clipped_stats

#import pyximport; pyximport.install()
import psf_kernels
from numba import autojit

@autojit
def compute_sky(z, sigma=2, mode='mean'):
    #TODO: Mode can be 'plane' too, but need to be implemented.
    '''
    mode:
        mean: compute de mean of 33% lower values
        sigma_clip: compute the sigma_clipped stats and do the median of the
                    values between the lower value and n*sigma.
    '''
    if mode == 'mean':
        z = np.ravel(z)
        return np.mean(z[np.argsort(z)[:int(len(z)/3)]])
    elif mode == 'sigma_clip':
        mean, median, rms = sigma_clipped_stats(z)

        newz = np.ravel(z)
        return np.nanmedian(newz[newz < np.min(z) + sigma*rms])
    else:
        raise ValueError('Sky compute mode %s unrecognized.' % str(mode))

@autojit
def xy2r(x, y, data, xc, yc):
    r = np.sqrt((x-xc)**2 + (y-yc)**2)
    return np.ravel(r), np.ravel(data)

@autojit
def extract_data(data, indices, box_size, position):
    x, y = position
    dx = dy = float(box_size)/2

    x_min = max(int(x-dx), 0)
    x_max = int(x+dx)+1
    y_min = max(int(y-dy), 0)
    y_max = int(y+dy)+1

    d = data[y_min:y_max, x_min:x_max]
    xi = indices[1][y_min:y_max, x_min:x_max]
    yi = indices[0][y_min:y_max, x_min:x_max]
    return d, xi, yi

def generic_radial_fit(data, positions, box_size,
                       dtype, function, flux_function, nparams,
                       precalc_sky=True, sky_method='mean',
                       compute_errors=True):

    results = np.zeros(len(positions), dtype)
    if compute_errors:
        errors = np.zeros(len(positions), dtype)

    indices = np.indices(data.shape)
    for i in range(len(positions)):
        xp, yp = positions[i]
        d, xi, yi = extract_data(data, indices, box_size, (xp, yp))
        r, f = xy2r(xi, yi, d, xp, yp)

        if precalc_sky:
            sky = compute_sky(d, sky_method)
        else:
            sky = 0.0

        try:
            guess = tuple([1]*(nparams-1) + [sky])
            params, p_errors = curve_fit(function, r, f, p0=guess)
            p_errors = tuple([j for j in np.diag(p_errors)])
        except:
            nantuple = tuple([np.nan]*nparams)
            params, p_errors = nantuple, nantuple

        flux, flux_error = flux_function(params, errors=p_errors)
        r = tuple([xp, yp, flux, flux_error] + list(params) + list(p_errors))
        results[i] = np.array(r, dtype=dtype)

    return results

def generic_spatial_fit(data, positions, box_size,
                        dtype, function, flux_function, nparams,
                        precalc_sky=True, sky_method='mean',
                        compute_errors=True):
    results = np.zeros(len(positions), dtype)
    if compute_errors:
        errors = np.zeros(len(positions), dtype)

    indices = np.indices(data.shape)
    for i in range(len(positions)):
        xp, yp = positions[i]
        d, xi, yi = extract_data(data, indices, box_size, (xp, yp))

        if precalc_sky:
            sky = compute_sky(d, sky_method)
        else:
            sky = 0.0

        try:
            guess = tuple([xp, yp] + [1]*(nparams-3) + [sky])
            params, p_errors = curve_fit(function, (xi, yi), np.ravel(d), p0=guess)
            p_errors = tuple([j for j in np.diag(p_errors)])
        except:
            nantuple = tuple([np.nan]*nparams)
            params, p_errors = nantuple, nantuple

        (xp, yp), (xp_err, yp_err) = params[0:2], p_errors[0:2]
        params, p_errors = params[2:] , p_errors[2:]

        flux, flux_error = flux_function(params, errors=p_errors)
        r = tuple([xp, yp, flux, flux_error] + list(params) + list(p_errors))
        results[i] = np.array(r, dtype=dtype)

    return results

#Gaussian functions#################################

def gaussian_radial(x, amplitude, sigma, sky):
    return psf_kernels.gaussian_radial(x, sigma, amplitude, sky)

def gaussian_spatial(xy, x_0, y_0, amplitude, sigma_x, sigma_y, theta, sky):
    x = xy[0]
    y = xy[1]
    return np.array(psf_kernels.gaussian_spatial(np.ravel(x), np.ravel(y), x_0, y_0, sigma_x, sigma_y, theta, amplitude)) + sky

def gaussian_flux_compute(params, errors=None):
    try:
        if len(params) == 5:
            flux = 2*np.pi*params[0]*np.abs(params[1]*params[2])
        elif len(params) == 3:
            flux = 2*np.pi*params[0]*np.abs(params[1]*params[1])
        else:
            raise ValueError('The imput params doesn\'t correspond to a gaussian fit.')
    except:
        flux = np.nan

    try:
        if errors is not None:
            if len(params) == 3:
                flux_error = flux*(errors[0]/params[0] + 2*errors[1]/params[1])
            elif len(params) == 5:
                flux_error = flux*(errors[0]/params[0] + errors[1]/params[1] + errors[2]/params[2])
        else:
            flux_error = 0.0
    except:
        flux_error = np.nan

    return flux, flux_error

def fit_gaussian_radial(data, positions, box_size,
                        precalc_sky=True, sky_method='mean',
                        compute_errors=True):
    dtype = np.dtype(zip(['x','y','flux','flux_error',
                          'amplitude','sigma','sky',
                          'amplitude_err', 'sigma_err', 'sky_err'],
                         ['f8']*10))

    return generic_radial_fit(data, positions, box_size,
                              dtype, gaussian_radial, gaussian_flux_compute, 3,
                              precalc_sky=precalc_sky, sky_method=sky_method,
                              compute_errors=compute_errors)

def fit_gaussian_spatial(data, positions, box_size,
                         precalc_sky=True, sky_method='mean',
                         compute_errors=True):
    dtype = np.dtype(zip(['x','y','flux','flux_error',
                          'amplitude','sigmax','sigmay','theta','sky',
                          'amplitude_err', 'sigmax_err','sigmay_err','theta_err', 'sky_err'],
                         ['f8']*14))

    return generic_spatial_fit(data, positions, box_size,
                               dtype, gaussian_spatial, gaussian_flux_compute, 7,
                               precalc_sky=precalc_sky, sky_method=sky_method,
                               compute_errors=compute_errors)

#Moffat functions#################################

def moffat_radial(r, amplitude, gamma, alpha, sky):
    return psf_kernels.moffat_radial(r, gamma, alpha, amplitude, sky)

def moffat_spatial(xy, x0, y0, amplitude, gamma, alpha, sky):
    x = xy[0]
    y = xy[1]
    return np.array(psf_kernels.moffat_spatial(np.ravel(x), np.ravel(y), x0, y0, gamma, alpha, amplitude)).ravel() + sky

def moffat_flux_compute(params, errors=None):
    try:
        r_max = params[1]*np.sqrt(10**(4/params[2]) - 1) #I(r_max) = 10^(-4)I(0)
        flux, _ = integrate.nquad(psf_kernels.moffat_integrate, [[0, r_max],[0, 2*np.pi]], args=params)
        flux_error = np.nan
        return flux, flux_error
    except:
        return np.nan, np.nan

def fit_moffat_radial(data, positions, box_size,
                       precalc_sky=True, sky_method='mean',
                       compute_errors=True):
    dtype = np.dtype(zip(['x','y','flux','flux_error',
                          'amplitude','gamma','alpha','sky',
                          'amplitude_err','gamma_err','alpha_err','sky_err'],
                         ['f8']*12))

    return generic_radial_fit(data, positions, box_size,
                              dtype, moffat_radial, moffat_flux_compute, 4,
                              precalc_sky=precalc_sky, sky_method=sky_method,
                              compute_errors=compute_errors)

def fit_moffat_spatial(data, positions, box_size,
                       precalc_sky=True, sky_method='mean',
                       compute_errors=True):
    dtype = np.dtype(zip(['x','y','flux','flux_error',
                          'amplitude','gamma','alpha','sky',
                          'amplitude_err','gamma_err','alpha_err','sky_err'],
                         ['f8']*12))

    return generic_spatial_fit(data, positions, box_size,
                               dtype, moffat_spatial, moffat_flux_compute, 6,
                               precalc_sky=precalc_sky, sky_method=sky_method,
                               compute_errors=compute_errors)
