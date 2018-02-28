import numpy as np
from scipy.spatial import cKDTree
from astropy.modeling.fitting import LevMarLSQFitter, LinearLSQFitter
from astropy.modeling.functional_models import Const1D

from ..logging import log as logger
from ..math.polarimetry_models import HalfWaveModel, QuarterWaveModel


def estimate_dxdy(x, y, steps=[100, 30, 5, 3], bins=30):
    """Estimate the x and y distance between the pairs of points"""
    dxa = []
    dya = []

    def _find_max(d, steps=[100, 30, 5, 3], bins=30):
        dx = 0
        for lim in (np.max(d), *steps):
            lo, hi = (dx-lim, dx+lim)
            lo, hi = (lo, hi) if (lo < hi) else (hi, lo)
            histx = np.histogram(d, bins=bins, range=[lo, hi])
            mx = np.argmax(histx[0])
            dx = (histx[1][mx]+histx[1][mx+1])/2
        return dx

    for i in range(len(x)):
        for j in range(len(x)):
            if x[i] < x[j]:
                dya.append(y[i] - y[j])
                dxa.append(x[i] - x[j])

    return (_find_max(dxa, steps=steps, bins=bins),
            _find_max(dya, steps=steps, bins=bins))


def match_pairs(x, y, dx, dy, tolerance=1.0):
    """Match the pairs of ordinary/extraordinary points (x, y)."""
    dt = np.dtype([('o', int), ('e', int)])
    results = np.zeros(len(x), dtype=dt)
    npairs = 0

    p = list(zip(x, y))
    kd = cKDTree(p)

    for i in range(len(p)):
        px = p[i][0]-dx
        py = p[i][1]-dy
        d, j = kd.query((px, py), k=1, eps=tolerance,
                        distance_upper_bound=tolerance, n_jobs=-1)
        if d <= tolerance:
            results[npairs]['o'] = i
            results[npairs]['e'] = j
            npairs = npairs+1
            kd = cKDTree(p)

    return results[:npairs]


def __trashed_estimate_normalize(o, e, positions, n_consecutive):
    """Estimate the normalization of a given set of data.
    Trashed, do not use!
    """
    data_o = [[]]*n_consecutive
    data_e = [[]]*n_consecutive

    # First, we separate the data in the positions, relative to consecutive
    for i, oi, ei in zip(positions, o, e):
        index = int(i/n_consecutive)
        data_o[index].append(oi)
        data_e[index].append(ei)

    # check if all positions have a value
    for i in data_o:
        if i == []:
            logger.warn('Could not calculate polarimetry normalization. '
                        'Not all needed positions are available. Using k=1.')
            return 1

    # Now we use as each consecutive value the mean of the values in each index
    for i in range(n_consecutive):
        data_o[index] = np.nanmean(data_o[index])
        data_e[index] = np.nanmean(data_e[index])

    # Now, assuming the k will multiply e
    k = np.sum(data_o)/np.sum(data_e)
    logger.debug('Polarimetry normalization estimated as k={}'.format(k))
    return k


def estimate_normalize(o, e, positions):
    """Estimate the normalization of a given set of data."""
    # We can estimate the k fitting a constant function to the data!
    # y(phi) = k = e(phi)/o(phi)
    fitter = LinearLSQFitter()
    fit = fitter(Const1D(), positions, e/o)

    return fit.parameters[0]


def calculate_polarimetry(o, e, positions, rotation_interval,
                          retarder='half', o_err=None, e_err=None,
                          normalize=True):
    """Calculate the polarimetry."""

    if retarder == 'half':
        model = HalfWaveModel()
    elif retarder == 'quarter':
        model = QuarterWaveModel()
    else:
        raise ValueError('retarder {} not supported.'.format(retarder))

    o = np.array(o)
    e = np.array(e)
    positions = np.array(positions)

    if normalize:
        k = estimate_normalize(o, e, positions)
        z = (o-(e*k))/(o+(e*k))
    else:
        z = (o-e)/(o+e)
    psi = positions*rotation_interval
    if o_err is None or e_err is None:
        z_erro = None
        th_error = None
    else:
        # Assuming individual z errors from propagation
        o_err = np.array(o_err)
        e_err = np.array(e_err)
        oi = 2*o/((o+e)**2)
        ei = -2*e/((o+e)**2)
        z_erro = np.sqrt((oi**2)*(o_err**2) + ((ei**2)*(e_err**2)))
        # Theor_sigma will be the sum of z errors over the sqrt of pos
        th_error = np.sum(z_erro)/np.sqrt(len(positions))

    fitter = LevMarLSQFitter()
    if z_erro is None:
        m_fitted = fitter(model, psi, z)
    else:
        m_fitted = fitter(model, psi, z, weights=1/z_erro)
    info = fitter.fit_info

    result = {}
    # The errors of parameters are assumed to be the sqrt of the diagonal of
    # the covariance matrix
    for i, j, k in zip(m_fitted.param_names, m_fitted.parameters,
                       np.sqrt(np.diag(info['param_cov']))):
        result[i] = {'value': j, 'sigma': k}

    q, u = result['q']['value'], result['u']['value']
    q_err, u_err = result['q']['sigma'], result['u']['sigma']
    p = np.sqrt(q**2 + u**2)
    p_err = np.sqrt(((q/p)**2)*(q_err**2) + ((u/p)**2)*(u_err**2))
    result['p'] = {'value': p, 'sigma': p_err}
    result['sigma_theor'] = th_error
    if z_erro is None:
        result['z'] = {'value': z, 'sigma': np.array([np.nan]*len(z))}
    else:
        result['z'] = {'value': z, 'sigma': z_erro}

    return result
