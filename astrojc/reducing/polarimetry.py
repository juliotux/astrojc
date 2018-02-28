import numpy as np
from scipy.spatial import cKDTree
from astropy.modeling.fitting import LevMarLSQFitter

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


def estimate_normalize(o, e, positions, n_consecutive):
    """Estimate the normalization of a given set of data."""
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


def calculate_polarimetry(o, e, positions, rotation_interval, z=None,
                          retarder='half', o_err=None, e_err=None):
    """Calculate the polarimetry."""
    # if o_err is None or e_err is None:
    #     do_th_error = False
    # else:
    #     # temporary, not compute theoretical error
    #     do_th_error = False
    #
    # if retarder == 'half':
    #     model = HalfWaveModel()
    #     n_consecutive = 4
    # elif retarder == 'quarter':
    #     model = QuarterWaveModel()
    #     n_consecutive = 8
    # else:
    #     raise ValueError('retarder {} not supported.'.format(retarder))
    #
    # if z is None:
    #     k = estimate_normalize(o, e, positions, n_consecutive)
    #     z = (np.array(o)-np.array(e)*k)/(np.array(o)+np.array(e)*k)
    # psi = np.array(positions)*rotation_interval
    #
    # fitter = LevMarLSQFitter()
    # m_fitted = fitter(model, psi, z)
    #
    # return m_fitted

    # result = OrderedDict()
    # errors = OrderedDict()
    # try:
    #     popt, pcov = curve_fit(func, psi, z)
    #     for i in range(len(args)):
    #         result[args[i]] = popt[i]
    #         errors[args[i]] = np.sqrt(np.diag(pcov))[i]
    # except RuntimeError:
    #     for i in args:
    #         result[i] = np.nan
    #         errors[i] = np.nan
    #
    # result['p'] = np.sqrt(result['q']**2 + result['u']**2)
    # errors['p'] = errors['q'] + errors['u']
    #
    # result['theta'] = np.arctan(result['q']/result['u'])
    # errors['theta'] = np.nan  # read references about this calculation
    #
    # if do_th_error:
    #     raise NotImplementedError('Implement the theoretical error of'
    #                               ' polarimetry.')
    # else:
    #     th_error = 0.0

    # TODO: transform errors in sigma like pccdpack
    # return result, errors, th_error
