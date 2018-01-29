import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import curve_fit
from collections import OrderedDict


def estimate_dxdy(x, y, steps=[100, 30, 5, 3], bins=30):
    """Estimate the x and y distance between the pairs of points"""
    dxa = []
    dya = []

    def _find_max(d, steps=[100, 30, 5, 3], bins=30):
        dx = 0
        for lim in (np.max(d), *steps):
            histx = np.histogram(d, bins=bins, range=[dx-lim, dx+lim])
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


def _quarter(psi, q, u, v, k):
    '''Z= Q*cos(2psi)**2 + U*sin(2psi)*cos(2psi) - V*sin(2psi)'''
    # I introduced a k constant for handle normalization in the fitting
    psi2 = 2*psi
    z = q*(np.cos(psi2)**2) + u*np.sin(psi)*np.cos(psi2) - v*np.sin(psi2) + k
    return z


def _half(psi, q, u, k):
    '''Z(I)= Q*cos(4psi(I)) + U*sin(4psi(I))'''
    # I introduced a k constant for handle normalization in the fitting
    return q*np.cos(4*psi) + u*np.sin(4*psi) + k


def calculate_polarimetry(o, e, psi, retarder='half', o_err=None, e_err=None):
    """Calculate the polarimetry."""
    if o_err is None or e_err is None:
        do_th_error = False
    else:
        do_th_error = True

    if retarder == 'half':
        func = _half
        args = ['q', 'u', 'k']
    elif retarder == 'quarter':
        func = _quarter
        args = ['q', 'u', 'v', 'k']
    else:
        raise ValueError('retarder {} not supported.'.format(retarder))

    z = np.subtract(o, e)/np.sum(o, e)
    result = OrderedDict()
    errors = OrderedDict()
    try:
        popt, pcov = curve_fit(func, psi, z)
        for i in range(len(args)):
            result[args[i]] = popt[i]
            errors[args[i]] = np.sqrt(np.diag(pcov))[i]
    except RuntimeError:
        for i in args:
            result[i] = np.nan
            errors[i] = np.nan

    if do_th_error:
        raise NotImplementedError('Implement the theoretical error of'
                                  ' polarimetry.')
    else:
        th_error = 0.0

    # TODO: transform errors in sigma like pccdpack
    return result, errors, th_error
