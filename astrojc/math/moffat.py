from ._numba_helper import *

def _moffat_r(r, gamma, alpha, flux, sky):
    return sky + flux*((alpha-1)/(pi*gamma**2))*(1+(r/gamma)**2)**(-alpha)

def _moffat_1d(x, x0, gamma, alpha, flux, sky):
    return sky + flux*((alpha-1)/(pi*gamma**2))*(1 + ((x - x0)**2)/gamma**2)**(-alpha)

def _moffat_2d(x, y, x0, y0, gamma, alpha, flux, sky):
    return sky + flux*((alpha-1)/(pi*gamma**2))*(1 + ((x - x0)**2 + (y - y0)**2)/gamma**2)**(-alpha)

if use_jit:
    moffat_r = vectorize([i+str(tuple([i]*5)).replace('\'','') for i in ('float32', 'float64')],
                         target=numba_target)(_moffat_r)
    moffat_1d = vectorize([i+str(tuple([i]*6)).replace('\'','') for i in ('float32', 'float64')],
                          target=numba_target)(_moffat_1d)
    moffat_2d = vectorize([i+str(tuple([i]*8)).replace('\'','') for i in ('float32', 'float64')],
                          target=numba_target)(_moffat_2d)
else:
    moffat_r = _moffat_r
    moffat_1d = _moffat_1d
    moffat_2d = _moffat_2d

__all__ = ['moffat_r', 'moffat_1d', 'moffat_2d']
