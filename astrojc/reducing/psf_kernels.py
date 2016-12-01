'''
Kernels for psf fitting in the psf_fitting module. May be optimized with numba.
'''

try:
    from numba import vectorize
    from math import sin, cos, exp
    _numbajit = True
    _target = 'parallel'
except:
    from numpy import sin, cos, exp
    _numbajit = False
    _target = None

def _moffat_radial(r, gamma, alpha, amplitude, sky):
    return sky + amplitude*(1+(r/gamma)**2)**(-alpha)

def _moffat_integrate(r, theta, amplitude, gamma, alpha, sky):
    return r*amplitude*(1+(r/gamma)**2)**(-alpha)

def _moffat_spatial(x, y, x0, y0, gamma, alpha, amplitude):
    return amplitude*(1 + ((x - x0)**2 + (y - y0)**2)/gamma**2)**(-alpha)

def _gaussian_radial(r, sigma, amplitude, sky):
    return sky + amplitude*exp(-0.5*(r/sigma)**2)

def _gaussian_spatial(x, y, x0, y0, sigma_x, sigma_y, theta, amplitude):
    cost2 = cos(theta)**2
    sint2 = sin(theta)**2
    sin2t = sin(2*theta)
    sigx2 = 2*sigma_x**2
    sigy2 = 2*sigma_y**2
    a = (cost2/sigx2) + (sint2/sigy2)
    b = -(sin2t/(2*sigx2)) + (sin2t/(2*sigy2))
    c = (sint2/sigx2) + (cost2/sigy2)
    xi = x - x0
    yi = y - y0
    return amplitude * exp(-(a*xi**2 + 2*b*xi*yi + c*yi**2))

if _numbajit:
    moffat_radial = vectorize('float64(float64,float64,float64,float64,float64)',
                              target=_target)(_moffat_radial)
    moffat_spatial = vectorize('float64(float64,float64,float64,float64,float64,float64,float64)',
                               target=_target)(_moffat_spatial)
    moffat_integrate = vectorize('float64(float64,float64,float64,float64,float64,float64)',
                                 target=_target)(_moffat_integrate)
    gaussian_radial = vectorize('float64(float64,float64,float64,float64)',
                                 target=_target)(_gaussian_radial)
    gaussian_spatial = vectorize('float64(float64,float64,float64,float64,float64,float64,float64,float64)',
                                 target=_target)(_gaussian_spatial)
else:
    moffat_radial = _moffat_radial
    moffat_spatial = _moffat_spatial
    moffat_integrate = _moffat_integrate
    gaussian_radial = _gaussian_ragial
    gaussian_spatial = _gaussian_spatial
                                     
