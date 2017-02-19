'''
Module to help to fit data without creating a lot of models.
'''

from astropy.modeling import models, fitting
from ..log import log

def _simple_fitting_1d(model, x, y):
    fitter = fitting.LevMarLSQFitter()
    return fitter(model, x, y)

def fit_moffat_1d(x, y, x_0=0, gamma=2.0, alpha=1.5, sky=0):
    model = models.Moffat1D(amplitude=max(y), x_0=x_0, gamma=gamma, alpha=alpha) + \
            models.Const1D(amplitude=sky)
    return _simple_fitting_1d(model, x, y)

def fit_gaussian_1d(x, y, x_0, sigma=1.0, sky=0):
    model = models.Gaussian1D(amplitude=max(y), x_0=x_0, sigma=sigma) + \
            models.Const1D(amplitude=sky)
    return _simple_fitting_1d(model, x, y)
    
