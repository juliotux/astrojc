'''Dereddening of a SED or spectrum, using literature reddeing curves.'''

#TODO: Search and implement from Cardelli 1989, Apj, 345, 245

import numpy as np

# http://adsabs.harvard.edu/abs/2003ApJ...594..279G
ref = {'SMC':{'Wavelength':np.array([0.116,0.119,0.123,0.127,0.131,0.136,0.14,0.145,
                                     0.151,0.157,0.163,0.17,0.178,0.186,0.195,0.205,
                                     0.216,0.229,0.242,0.258,0.276,0.296,0.37,0.44,
                                     0.55,0.65,0.81,1.25,1.65,2.198])*10000,
              'Ax/Av':np.array([6.992,6.436,6.297,6.074,5.795,5.575,5.272,5,4.776,
                                4.472,4.243,4.013,3.866,3.637,3.489,3.293,3.161,2.947,
                                2.661,2.428,2.22,2,1.672,1.374,1,0.801,0.567,0.131,
                                0.169,0.016])},
       'LMC':{'Wavelength':np.array([0.119,0.123,0.127,0.131,0.136,0.14,0.145,0.151,
                                     0.157,0.163,0.17,0.178,0.186,0.195,0.205,0.216,
                                     0.229,0.242,0.258,0.276,0.296,0.37,0.44,0.55,
                                     1.25,1.65,2.198])*10000,
              'Ax/Av':np.array([3.467,3.366,3.374,3.231,3.118,2.983,2.874,2.787,
                                2.668,2.607,2.598,2.566,2.565,2.646,2.846,2.967,
                                2.771,2.391,2.149,1.969,1.786,1.518,1.293,1,0.257,
                                0.186,0.03])}}

def dered_flux(wavelength, flux, location):
    corr_index = np.interp(wavelength, ref[location]['Wavelength'], ref[location]['Ax/Av'])
    return flux*(10**(corr_index/2.5))

def dered_mag(wavelength, mag, location):
    corr_index = np.interp(wavelength, ref[location]['Wavelength'], ref[location]['Ax/Av'])
    return mag - corr_index
