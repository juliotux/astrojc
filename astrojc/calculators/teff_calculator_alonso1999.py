#Alonso 1999 Teff calculator

import numpy as np

_equations = {1:[0.6388, 0.4065, -0.1117, -2.308e-3, -7.782e-2, -1.2e-2, 0.023],
              2:[0.8323, 9.374e-2, 1.184e-2, 2.351e-2, -0.1392, -1.944e-2, 0.020],
              3:[0.5716, 0.5404, -6.126e-2, -4.862e-2, -0.1777e-2, -7.969e-3, 0.020],
              4:[0.6177, 0.4354, -4.025e-3, 5.204e-2, -0.1127, -1.385e-2, 0.024]}

def calc_teff(equation, color, Feh, color_error=0):
    eq = _equations[equation]
    a = eq[:6]
    st = eq[6]
    X = color

    t=a[0] + a[1]*X + a[2]*(X**2) + a[3]*X*Feh + a[4]*Feh + a[5]*(Feh**2)
    st = st + 3*(color_error/color)
    Teff = 5040/t

    teff_p = 5040/(t+st)
    teff_m = 5040/(t-st)
    sigma = np.mean(np.abs([Teff-teff_p, Teff-teff_m]))
    return Teff, sigma
