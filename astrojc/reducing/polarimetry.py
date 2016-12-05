from scipy.spatial import cKDTree
from scipy.optimize import curve_fit
import numpy as np

def _estimate_d(x):
    #TODO: Think in use a autocorrelation instead histograms to determine dx
    dxa = []
    for i in range(len(x)):
        for j in range(len(x)):
            if i < j:
                dxa.append(x[i] - x[j])

    dx = 0
    for lim in (100, 30, 5, 3):
        histx = np.histogram(dxa, bins=30, range=[dx-lim, dx+lim])
        mx = np.argmax(histx[0])
        dx = (histx[1][mx]+histx[1][mx+1])/2

    return dx

def estimate_dxdy(x, y):
    return _estimate_d(x), _estimate_d(y)

def match_pairs(positions, dx, dy, tolerance=1.0):
    dt = np.dtype([('o', int), ('e', int)])
    results = np.zeros(len(positions), dtype=dt)
    npairs = 0

    p = np.copy(positions)
    kd = cKDTree(p)

    for i in range(len(p)):
        px = p[i][0]-dx
        py = p[i][1]-dy
        d,j = kd.query((px, py), k=1, eps=tolerance, distance_upper_bound=tolerance, n_jobs=-1)
        if d <= tolerance:
            p[i] = (np.nan, np.nan)
            p[j] = (np.nan, np.nan)
            results[npairs]['o'] = i
            results[npairs]['e'] = j
            npairs = npairs+1
            kd = cKDTree(p)

    return results[:npairs]

def _quarter(psi, q, u, v):
    'Z(I)= Q*cos(2psi(I))**2 + U*sin(2psi(I))*cos(2psi(I)) - V*sin(2psi(I))'
    psi2 = 2*psi
    return q*(np.cos(psi2)**2) + u*np.sin(psi)*np.cos(psi2) - v*np.sin(psi2)

def _half(psi, q, u):
    'Z(I)= Q*cos(4psi(I)) + U*sin(4psi(I))'
    return q*np.cos(4*psi) + u*np.sin(4*psi)

def calculate_polarimetry_quater(fo, fe, psis):
    z = (fe-fo)/(fe+fo)
    return curve_fit(_quarter, psis, z)

def calculate_polarimetry_half(fo, fe, psis):
    k=1
    z = (fo - fe*k)/(fo + fe*k)
    try:
        return curve_fit(_qhalf, psis, z)
    except:
        return (np.nan, np.nan), [[np.nan, np.nan],[np.nan, np.nan]]
