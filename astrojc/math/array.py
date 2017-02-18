from ._numba_helper import *
import numpy as np

def _xy2r(x, y, data, xc, yc):
    r = np.sqrt((x-xc)**2 + (y-yc)**2)
    return np.ravel(r), np.ravel(data)

def trim_array(data, indices, box_size, position):
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

if use_jit:
    xy2r = autojit(_xy2r)
    #trim_array = autojit(_trim_array)
else:
    xy2r = _xy2r
    #trim_array = _trim_array
