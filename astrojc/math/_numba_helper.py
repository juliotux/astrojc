try:
    from numba import vectorize, autojit
    from math import sin, cos, exp, pi, sqrt
    use_jit = True
    numba_target = 'parallel'
except:
    from numpy import sin, cos, exp, pi, sqrt
    use_jit = False
    numba_target = None
