import numpy as np

__all__ = ['uvby2UBV']

def uvby2UBV(u, v, b, y):
    '''
    Convert Stromgren uvby to Johnson UBV magnitudes, using the polinomial
    relations from Harmanec & Bozic, 2001, A&A, 369, 1140.
    '''
    V=y
    B=V + 1.41694*(b-v) + 0.07010*(u-b) + 0.57145*(b-v)**2 - 0.60399*(b-v)**3 - 0.10118
    U=B + 0.66567*(u-b) - 0.09718*(b-v) + 0.24407*(b-v)**2 + 0.29340*(b-v)**3 - 0.91958
    return (U,B,V)
