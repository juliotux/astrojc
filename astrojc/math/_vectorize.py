"""Create convenient functions to vectorize functions, even without numba"""
import numpy as np
from ..logging import log as logger

try:
    from numba import vectorize
    vectorize_target = 'parallel'
    logger.debug('Numba found, using parallel Moffat kernel.')
except ModuleNotFoundError:
    # Handle numba not installed
    import inspect

    def vectorize(arg, target):
        def decorator(func):
            sign = inspect.signature(func)
            ninp = len(sign.parameters)
            return np.frompyfunc(func, ninp, 1)
        return decorator

    logger.info('Numba not found, using a serialized version of Moffat func.')
    vectorize_target = None
