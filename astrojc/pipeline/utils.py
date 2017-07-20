'''General utils for the pipeline.'''

import inspect
from os import path, makedirs
import errno


__all__ = ['check_class', 'EmptyValue', 'empty_value', 'get_all_subclasses',
           'mkdir_p']


def mkdir_p(fname):
    '''
    Function to simulate 'mkdir -p' bash function, with error handling.
    '''
    try:
        makedirs(fname)
    except OSError as exc:
        if exc.errno == errno.EEXIST and path.isdir(fname):
            pass
        else:
            raise exc


class EmptyValue(object):
    def __str__(self):
        return 'Empty Value'

empty_value = EmptyValue()


def check_class(cls, name):
    '''Return true if `name` is a baseclass of the given class.'''
    try:
        mros = inspect.getmro(cls)
    except:
        try:
            mros = inspect.getmro(cls.__class__)
        except:
            return False
    return name in [i.__name__ for i in mros]


def get_all_subclasses(cls):
    all_subclasses = []

    for subclass in cls.__subclasses__():
        all_subclasses.append(subclass)
        all_subclasses.extend(get_all_subclasses(subclass))

    return all_subclasses
