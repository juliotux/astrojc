from os import path, makedirs
import errno

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

