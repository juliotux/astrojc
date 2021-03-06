'''
Simple class to emulate signals in python.
'''

class NotConnectedException(Exception):
    pass

class MySignal(object):
    '''
    This class implements a simple way to handle signals in Python, without use
    qt, but with similar (simpler) behavior. No slot is supported.
    '''
    def __init__(self, raise_error=True):
        self._emit_function = None
        self.raise_error = raise_error

    @property
    def is_connected(self):
        return self._emit_function != None

    def connect(self, function):
        '''
        Connect the function that will be executed when the signal is emitted.
        '''
        self._emit_function = function

    def disconnect(self):
        '''
        Disconnect the function exectution when the signal is emitted.
        '''
        self._emit_function = None

    def emit(self, *args, **kwargs):
        '''
        Execute the connected functions, passing to them the args and kwargs.
        '''
        if self._emit_function is not None:
            return self._emit_function(*args, **kwargs)
        elif self.raise_error:
            raise NotConnectedException('This signal is not connected with any '
                                        'emit function.')
        else:
            return

    def __call__(self, *args, **kwargs):
        return self.emit(*args, **kwargs)
