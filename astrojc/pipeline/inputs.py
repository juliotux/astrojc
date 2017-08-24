from collections import OrderedDict
import copy

from .utils import check_class, empty_value


class LinkError(Exception):
    pass


class Input(object):
    '''Class for common input docks.'''
    def __init__(self, name, node, active=True,
                 default_value=empty_value, optional=False,
                 link_callback=None, unlink_callback=None):
        if '.' in name or ' ' in name:
            raise ValueError('Input names can contain only numbers and letters,'
                             ' not spaces or dots.')
        self.name = name
        self.node = node
        self.active = active
        self.optional = optional
        self._default_value = default_value
        self.link_callback = link_callback
        self.unlink_callback = unlink_callback
        self._value = self._default_value

        self._connector = None

    @property
    def is_linked(self):
        return not self._connector is None

    def is_linked_to(self, conn):
        return conn == self._connector

    @property
    def value(self):
        return self.get_value()

    @property
    def value_ready(self):
        if self.optional:
            return True
        return self.value != empty_value

    @property
    def have_value(self):
        return self.get_value() not in (empty_value, None)

    def set_value(self, value):
        print('Setting value {} for this node.'.format(value))
        self._value = value

    def get_value(self):
        try:
            return self._connector.value
        except:
            try:
                return self._value
            except:
                return self._default_value

    def reset(self, **kwargs):
        '''Reset the input properties, except by the output.'''
        self.name = kwargs.pop('name', self.name)
        self.node = kwargs.pop('node', self.node)
        self.active = kwargs.pop('active', self.active)
        self.return_copy = kwargs.pop('return_copy', self.return_copy)
        self._default_value = kwargs.pop('default_value', self._default_value)
        self._value = kwargs.pop('default_value', self._default_value)

    def link(self, conn):
        '''Link this input to an Connector.'''
        if not check_class(conn, 'Connection') or self._connector is not None:
            raise LinkError('You can connect just one Connector to this'
                            ' node input.')

        self._connector = conn
        if not conn.is_linked_to(self):
            conn.set_input(self)

        try:
            self.link_callback()
        except TypeError:
            pass

    def unlink(self, conn):
        '''Unlink this input from an Connector'''
        if not self.is_linked_to(conn):
            return

        self._connector = None
        if conn.is_linked_to(self):
            conn.set_input(None)

        try:
            self.unlink_callback(dock)
        except TypeError:
            pass

    def unlink_all(self):
        self.unlink(self._connector)


###############################################################################
# Implementation of input data types
###############################################################################

class FloatInput(Input):
    '''Input for single float values.'''
    def __init__(self, name, node, *args, **kwargs):
        Input.__init__(self, name, node, *args, **kwargs)
        self.dtype = 'float'


class IntInput(Input):
    '''Input for single integer values.'''
    def __init__(self, name, node, *args, **kwargs):
        Input.__init__(self, name, node, *args, **kwargs)
        self.dtype = 'int'


class StringInput(Input):
    '''Input for single string values.'''
    def __init__(self, name, node, *args, **kwargs):
        Input.__init__(self, name, node, *args, **kwargs)
        self.dtype = 'str'


class TextInput(Input):
    '''Input for multiline text values.'''
    def __init__(self, name, node, *args, **kwargs):
        Input.__init__(self, name, node, *args, **kwargs)
        self.dtype = 'str'


class BoolInput(Input):
    '''Input for single boolean values.'''
    def __init__(self, name, node, *args, **kwargs):
        Input.__init__(self, name, node, *args, **kwargs)
        self.dtype = 'bool'


class FileInput(Input):
    '''Input for single file name or file buffer.'''
    def __init__(self, name, node, *args, **kwargs):
        Input.__init__(self, name, node, *args, **kwargs)
        self.dtype = 'file'


class MultFileInput(Input):
    '''Input for multiple (or single) file names or file buffers.'''
    def __init__(self, name, node, *args, **kwargs):
        Input.__init__(self, name, node, *args, **kwargs)
        self.dtype = 'list(file)'
