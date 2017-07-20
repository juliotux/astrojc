from collections import OrderedDict
import copy

from .utils import check_class, empty_value


class LinkError(Exception):
    pass


class Output(object):
    '''Class for common output docks.'''
    def __init__(self, name, node, active=True,
                 default_value=empty_value, return_copy=False,
                 link_callback=None, unlink_callback=None):
        if '.' in name or ' ' in name:
            raise ValueError('Output names can contain only numbers and letters,'
                             ' not spaces or dots.')
        self.name = name
        self.node = node
        self.active = active
        self.return_copy = return_copy
        self._default_value = default_value
        self.link_callback = link_callback
        self.unlink_callback = unlink_callback
        self._connectors = []
        self.set_value(self._default_value)

    @property
    def is_linked(self):
        return len(self._connectors) > 0

    def is_linked_to(self, conn):
        return conn in self._connectors

    @property
    def value(self):
        return self.get_value()

    def set_value(self, value):
        '''Sets the value of the output.'''
        self._value = value

    def get_value(self):
        if not self.return_copy:
            return self._value
        else:
            return copy.copy(self._value)

    def reset(self, **kwargs):
        '''Reset the output properties, except by the inputs.'''
        self.name = kwargs.pop('name', self.name)
        self.node = kwargs.pop('node', self.node)
        self.active = kwargs.pop('active', self.active)
        self.return_copy = kwargs.pop('return_copy', self.return_copy)
        self._default_value = kwargs.pop('default_value', self._default_value)

    def link(self, conn):
        '''Link this output to a connector.'''
        if self.is_linked_to(conn):
            return
        if not check_class(conn, 'Connection'):
            raise LinkError('{} is not a valid Connector class.'.format(conn))

        self._connectors.append(conn)
        if not conn.is_linked_to(self):
            conn.set_input(self)

        try:
            self.link_callback()
        except TypeError:
            pass

    def unlink(self, conn):
        '''Unlink this output from a Connector'''
        if not self.is_linked_to(conn):
            return

        self._connectors.remove(conn)
        if conn.is_linked_to(self):
            conn.set_input(None)

        try:
            self.unlink_callback(dock)
        except TypeError:
            pass

    def unlink_all(self):
        '''Unlink all the linked docks.'''
        for i in self._connectors:
            self.unlink(i)
