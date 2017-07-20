'''Class for handle nodes of the pipeline. This is not intended to be changed
during the execution and not intended to be executed in realtime, but only when
you press a 'Play' button or excetue a play command.'''

from collections import OrderedDict
import copy
import logging

from .utils import check_class, empty_value
from .inputs import *
from .outputs import *

__all__ = ['Node', 'NodeNotReadyError']


class NodeNotReadyError(Exception):
    pass


class Node(object):
    '''Class for common pipeline nodes.'''
    def __init__(self, name, x=None, y=None, parent=None,
                 resizable=True, deletable=True, context=None, enable=True):
        if '.' in name or ' ' in name:
            raise ValueError('Node names can contain only numbers and letters,'
                             ' not spaces or dots.')
        self._name = name
        self._inputs = OrderedDict()
        self._outputs = OrderedDict()
        self._parent = parent
        self._resizable = resizable
        self._deletable = deletable
        self._ctx = context
        self._enable = enable

        self._properties = OrderedDict()

        try:
            self.log = self._ctx.log
        except:
            self.log = logging.getLogger(self.name)

        self._state = 'idle' # Node running state: 'idle', 'waiting', 'running',
                             #                     'stoped', 'done', 'error'

        self.setup()

    @property
    def name(self):
        if self._name is None:
            return self.__class__.__name__
        return self._name

    @property
    def parent(self):
        return self._parent

    def set_name(self, name):
        self._name = name

    def _check_ready(self):
        #TODO: redo this
        #if not self._enable:
        #    return False
        for i,v in self._inputs.items():
            if not v.value_ready and not v.optional:
                self.log.debug('{} node not ready due {} input.'
                               ''.format(self.name, i))
                return False
        return True

    @property
    def is_ready(self):
        '''Returns True if the node is ready to be executed.'''
        return self._check_ready()

    def setup(self):
        '''Function called for setup a node, to include inputs, outputs and properties.'''

    def run(self):
        '''Function to run this node when ready.'''

    def __run_node__(self):
        '''Protected function to run a node securely.'''
        if not self.is_ready:
            raise NodeNotReadyError()

        try:
            self.run()
        except:
            self._state = 'error'
            raise

        self._state = 'done'

    def get_state(self):
        return self._state

    def set_state(self, state):
        if state in ['idle', 'waiting', 'running', 'stoped', 'done', 'error']:
            self._state = state
        else:
            raise ValueError('State not recognized.')

    def add_input(self, name, input_class=Input, default_value=empty_value,
                  active=True, optional=False):
        '''Adds an input to this node.'''
        if not check_class(input_class, 'Input'):
            raise ValueError('This input is not a valid input class.')
        try:
            inp = input_class(name=name, node=self, default_value=default_value,
                              active=active, optional=optional)
        except Exception as e:
            raise ValueError('Could not create this input due: {}'.format(e))

        self._inputs[name] = inp
        return inp

    def remove_input(self, name):
        '''Remove an input from this node.'''
        self._inputs[name].unlink_all()
        del self._inputs[name]

    def get_inputs(self):
        '''Return a list of the inputs.'''
        return list(self._inputs)

    def add_output(self, name, output_class=Output, default_value=empty_value,
                   active=True, return_copy=False):
        '''Adds an output to this node.'''
        if not check_class(output_class, 'Output'):
            raise ValueError('This input is not a valid input class.')
        try:
            out = output_class(name=name, node=self, default_value=default_value,
                               active=active, return_copy=return_copy)
        except Exception as e:
            raise ValueError('Could not create this output due: {}'.format(e))

        self._outputs[name] = out
        return out

    def remove_output(self, name):
        '''Remove an output from this node.'''
        self._outputs[name].unlink_all()
        del self._outputs[name]

    def get_outputs(self):
        '''Return a list of the outputs.'''
        return list(self._outputs)

    def unlink_all(self):
        '''Unlink all inputs and outputs.'''
        for i,v in self._inputs.items():
            v.unlink_all()
        for o,v in self._outputs.items():
            v.unlink_all()

    def close(self):
        '''Closes the node.'''
        self.unlink_all()

    def get_dock(self, key):
        '''Gets the object of a dock with a given key.'''
        keys = key.split('.')
        keyname = keys[0]
        name = ".".join(keys[1:])
        if keyname == 'input':
            if name in self._inputs.keys():
                return self._inputs[name]
            else:
                raise ValueError('Input name {} not identified.'.format(name))
        elif keyname == 'output':
            if name in self._outputs.keys():
                return self._outputs[name]
            else:
                raise ValueError('Output name {} not identified.'.format(name))
        else:
            raise ValueError('Only inputs and outputs can be obtained with '
                             'this method.')

    def vars(self):
        '''Return a dict with the keys of inputs, outputs, props and nodes.'''
        items = OrderedDict()

        for i,p in self._inputs.items():
            items['.'.join(['input', str(i)])] = p.value

        for i,p in self._outputs.items():
            items['.'.join(['output', str(i)])] = p.value

        for i,p in self._properties.items():
            items['.'.join(['property', str(i)])] = p

        return items

    def __getitem__(self, key):
        '''Get the value of input, output or property using keys.'''
        if key == '':
            return self
        keys = key.split('.')
        keyname = keys[0]
        if len(keys) > 1:
            arg = keys[1:]
        else:
            arg = []
        name = ".".join(arg)
        if keyname in ('input', 'output'):
            return self.get_dock(key).value
        elif keyname == 'property':
            if name in self._properties.keys():
                return self._properties[name]
            else:
                raise ValueError('Property name {} not identified.'.format(name))
        else:
            raise ValueError("{} key is not valid, please use keys in format "
                             "'input.name', 'output.name' or 'property.name'."
                             "".format(key))

    def __setitem__(self, key, value):
        '''Used for set property or output. Input not supported.'''
        keys = key.split('.')
        keyname = keys[0]
        name = ".".join(keys[1:])
        if keyname == 'input':
            if name in self._inputs.keys():
                self.get_dock(key).set_value(value)
            else:
                raise ValueError('Input name {} not identified.'.format(name))
        elif keyname == 'output':
            if name in self._outputs.keys():
                self.get_dock(key).set_value(value)
            else:
                raise ValueError('Output name {} not identified.'.format(name))
        elif keyname == 'property':
            if name in self._properties.keys():
                self._properties[name] = value
            else:
                raise ValueError('Property name {} not identified.'.format(name))
        else:
            raise ValueError("{} key is not valid, please use keys in format "
                             "'input.name', 'output.name' or 'property.name'."
                             "".format(key))
