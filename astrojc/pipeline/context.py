import os
import sys

from collections import OrderedDict
import logging

from astropy.logger import log

from .nodes import Node
from .connections import Connection
from .utils import check_class


__all__ = ['Context', 'Processor']


class Processor(object):
    '''Processes a list of nodes.'''
    def __init__(self, nodes):
        self._nodes = nodes
        self._done_nodes = []
        self._nodes_to_run = []

    def run(self):
        self._nodes_to_run = list(self._nodes.values())

        while not len(self._nodes_to_run) == 0:
            previous_done = len(self._done_nodes)
            for i in self._nodes_to_run:
                if i.is_ready:
                    log.info('Running node {}'.format(i.name))
                    i.__run_node__()
                    self._done_nodes.append(i)
                    self._nodes_to_run.remove(i)
            if previous_done == len(self._done_nodes):
                log.error('Infinite loop with not ready nodes. Breaking.')
                break

class Context(object):
    '''Handles a full context fo nodes, node groups, variables and runners.'''
    def __init__(self, name, nodes=[], runner=None,
                 properties = OrderedDict()):
        if '.' in name or ' ' in name:
            raise ValueError('Context names can contain only numbers and letters,'
                             ' not spaces or dots.')
        self._name = name
        self._nodes = OrderedDict()
        self._runner = runner
        self._properties = properties
        self._connectors = []

        for n in nodes:
            self.add_node(n)

        self.log = log

    @property
    def name(self):
        return self._name

    def get_property(self, key, value=None):
        '''Get a property from this context or a node inside this context.'''
        #TODO: Include get a property of nodes and spliting properties in groups
        #      with points.
        return self._properties.get(key, value)

    def set_property(self, key, value):
        '''Set the value of a property in this context.'''
        self._properties[key] = value

    def add_node(self, node):
        '''Adds an existent node instance to the context.'''
        if not check_class(node, 'Node'):
            self.log.debug('Node class: {}'.format(node.__class__.__class__))
            self.log.error('You tried to add a non Node instance as a node.')
            return

        self.log.info('Adding the Node {} to {} context.'.format(node.name, self._name))
        self._nodes[node.name] = node
        self._nodes[node.name]._ctx = self

    def create_node(self, name, node_class, **node_kwargs):
        '''Creates a new node in this context.'''
        if '.' in name:
            n_name = name.split('.')[0]
            if n_name not in self._nodes.keys():
                raise ValueError('Supernode name {} not valid.'.format(n_name))
            return self._nodes[n_name].create_node('.'.join(name.split('.')[1:]),
                                                   node_class)

        try:
            node = node_class(name, context=self, **node_kwargs)
        except Exception as e:
            self.log.error('We could not create the node with this class')
            node = None
            raise e
        if not check_class(node, 'Node'):
            self.log.error('The node you tryied to create has not a valid Node'
                           ' class.')
            return

        try:
            self.add_node(node)
        except:
            self.log.error('Could not add the Node {} to {} context.'.format(name, self._name))
            raise

        return node

    def del_node(self, node=None, name=None):
        if name is not None:
            if name in self._nodes.keys():
                self._nodes[name].close()
                del self._nodes[name]
                self.log.info('Removing Node {} from this context.'.format(name))
                return
            else:
                self.log.info('Node {} not in this context.'
                              'Ignoring it.'.format(name))
                return
        if node is not None:
            if not check_class(node_class, 'Node'):
                self.log.error('The node you are trying to remove is not a '
                               'Node instance.')
                return
            if node not in self._nodes.values():
                self.info('Node not in this context. Ignoring it.')
                return
            node.close()
            for i,v in self._nodes.items():
                if node == v:
                    del self._nodes[i]
                    self.log.info('Removing Node {} from this context.'.format(i))
                    return

    def run(self):
        self.log.info('Running context {}'.format(self.name))
        p = Processor(self._nodes)
        p.run()

    def get_dock(self, key):
        '''Return the dock object from a given key.'''
        keys = key.split('.')
        name = keys[0]
        if len(keys) > 1:
            arg = keys[1:]
        else:
            arg = []
        if name in self._nodes.keys():
            try:
                return self._nodes[name].get_dock(".".join(arg))
            except:
                raise ValueError('{} key doesn\'t correspond to any dock '
                                 'from this context.'.format(key))

    def link(self, input_key, output_key):
        '''Creates a link between one input and one output.'''
        i_node = self[input_key.split('.')[0]]
        i_dock = i_node.get_dock('.'.join(input_key.split('.')[1:]))
        o_node = self[output_key.split('.')[0]]
        o_dock = o_node.get_dock('.'.join(output_key.split('.')[1:]))

        self._connectors.append(Connection(i_dock, o_dock))
        self.log.info('Linked {} to {}'.format(input_key, output_key))

    def unlink(self, input_key, output_key):
        '''Removes the link between one input and one output.'''
        i_node = self[input_key.split('.')[0]]
        i_dock = i_node.get_dock('.'.join(input_key.split('.')[1:]))
        o_node = self[output_key.split('.')[0]]
        o_dock = o_node.get_dock('.'.join(output_key.split('.')[1:]))

        for i in [c for c in self._connectors if c._input == i_dock and
                  c._output == o_dock]:
            i.unlink()
            del i

        self.log.info('Unlinked {} from {}'.format(input_key, output_key))

    def vars(self):
        '''Return a dict with the keys of inputs, outputs and props'''
        items = OrderedDict()

        # Iterate over nodes
        for i,n in self._nodes.items():
            for j, v in n.vars().items():
                items['.'.join([str(i), str(j)])] = v

        for i,p in self._properties.items():
            items['.'.join(['property', str(i)])] = p

        return items


    def __getitem__(self, key):
        '''Get a node or a property from this context.'''
        if key[0:9] == 'property.':
            return self._properties[key[9:]]
        else:
            keys = key.split('.')
            name = keys[0]
            if len(keys) > 1:
                arg = keys[1:]
            else:
                arg = []
            if name in self._nodes.keys():
                return self._nodes[name][".".join(arg)]
            else:
                raise ValueError('{} key doesn\'t correspond to any property '
                                 'or node from this context.'.format(key))

    def __setitem__(self, key, value):
        '''Used for set a property, and node items (property, input or output)'''
        if key[0:9] == 'property.':
            self._properties[key[9:]] = value
        else:
            keys = key.split('.')
            name = keys[0]
            if len(keys) > 1:
                arg = keys[1:]
            else:
                arg = []
            if name in self._nodes.keys():
                self._nodes[name][".".join(arg)] = value
            else:
                raise ValueError('{} key doesn\'t correspond to any property '
                                 'or node from this context.'.format(key))
