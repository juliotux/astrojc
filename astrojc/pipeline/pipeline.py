from os import path
import os
import glob
import imp

import json
from collections import OrderedDict
import copy
import time
import inspect
from six import string_types

from .context import Context
from .nodes import Node
from .utils import check_class, get_all_subclasses, mkdir_p


nodes_dirs = [path.expanduser('~/.astrojc/nodes'),
              path.expanduser('~/.local/share/astrojc/nodes'),
              '/usr/local/share/astrojc/nodes',
              '/usr/share/astrojc/nodes']

def read_json(fname):
    f = open(fname, 'r').read()
    j = json.JSONDecoder(object_pairs_hook=OrderedDict).decode(f)
    return j

class Pipeline(object):
    '''General class to manipulate pipelines.'''
    def __init__(self, filename=None, nodes_dir=None):
        self._name = 'pipeline_context'
        self._version = None
        self._ctx = Context('pipeline_context')
        self._nodes_dir = nodes_dirs
        self._products = OrderedDict()
        self._connections = []
        self._redirects = OrderedDict()
        self._valid_nodes = OrderedDict()

        if nodes_dir is not None:
            d = path.realpath(path.join(path.dirname(filename), nodes_dir))
            self._nodes_dir.append(d)


        self._path = None
        self._prod_path = None

        if filename is not None:
            #self._nodes_dir.append(self._path)
            #self._nodes_dir.append(path.join(self._path, 'nodes'))
            self.gen_from_file(filename)

        self.update_node_list()

        self._variables = {'%PRODUCT_FILE%' : '',
                           '%PRODUCT_NAME%' : '',
                           '%RAW_DIR%' : '',
                           '%PROD_DIR%' : '',
                           '%CALIB_DIR%' : '',
                           '%PROD_TYPE%' : '',
                           '%TIME%' : time.asctime().replace(' ', '_'),
                           '%PIPELINE%' : self._name + '_' + str(self._version)}
        self._save_variables = copy.copy(self._variables)

        self._root_state = self._ctx.vars()
        self._prod_state = self._ctx.vars()
        self._local_state = self._ctx.vars()

    def _sanitize_values(self, values):
        '''Sanitize the variables names in values.'''
        if isinstance(values, (list, tuple)):
            for i in range(len(values)):
                values[i] = self._sanitize_values(values[i])
        elif isinstance(values, (dict, OrderedDict)):
            for i,v in values.items():
                values[i] = self._sanitize_values(v)
        else:
            try:
                for k in self._variables.keys():
                    values = values.replace(k, str(self._variables[k]))
            except:
                return values
        return values

    def gen_from_file(self, fname):
        '''Generates the pipeline from a config file.'''
        j = read_json(fname)
        self._name = j['name']
        self._version = j['version']
        self._ctx = Context(self._name)

        self.update_node_list()

        for i,v in j['nodes'].items():
            if v not in self._valid_nodes.keys():
                raise ValueError('{} is not a valid node.'.format(v))
            try:
                self._ctx.create_node(i, self._valid_nodes[v])
            except:
                raise

        for i,v in j['connections'].items():
            self._connections.append((i, v))
            try:
                self._ctx.link(i, v)
            except:
                try:
                    self._ctx.link(v, i)
                except:
                    raise ValueError('Could not connect {} to {}'.format(i, v))

        for i,v in j['redirects'].items():
            self._redirects[i] = v

        self._path = path.dirname(path.realpath(fname))

        self._variables = {'%PRODUCT_FILE%' : '',
                           '%PRODUCT_NAME%' : '',
                           '%RAW_DIR%' : '',
                           '%PROD_DIR%' : '',
                           '%CALIB_DIR%' : '',
                           '%PROD_TYPE%' : '',
                           '%TIME%' : time.asctime().replace(' ', '_'),
                           '%PIPELINE%' : self._name + '_' + str(self._version)}
        self._save_variables = copy.copy(self._variables)

        self._root_state = self._ctx.vars()

    def load_defaults_file(self, fname):
        '''Loads the default values from a config file.'''
        try:
            j = read_json(fname)
        except:
            raise ValueError('Could not load json file.')
        if j['pipeline_name'] != self._name:
            raise ValueError('This product is not for this pipeline.')

        for i,v in j['values'].items():
            self.update_parameter(i, v)
        self._root_state = self._ctx.vars()

    def load_products_file(self, fname):
        '''Generates the product list from a config file.'''
        try:
            j = read_json(fname)
        except:
            raise ValueError('Could not load json file.')
        if j['pipeline_name'] != self._name:
            raise ValueError('This product is not for this pipeline.')

        for i in j['products'].keys():
            a = copy.copy(j['values'])
            a.update(j['products'][i]['values'])
            self._products[i] = a
            self._products[i].status = 'ready'

        self._prod_path = path.dirname(path.realpath(fname))
        self._prod_state = self._ctx.vars()

    def update_parameter(self, key, value):
        '''Updates a parameter in the pipeline.'''
        self._ctx[key] = value
        self._local_state = self._ctx.vars()

    def update_node_list(self):
        '''Searches in all paths the valid nodes and put them in the dict'''
        self._valid_nodes = OrderedDict()
        for i in self._nodes_dir:
            self._ctx.log.debug('Scanning {} to search nodes.'.format(i))
            for f in glob.glob(path.join(i, '*.py')):
                if path.basename(f)[:2] == '__':
                    continue
                self._ctx.log.debug('Finding nodes in {}'.format(f))
                _tmp_mod = imp.load_source('nodes', f)
                for name, obj in inspect.getmembers(_tmp_mod):
                    if inspect.isclass(obj) and obj in get_all_subclasses(Node):
                        self._ctx.log.info('Adding {} to valid node list.'.format(name))
                        self._valid_nodes[name] = obj

    def run_product(self, product):
        '''Run a single product.'''
        self._ctx.log.info('Reseting the pipeline context.')
        # Reset the context
        self._variables = copy.copy(self._save_variables)

        for i,v in self._root_state.items():
            self._ctx[i] = v

        # Reset the product
        for i,v in self._prod_state.items():
            self._ctx[i] = v

        # Reset the the custom local state
        for i,v in self._local_state.items():
            self._ctx[i] = v

        self._ctx.log.info('Setting {} product variables.'.format(product))

        _tmp_prod = copy.copy(self._products[product])
        _tmp_prod.update(self._variables)
        _keys = list(_tmp_prod.keys())

        _tmp_prod['%PRODUCT_NAME%'] = product

        for i in _keys:
            tmp = i
            if tmp in self._redirects.keys():
                while tmp in self._redirects.keys():
                    self._ctx.log.debug('Redirecting {} to {}.'.format(tmp, self._redirects[tmp]))
                    _tmp_prod[self._redirects[tmp]] = _tmp_prod[tmp]
                    tmp = self._redirects[tmp]

        for i,v in _tmp_prod.items():
            if i in self._variables.keys():
                self._variables[i] = v

        self._variables['%RAW_DIR%'] = path.realpath(path.join(self._prod_path,
                                                     self._variables['%RAW_DIR%']))
        self._variables['%CALIB_DIR%'] = path.realpath(path.join(self._prod_path,
                                                       self._variables['%CALIB_DIR%']))
        self._variables['%PROD_DIR%'] = path.realpath(path.join(self._prod_path,
                                                      self._variables['%PROD_DIR%']))
        if _tmp_prod['type'] == 'calib':
            self._variables['%PRODUCT_FILE%'] = path.join(self._variables['%CALIB_DIR%'],
                                                          self._variables['%PRODUCT_NAME%'])
        else:
            self._variables['%PRODUCT_FILE%'] = path.join(self._variables['%PROD_DIR%'],
                                                          self._variables['%PRODUCT_NAME%'])
        self._variables['%PROD_TYPE%'] = self._products[product]['type']

        _tmp_prod.update(self._variables)

        _keys = list(_tmp_prod.keys())
        for i in _keys:
            tmp = i
            if tmp in self._redirects.keys():
                while tmp in self._redirects.keys():
                    self._ctx.log.debug('Redirecting {} to {}.'.format(tmp, self._redirects[tmp]))
                    _tmp_prod[self._redirects[tmp]] = _tmp_prod[tmp]
                    tmp = self._redirects[tmp]

        for i,v in _tmp_prod.items():
            if i in self._ctx.vars().keys():
                self._ctx[i] = v

        #sanitize values
        for i,v in self._ctx.vars().items():
            self._ctx[i] = self._sanitize_values(v)

        for i in ('%PROD_DIR%', '%CALIB_DIR%'):
            mkdir_p(self._variables[i])

        self._ctx.run()

        self._products.pop(product)


    def run(self, products='all'):
        '''Run the pipeline for a list of product names or 'all' prod. loaded'''

        if len(products) == 0:
            self._ctx.error('No products found to run!')
            return

        if products == 'all':
            return self.run(products=[i for i in self._products.keys()])

        for i in products:
            self.run_product(i)
