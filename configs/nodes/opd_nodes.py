import os
from os import path
import sys
sys.path.append('../')

import numpy as np
from copy import copy
import json

from astropy.io import fits, ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astroscrappy import detect_cosmics
from ccdproc import CCDData, block_reduce, trim_image, combine, ccd_process

from tempfile import NamedTemporaryFile
from shutil import rmtree
import functools
import glob
from six import string_types

from collections import OrderedDict

from astrojc.pipeline.core.inputs import *
from astrojc.pipeline.core.outputs import *
from astrojc.pipeline.core.nodes import Node, empty_value
from astrojc.pipeline.core.context import Processor, Context


__all__ = ['ForSuperNode', 'ReadCCDData', 'WriteCCDData', 'CombineImages',
           'BinAndTrim', 'LACosmic', 'CCDProcess']


################################################################################


class ForSuperNode(Node):
    '''This is a temporary node that will act as node group with for function.'''
    '''This time, the nodes have to be created outside this.'''
    def __init__(self, name, subcontext = None, **kwargs):
        Node.__init__(self, name, **kwargs)
        if isinstance(subcontext, Context):
            self._subctx = subcontext
        elif subcontext == None:
            self._subctx = Context(name)
        else:
            raise ValueError('Invalid subcontext.')

        self._setup_subctx()

    def _setup_subctx(self):
        self._subctx.log = self._ctx.log
        self.add_node = self._subctx.add_node
        self.del_node = self._subctx.del_node
        self.create_node = self._subctx.create_node

    def setup(self):
        self.add_input('[i]')
        self.add_output('i')
        self.add_input('[j]', optional = True)
        self.add_output('j')
        self.add_input('[k]', optional = True)
        self.add_output('k')

        self.add_input('result', optional = True)
        self.add_output('[results]')

    def run(self):
        iter_size = len(self['input.[i]'])

        have_j = False
        have_k = False

        if self.get_dock('input.[j]').have_value:
            if len(self['input.[j]']) != iter_size:
                raise ValueError('[i] and [j] list do not match in size.')
            else:
                have_j = True

        if self.get_dock('input.[k]').have_value:
            if len(self['input.[k]']) != iter_size:
                raise ValueError('[i] and [k] list do not match in size.')
            else:
                have_k = True

        results = [None] * iter_size

        for i in range(iter_size):
            self.log.debug('Running {} {}/{}'.format(self.name, i+1, iter_size))
            self['output.i'] = self['input.[i]'][i]
            if have_j:
                self['output.j'] = self['input.[j]'][i]
            if have_k:
                self['output.k'] = self['input.[k]'][i]

            self._subctx.run()

            results[i] = self['input.result']
            self['output.i'] = empty_value
            self['output.j'] = empty_value
            self['output.k'] = empty_value
            for j in self._subctx._nodes.values():
                j.set_state('idle')

        self['output.[results]'] = results
        for i in self._subctx._nodes.values():
            i.set_state('done')

    def _check_ready(self):
        self['output.i'] = self['input.[i]']
        self['output.j'] = self['input.[j]']
        self['output.k'] = self['input.[k]']
        return Node._check_ready(self)
        ## I think this is a source of bug and useless
        #for i in self._nodes.values():
        #    if not i.is_ready:
        #        return False

    def get_dock(self, key):
        '''Gets the object of a dock with a given key.'''
        try:
            return Node.get_dock(self, key)
        except:
            try:
                return self._subctx.get_dock(key)
            except:
                raise ValueError('{} key doesn\'t correspond to any '
                                 'dock in this supernode.'.format(key))

    def items(self):
        items = Node.vars(self)

        for i,v in self._subctx.vars():
            items[i] = v

        return items

    def __getitem__(self, key):
        try:
            return Node.__getitem__(self, key)
        except:
            try:
                return self._subctx[key]
            except:
                raise ValueError('{} key doesn\'t correspond to any property '
                                 'or node from this supernode.'.format(key))

    def __setitem__(self, key, value):
        try:
            Node.__setitem__(self, key, value)
        except:
            try:
                self._subctx[key] = value
            except:
                raise ValueError('{} key doesn\'t correspond to any property '
                                 'or node from this supernode.'.format(key))


################################################################################


class ReadCCDData(Node):
    '''Read a bunch of fits files.'''
    def setup(self):
        self.add_input('file_name', StringInput)
        self.add_input('hdu_number', IntInput, 0)
        self.add_input('path', StringInput, '')
        self.add_output('image')

    def run(self):
        fname = path.join(self['input.path'], self['input.file_name'])
        hdu_n = self['input.hdu_number']
        self['output.image'] = CCDData.read(fname, hdu=hdu_n, unit='adu')
        self.log.info('File {}[{}] loaded.'.format(fname, hdu_n))


class WriteCCDData(Node):
    '''Write a CCDData instance to a fits file'''
    def setup(self):
        self.add_input('image')
        self.add_input('file_name', StringInput)
        #self.add_input('save_jpeg', 'bool', False)

    def run(self):
        self.log.info('Saving data to file {}'.format(self['input.file_name']))
        image = self['input.image']
        image.to_hdu().writeto(self['input.file_name'], overwrite=True)


class CombineImages(Node):
    def setup(self):
        self.add_input('image_list')
        self.add_input('method', StringInput, 'median')
        self.add_input('weights', default_value=None, optional=True)
        self.add_input('normalize', BoolInput, False)
        self.add_output('image')
        self._properties['mem_limit'] = 1e7
        self._properties['scale'] = False
        #TODO: include all arguments

    def run(self):
        method = self['input.method']
        images = self['input.image_list']
        weights = self['input.weights']
        mem_limit = self['property.mem_limit']
        normalize = self['input.normalize']
        self.log.debug('Combining {} images with {} method.'
                       ''.format(len(images), method))

        scaling = lambda arr: 1/np.nanmean(arr)
        if self['property.scale']:
            scale = scaling
        else:
            scale = None

        if method in ('median', 'average', 'sum'):
            result = combine(images, scale=scale, method=method,
                             mem_limit=mem_limit)
        else:
            raise ValueError('Method {} unrecognized.'.format(method))

        if normalize:
            result.data *= scaling(result.data)

        result.meta['hierarch ccdproc combiner n_images'] = len(images)
        result.meta['hierarch ccdproc combiner method'] = method
        result.meta['hierarch ccdproc combiner normalize'] = normalize
        self['output.image'] = result

class BinAndTrim(Node):
    '''Pre-process functions for pre-binning and trimming in images.'''
    def setup(self):
        self.add_input('image')
        self.add_output('image')
        self.add_input('trim_rect', optional=True)
        self.add_input('vbin', IntInput, 1, optional=True)
        self.add_input('hbin', IntInput, 1, optional=True)

    def run(self):
        image = self['input.image']
        if self.get_dock('input.trim_rect').have_value:
            self.log.dedug('Trimming image with rect: {}'
                           ''.format(self['input.trim_rect']))
            image = trim_image(image, fits_section=self['input.trim_rect'],
                               add_keyword = {'hierarch ccdproc trimmed_rect' :
                                              self['input.trim_rect']})
        self.log.debug('Binning the data with block size of ({}, {})'
                       ''.format(self['input.vbin'], self['input.hbin']))
        image = block_reduce(image, (self['input.vbin'], self['input.hbin']))
        image.meta['hierarch ccdproc binning'] = '{}x{}'.format(self['input.vbin'],
                                                                self['input.hbin'])
        self['output.image'] = image


class LACosmic(Node):
    '''Remove cosmic rays with L.A. Cosmic algorith from astroscrappy.'''
    def setup(self):
        self.add_input('image')
        self.add_output('image')
        self.add_output('cosmics_map')

    def run(self):
        im = self['input.image']
        cmx, imc = detect_cosmics(im.data)
        self['output.cosmics_map'] = cmx
        self['output.image'] = CCDData(imc, meta=im.meta, unit=im.unit)


class CCDProcess(Node):
    '''Node for bias/flat corrections.'''
    def setup(self):
        self.add_input('image')
        self.add_output('image')
        self.add_input('gain', FloatInput, optional=True)
        self.add_input('oscan', optional=True)
        self.add_input('bias', optional = True)
        self.add_input('dark', optional = True)
        self.add_input('flat', optional = True)
        self.add_input('badpix', optional = True)
        self._properties['kwargs'] = {}

    def run(self):
        keywords = {}
        kwargs = {}
        kwargs.update(self._properties['kwargs'])
        if self.get_dock('input.gain').have_value:
            kwargs['gain'] = self['input.dock']
            keywords['hierarch ccdproc gain'] = kwargs['gain']

        if self.get_dock('input.oscan').have_value:
            kwargs['oscan'] = self['input.oscan']
            keywords['hierarch ccdproc oscan'] = kwargs['oscan']

        if self.get_dock('input.bias').have_value:
            inp = self['input.bias']
            if isinstance(inp, string_types):
                inp = CCDData.read(inp, hdu=0, unit='adu')
            kwargs['master_bias'] = inp
            keywords['hierarch ccdproc bias'] = kwargs['master_bias']

        if self.get_dock('input.dark').have_value:
            inp = self['input.dark']
            if isinstance(inp, string_types):
                inp = CCDData.read(inp, hdu=0, unit='adu')
            kwargs['dark_frame'] = inp
            keywords['hierarch ccdproc dark'] = kwargs['dark_frame']

        if self.get_dock('input.flat').have_value:
            inp = self['input.flat']
            if isinstance(inp, string_types):
                inp = CCDData.read(inp, hdu=0, unit='adu')
            kwargs['master_flat'] = inp
            keywords['hierarch ccdproc flat'] = kwargs['master_flat']

        if self.get_dock('input.badpix').have_value:
            inp = self['input.badpix']
            if isinstance(inp, string_types):
                inp = CCDData.read(inp, hdu=0, unit='adu')
            kwargs['bad_pixel_mask'] = inp
            keywords['hierarch ccdproc badpix'] = kwargs['bad_pixel_mask']

        image = self['input.image']
        image = ccd_process(image, add_keyword=keywords, **kwargs)
        self['output.image'] = image
