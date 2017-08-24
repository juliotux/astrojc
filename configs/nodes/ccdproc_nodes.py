import os
from os import path

import numpy as np
from copy import copy
import json
from collections import OrderedDict
from six import string_types

from astropy.io import fits
from astropy import units as u
from astroscrappy import detect_cosmics
from ccdproc import CCDData, block_reduce, trim_image, combine, ccd_process
from astrojc.pipeline import *


__all__ = ['ReadCCDData', 'WriteCCDData', 'CombineImages', 'BinAndTrim',
           'LACosmic', 'CCDProcess']


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
        self.add_input('path', StringInput, optional=True)
        #self.add_input('save_jpeg', 'bool', False)

    def run(self):
        image = self['input.image']
        if self.get_dock('input.path').have_value:
            fpath = self['input.path']
        else:
            fpath = path.dirname(self['input.file_name'])
        name = path.join(fpath, path.basename(self['input.file_name']))
        self.log.info('Saving data to file {}'.format(name))
        #image.to_hdu().writeto(name, overwrite=True)


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
