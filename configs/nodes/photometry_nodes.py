import json

from astrojc.pipeline import *
from astrojc.reducing import photutils_wrapper as phtw


class PhotutilsFindSource(Node):
    def setup(self):
        self.add_input('image')
        self.add_input('backgroung')
        self.add_input('snr', FloatInput, 10, optional=True)
        self._properties['fwhm'] = 1.5

        self.add_output('x')
        self.add_output('y')
        self.add_output('flux')
        self.add_output('peak')
        self.add_output('sharpness')
        self.add_output('roundness')
        self.add_output('sky')

    def run(self):
        result = phtw.find_source(data, fwhm=self['property.fwhm'],
                                  bkg=self['input.background'],
                                  snr=self['input.snr'])

        for i in ['x', 'y', 'flux', 'peak', 'sky', 'sharpness', 'roundness']:
            self['output.{}'.format(i)] = result[i]

