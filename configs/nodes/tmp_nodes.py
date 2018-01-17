from os import path

from astrojc.pipeline import *
from astrojc.io import mkdir_p


class TempPCCDPACKScript(Node):
    def setup(self):
        self.add_input('files')
        self.add_input('path')

    def run(self):
        print(self['input.path'], self['input.files'])
