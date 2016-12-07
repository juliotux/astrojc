#line identification i/o library

import pandas as pd
import numpy as np

__all__ = ['lineID','specID']

greek = lambda greek_code : chr(greek_code).encode('utf-8')

latex_dict = {'$':'',
              '\\,':' ',
              '{\\sc i}':'I',
              '{\\sc ii}':'II',
              '{\\sc iii}':'III',
              '{\\sc iv}':'IV',
              '{\\sc v}':'V',
              '{\\sc vi}':'VI',
              '{\\sc vii}':'VII',
              '{\\sc viii}':'VIII',
              '{\\sc ix}':'IX',
              '{\\sc x}':'X',
              '\\alpha':greek(0x3b1),
              '\\beta':greek(0x3b2),
              '\\gamma':greek(0x3b3),
              '\\delta':greek(0x3b4),
              '\\epsilon':greek(0x3b5),
              '\\zeta':greek(0x3b6),
              '\\eta':greek(0x3b7),
              '\\theta':greek(0x3b8),
              '\\iota':greek(0x3b9),
              '\\kappa':greek(0x3ba),
              '\\lambda':greek(0x3bb),
              '\\mu':greek(0x3bc),
              '\\nu':greek(0x3bd),
              '\\xi':greek(0x3be),
              '\\omicron':greek(0x3bf),
              '\\pi':greek(0x3c0),
              '\\rho':greek(0x3c1),
              '\\sigma':greek(0x3c3),
              '\\tau':greek(0x3c4),
              '\\upsilon':greek(0x3c5),
              '\\phi':greek(0x3c6),
              '\\chi':greek(0x3c7),
              '\\psi':greek(0x3c8),
              '\\omega':greek(0x3c9)}

class lineID():
    def __init__(self, wavelength, line_id, mult=None, weq=None):
        '''
        Class to store a single line identification.
        '''
        self.wavelength = wavelength
        self.line_id = line_id
        self.mult = mult
        self.weq = weq

    def fix_latex(self):
        for i in latex_dict.keys():
            self.line_id = self.line_id.replace(i, latex_dict[i])

class specID():
    def __init__(self, fname=None, cmfgen=False, latex_fix=False):
        '''
        Class to stores the identifications of a huge number of lines in one
        spectrum.
        '''
        self.ids = []

        if not fname is None and not cmfgen:
            self.from_file(fname)
            if latex_fix:
                self.fix_latex()

    def from_file(self, fname):
        '''
        Reads the identification table from a file in the format:

        # comment
        lambda  id  mult

        The separation must be a tab space.
        '''
        mf = pd.read_csv(fname, sep='\t', header=None, names=['wave','id','m'])
        self.ids = [None] * len(mf.values)
        
        for i in range(len(mf.values)):
            if mf.values[i][2] == '-':
                self.ids[i] = lineID(mf.values[i][0],
                                     mf.values[i][1],
                                     None)
            else:
                self.ids[i] = lineID(mf.values[i][0],
                                     mf.values[i][1],
                                     mf.values[i][2])

    def add_id(self, wave, line_id, mult=None, weq=None):
        '''
        Add a single identification to table.
        '''
        self.ids.append(lineID(wave, line_id, mult, weq))

    
    def fix_rv(self, rv):
        '''
        Correct the wavelength with radial velocity (in km/s).
        '''
        for i in self.ids:
            i.wavelength = i.wavelength*(1-(rv/3E5))

    def fix_latex(self):
        for i in self.ids:
            i.fix_latex()

    wavelength = lambda self : np.array([i.wavelength for i in self.ids])
    line_id = lambda self : np.array([i.line_id for i in self.ids])
    mult = lambda self : np.array([i.mult for i in self.ids])
    weq = lambda self : np.array([i.weq for i in self.ids])