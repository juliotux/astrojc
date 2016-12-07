#plot bokeh line ids

from line_ident import specID
import lineid_plot

import numpy as np
import pandas as pd

def plot_line_ids(ax, spec, spec_ids, min_weq=None, **kwargs):
    '''
    Plots the line IDs in an graph if the equivalent width of the line is bigger
    then min_weq.
    '''
    #TODO: implantar sistema de min_weq
    lineid_plot.plot_line_ids(spec.dispersion, spec.flux, spec_ids.wavelength(), spec_ids.line_id(), ax=ax)
