'''Wrapper from ccdproc package.'''
#TODO: change this functions to self written functions, with better handling of memory and multiprocessing.

import ccdproc
from ccdproc import CCDData
from ccdproc.combiner import combine

def process_image(*args, **kwargs):
    '''ccdproc.ccd_process(ccd, oscan=None, trim=None, error=False, master_bias=None,
                           dark_frame=None, master_flat=None, bad_pixel_mask=None, gain=None,
                           readnoise=None, oscan_median=True, oscan_model=None, min_value=None,
                           dark_exposure=None, data_exposure=None, exposure_key=None,
                           exposure_unit=None, dark_scale=False, add_keyword=True)'''
    return ccdproc.ccd_process(*args, **kwargs)

def combine_bias(ccd_list, combine_method='median', **reductwargs):
    for i in range(len(ccd_list)):
        ccd_list[i] = process_image(ccd_list[i], **reductwargs)
    return combine(ccd_list, method=combine_method, mem_limit=mem_limit)

def combine_flat(ccd_list, master_bias, combine_method='median', **reductwargs):
    for i in range(len(ccd_list)):
        ccd_list[i] = process_image(ccd_list[i], master_bias=master_bias, **reductwargs)
        ccd_list[i].data = ccd_list[i].data/np.median(ccd_list[i].data)
    return combine(ccd_list, method=combine_method, mem_limit=mem_limit)

def read_fits(fname):
    return CCDData.read(fname, unit='adu')

def save_fits(ccd, fname):
    ccd.write(fname)
