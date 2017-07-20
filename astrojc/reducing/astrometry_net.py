import re
import tempfile
import functools
import os
import shutil
import subprocess
import textwrap

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy import units
import numpy as np

from ..logging import log

ASTROMETRY_COMMAND = '/usr/bin/solve-field'

class AstrometryNetUnsolvedField(subprocess.CalledProcessError):
    """ Raised if Astrometry.net could not solve the field """
    def __init__(self, path):
        self.path = path
    def __str__(self):
        return "{0}: could not solve field".format(self.path)

def create_xyls(fname, x, y, flux, imagew, imageh, header=None, dtype='f8'):
    '''
    Create and save the xyls file to run in astrometry.net

    Parameters:
    -----------
    fname : str
        The path to save the .xyls file.
    x : array_like
        X coordinates of the sources.
    y : array_like
        Y coordinates of the sources.
    flux : array_like
        Estimated flux of the sources.
    imagew : float or int
        Width of the original image. `IMAGEW` header field of .xyls
    imageh : float or int
        Height of the original image. `IMAGEH` header field of .xyls
    header : dict_like, optional
        Header of the orginal image. If None, no header will be added to the
        PrimaryHDU of the .xyls
    dtype : `numpy.dtype`, optional
        Data type of the fields. Default: 'f8'
    '''
    head = {'IMAGEW':imagew,
            'IMAGEH':imageh}
    xyls = np.array(list(zip(x, y, flux)), np.dtype([('x', dtype), ('y', dtype), ('flux', dtype)]))
    f = fits.HDUList([fits.PrimaryHDU(header=header),
                      fits.BinTableHDU(xyls,header=fits.Header(head))])
    log.info('Saving .xyls to ' + fname)
    f.writeto(fname)

def solve_field(path, stdout=None, stderr=None, **options):
    '''
    Function to call `solve_field` astrometry.net command, that reads the output
    and returns the solved WCS and the outrput_dir.
    '''
    #TODO: Better handling the solver output (stdout).
    #TODO: Better handling the errors.
    basename = os.path.basename(path)
    root, ext = os.path.splitext(basename)
    # Place all output files in this directory
    kwargs = dict(prefix = root + '_', suffix = '_astrometry.net')
    output_dir = tempfile.mkdtemp(**kwargs)

    solved_file = os.path.join(output_dir, root + '.solved')

    args = [ASTROMETRY_COMMAND, path, '--dir', output_dir]

    for key, value in options.items():
        ndashes = 1 if len(key) == 1 else 2
        args.append("{0}{1}".format(ndashes * '-', key))
        if value is not None:
            args.append(str(value))

    try:
        log.info('runing: ' + str(args))
        subprocess.check_call(args, stdout=stdout, stderr=stderr)

        # .solved file must exist and contain a binary one
        with open(solved_file, 'rb') as fd:
            if ord(fd.read()) != 1:
                raise AstrometryNetUnsolvedField(path)
        solved_wcs_file = os.path.join(output_dir, root + '.wcs')
        log.info('Loading solved header from %s' % solved_wcs_file)
        solved_header = fits.getheader(solved_wcs_file, 0)
        return solved_header, output_dir

    except subprocess.CalledProcessError as e:
        #shutil.rmtree(output_dir)
        raise e
    # If .solved file doesn't exist or contain one
    except (IOError, AstrometryNetUnsolvedField):
        #shutil.rmtree(output_dir)
        raise AstrometryNetUnsolvedField(path)

def _get_coordinates(header, rak, deck):
    '''
    Guess the center coordinates of field based in header keys.
    '''
    #The universal function
    ra  = str(header[rak]).replace(',','.')
    dec = str(header[deck]).replace(',','.')
    coords = functools.partial(SkyCoord, ra, dec)

    regexp = "\d{1,3}\.?\d"
    match_degrees = functools.partial(re.match, regexp)
    if match_degrees(ra) and match_degrees(dec):
        return coords(unit=(units.deg, units.deg))

    # Assume (at least for now) that it's in sexagesimal
    return coords(unit=(units.hourangle, units.deg))

def solve(path, output_file=None, ra=None, dec=None, pltscl=None,
          rak=None, deck=None, pltk=None, radius = 1,
          scamp_basename=None):
    '''
    Solve a image or .xyls file located in `path`, with solving optimized with
    several parameters, like plate scale, center coordinates ou header keys.

    Parameters:
    -----------
    path : str
        Path of the image or .xyls files to solve.
    output_file : str, optional
        Path to store the result file (`new-fits` argument of solve-field).
    ra : str or float, optional
        Aproximated Right Ascension of the center of the field in degrees.
    dec : str or float, optional
        Aproximated Declination of the center of the field in degrees.
    pltscl: str or float, optional
        Aproximated plate scale of the field, in arcsec. We assume it has 10% of
        accuracy.
    rak : str, optional
        Header key for extract the right ascension of the center of the field
        from header.
    deck : str, optional
        Header key for extract the declination of the center of the field
        from header.
    pltk : str, optional
        Header key for extract the plate scale of the center of the field
        from header.
    scamp_basename : str, optional
        Basename to prepare the files for posterior `scamp` solving of the
        filed. It will passe the parameters `scamp-ref`, `scamp-config` and
        `scamp` to `solve-filed`, generating the files `scamp_basename.ref`
        `scamp_basename.scamp` and `scamp_basename.cat`.

    Returns:
    --------
    header : `astropy.wcs.WCS`
        A header containing the solved wcs of the field and another information
        from astrometry.net.

    '''
    # Path to the temporary FITS file containing the WCS header
    basename = os.path.basename(path)
    root, ext = os.path.splitext(basename)
    kwargs = dict(prefix = '{0}_astrometry_'.format(root), suffix = ext)
    with tempfile.NamedTemporaryFile(**kwargs) as fd:
        output_path = fd.name

    options = {
        'no-plot' : None,
        #'no-fits2fits' : None,
        'overwrite' : None,
        'depth' : "40,60,80,120"
        }

    if scamp_basename is not None:
        options['scamp-ref'] = scamp_basename + ".ref"
        options['scamp-config'] = scamp_basename + ".scamp"
        options['scamp'] = scamp_basename + ".cat"

    if output_file is not None:
        options['new-fits'] = output_file

    hdulist = None

    if None not in (ra, dec):
        log.info("Usign given RA and DEC coordinates: %s    %s" % (str(ra), str(dec)))
        options['ra'] = ra
        options['dec'] = dec
        options['radius'] = radius
    elif None not in (rak, deck):
        log.info("Figuring out field center coordinates")
        hdulist = fits.open(path)
        header = hdulist[0].header
        try:
            coords = _get_coordinates(header, rak, deck)
        except KeyError as e:
            log.warn("Cannot understand coordinates in FITS header")
            log.warn("Astrometry.net will try to solve the image blindly")
        else:
            options['ra']  = coords.ra.degree
            options['dec'] = coords.dec.degree
            options['radius'] = radius

    if pltscl is not None:
        log.info("Usign given plate scale: %s" % str(pltscl))
        options['scale-units'] = 'arcsecperpix'
        options['scale-high'] = 1.1*pltscl
        options['scale-low'] = 0.9*pltscl
    elif not pltk is None:
        log.info("Figuring out the plate scale from FITS header")
        if hdulist is None:
            hdulist = fits.open(path)
            header = hdulist[0].header
        try:
            pltscl = float(str(header[pltk]).replace(',','.')) #fix for OPD
            options['scale-units'] = 'arcsecperpix'
            options['scale-high'] = 1.1*pltscl
            options['scale-low'] = 0.9*pltscl
        except KeyError:
            log.warn("Cannot understand plate scale in FITS header")

    with open(os.devnull, 'wb') as fd:
        log.info("Running {0}".format(ASTROMETRY_COMMAND))
        solved_header, output_dir = solve_field(path, stdout=fd, stderr=fd, **options)

    log.info("Removing working directory")
    shutil.rmtree(output_dir)
    return solved_header
