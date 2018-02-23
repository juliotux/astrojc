import six
import numpy as np
from astropy.io import fits

from ..math.imarith import imarith, imcombine, check_hdu
from ..logging import log as logger


def subtract_bias(image, master_bias):
    """Subtract a master_bias frame from an image."""


def subtract_dark(image, dark_frame, dark_exposure=None, image_exposure=None,
                  exposure_key=None):
    """Subtract dark frame from an image, scaling by exposure."""


def subtract_overscan(image, overscan):
    """Subtract the overscan of one image."""


def trim_image(image, trim_section):
    """Trim an image to a section."""


def divide_flat(image, master_flat):
    """Divide a image by a maaster flat field frame."""


def gain_correct(image, gain=None, gain_key=None):
    """Process the gain correction of an image."""


def rebin(image, block_size, func=np.sum, readnoise_key=None):
    """Process rebinnig in one image. Like block_reduce."""


def lacosmic(image):
    """Remove cosmic rays with LAcosmic. From astroscrappy package."""


def process_image(image, master_biar=None, dark_frame=None, master_flat=None,
                  gain=None, gain_key=None, image_exposure=None, overscan=None,
                  dark_exposure=None, exposure_key=None):
    """Full process of an image."""
