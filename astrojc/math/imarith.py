import six
import copy
import numpy as np
from astropy.io import fits

from ..logging import log as logger

__all__ = ['check_hdu', 'extrema_clip', 'sigma_clip', 'minmax_clip',
           'imcombine', 'imarith']

imhdus = (fits.ImageHDU, fits.PrimaryHDU, fits.CompImageHDU,
          fits.StreamingHDU)


class IncompatibleHeadersError(ValueError):
    """When 2 image header are not compatible"""


def check_header_keys(image1, image2, keywords=[]):
    """Compare some header keys from 2 images and check if the have equal
    values."""
    image1 = check_hdu(image1)
    image2 = check_hdu(image2)
    for i in keywords:
        if i in image1.header.keys():
            if i in image2.header.keys():
                v1 = image1.header[i]
                v2 = image2.header[i]
                if v1 != v2:
                    raise IncompatibleHeadersError('Keyword {} have different '
                                                   'values for images 1 and 2:'
                                                   '\n{}\n{}'.format(i,
                                                                     v1, v2))
            else:
                raise IncompatibleHeadersError('Image 2 do not have Keyword {}'
                                               .format(i))
        else:
            raise IncompatibleHeadersError('Image 1 do not have Keyword {}'
                                           .format(i))
    return True


def check_hdu(data):
    """Check if a data is a valid ImageHDU type and convert it."""
    if not isinstance(data, imhdus):
        if isinstance(data, fits.HDUList):
            data = data[0]
        elif isinstance(data, six.string_types):
            data = fits.open(data)[0]
        elif isinstance(data, np.ndarray):
            data = fits.PrimaryHDU(data)
        else:
            raise ValueError('The given data is not a valid CCDData type.')
    return data


def extrema_clip(data, nlow=1, nhigh=1, axis=None):
    """Creates a mask clipping nlow and nhigh pixels."""
    # TODO: Not working. Use @vectorize?
    data = np.array(data)
    mask = np.zeros(data.shape, dtype=bool)

    if nlow is None:
        nlow = 0
    if nhigh is None:
        nhigh = 0

    if axis is None:
        mask = mask.ravel()
        argsorted = np.argsort(data, axis=axis)
        for i in range(-1*nhigh, nlow):
            mask[np.where(argsorted[i])] = True
        mask.reshape(data.shape)
    else:
        argsorted = np.argsort(data, axis=axis)
        mg = np.mgrid[[slice(ndim)
                       for i, ndim in enumerate(data.shape) if i != axis]]
        for i in range(-nhigh, nlow):
            where = tuple([argsorted[i, :, :].ravel()] +
                          [j.ravel() for j in mg])
            mask[where] = True
    return mask


def sigma_clip(data, sigma_clip_low=3, sigma_clip_high=3, func=np.nanmedian,
               dev_func=np.nanstd, axis=None):
    """Creates a mask of the sigma clipped pixels."""
    data = np.array(data)
    mask = np.zeros(data.shape, dtype=bool)

    base = func(data, axis=axis)
    dev = dev_func(data, axis=axis)

    if sigma_clip_low is not None:
        mask[np.where(base-data > abs(sigma_clip_low*dev))] = True
    if sigma_clip_high is not None:
        mask[np.where(data-base > abs(sigma_clip_high*dev))] = True

    return mask


def minmax_clip(data, min_clip=None, max_clip=None, axis=None):
    """Creates a mask of pixels clipped between min_clip and max_clip vals."""
    data = np.array(data)
    mask = np.zeros(data.shape, dtype=bool)

    if min_clip is not None:
        mask[np.where(data < min_clip)] = True

    if max_clip is not None:
        mask[np.where(data > max_clip)] = True

    return mask


_comb_funcs = {'average': np.nanmean,
               'median': np.nanmedian,
               'sum': np.nansum}


def imcombine(image_list, output_file=None, method='average', weights=None,
              scale=None, mem_limit=1e8, reject=None,
              nlow=1, nhigh=1, min_clip=None, max_clip=None, sigma_clip_low=3,
              sigma_clip_high=3, dtype=None):
    """Combine a set of images like IRAF imcombine.

    Methods:
        - 'average' : np.nanmean
        - 'median' : np.nanmedian
        - 'sum' : np.sum

    reject are:
        - 'extrema' : nlow and nhigh extrema pixels will be rejected
        - 'minmax' : data will be clipped between min_clip and max_clip
        - 'sigmaclip' : sigma clip will be used
    It is possible to use a combination of them as a list. Like:
    ['sigma', 'minmax']

    scale can be a function applyed for each image, or a list of scale values,
    one per image.

    Keeps the header of the first image.

    Partially using ccdproc's combine function.
    """
    if not isinstance(image_list, (list, tuple)):
        if isinstance(image_list, np.ndarray):
            image_list = image_list.tolist()
        elif isinstance(image_list, six.string_types) and ',' in image_list:
            image_list = image_list.split(',')
        else:
            raise ValueError("Unrecognized image_list type.")

    # query the information of first image and check images
    n_image = len(image_list)
    im = check_hdu(image_list[0])
    if dtype is None:
        dtype = im.data.dtype  # first image data type
    shape = im.data.shape  # first image shape
    xs, ys = shape
    im.data = im.data.astype(dtype)
    nbytes = im.data.nbytes  # size of first image
    del im
    for i in image_list:
        if check_hdu(i).data.shape != shape:
            raise ValueError('Images with different shapes. Aborting.')

    # check the reject argument
    if reject is not None:
        if isinstance(reject, six.string_types):
            reject = [reject]
        elif not isinstance(reject, (tuple, list)):
            raise ValueError('reject type unrecognized: {}'.format(reject))

        for i in reject:
            if i not in ['sigmaclip', 'minmax', 'extrema']:
                raise ValueError('reject type unsuported: {}'.format(i))
    else:
        reject = []

    # prepare the scale of the images
    if callable(scale):
        scalevalues = []
        for i in image_list:
            scalevalues.append(scale(check_hdu(i).data))
        scale = np.array(scalevalues)
        del scalevalues
    if scale is not None:
        try:
            scale = np.array(scale)
        except Exception as e:
            raise ValueError('Could not read the scale values: {}'.format(e))
        if len(scale.shape) > 1:
            raise ValueError('scale must be 1D array, containing one scale'
                             ' per image.')
        if scale.shape[0] != n_image:
            raise ValueError('the number of scale values and the number of'
                             ' images are different.')

    nchunks = int((nbytes * len(image_list)) / mem_limit) + 1
    if nchunks > 1:
        logger.info('Splitting each image into {0} chunks to limit memory'
                    ' usage to {1} bytes.'.format(nchunks, mem_limit))
    xstep = max(1, int(xs/nchunks))
    if nchunks > xs:
        ystep = max(1, int(ys/(1 + nchunks - int(xs/xstep))))
    else:
        ystep = 1

    combined = np.zeros(shape, dtype=dtype)

    # process every chunk
    for x in range(0, xs, xstep):
        for y in range(0, ys, ystep):
            xend, yend = min(xs, x + xstep), min(ys, y + ystep)
            data_list = []
            for i in image_list:
                data_list.append(check_hdu(i).data[x:xend, y:yend])

            # apply the scale
            if scale is not None:
                for n in range(n_image):
                    data_list[n] = data_list[n]*scale[n]

            data_list = np.array(data_list)
            mask = np.zeros(data_list.shape, dtype=bool)
            # apply the clipping
            if 'extrema' in reject:
                # mask = mask | extrema_clip(data_list, nlow, nhigh, axis=0)
                raise ValueError('Extrema clipping not implemented yet')
            if 'minmax' in reject:
                mask = mask | minmax_clip(data_list, min_clip, max_clip)
            if 'sigmaclip' in reject:
                mask = mask | sigma_clip(data_list, sigma_clip_low,
                                         sigma_clip_high, axis=0)

            # TODO: make it less dumb
            try:
                data_list[np.where(mask)] = np.nan
            except Exception:
                pass

            combined[x:xend, y:yend] = _comb_funcs[method](data_list, axis=0)

    hdu = fits.PrimaryHDU(combined, header=check_hdu(image_list[0]).header)
    hdu.header['hierarch combined number'] = n_image
    hdu.header['hierarch combined method'] = method
    hdu.header['hierarch combined reject'] = ','.join(reject)

    if output_file is not None:
        logger.info('Combined image saved at {}'.format(output_file))
        hdu.writeto(output_file)

    return hdu


_arith_funcs = {'+': np.add,
                '-': np.subtract,
                '/': np.divide,
                '*': np.multiply,
                '%': np.remainder,
                '**': np.power}


def imarith(operand1, operand2, operation, inplace=False):
    """Simple arithmetic operations using fits hdus.

    Supported operations: '+', '-', '*', '/', '%', '**'

    Keeps the header of the first image.
    """
    logger.debug('Operation {} between {} and {}'.format(operation, operand1,
                                                         operand2))
    if not isinstance(operand1, imhdus) and not isinstance(operand2, imhdus):
        raise ValueError("Both operand1 and operand2 are not valid "
                         "fits images.")

    if operation not in _arith_funcs.keys():
        raise ValueError("Operation {} not supported.".format(operation))

    nhdu = None
    try:
        hdu = check_hdu(operand1)
        data1 = hdu.data
        if inplace:
            nhdu = hdu
        else:
            nhdu = fits.PrimaryHDU(hdu.data, hdu.header)
    except ValueError:
        data1 = operand1

    try:
        hdu = check_hdu(operand2)
        data2 = hdu.data
        if inplace and nhdu is None:
            nhdu = hdu
        elif nhdu is None:
            nhdu = fits.PrimaryHDU(hdu.data, hdu.header)
    except ValueError:
        data2 = operand2

    try:
        nhdu.data = _arith_funcs[operation](data1, data2)
    except Exception as e:
        raise ValueError('Could not process the operation {} between {} and {}'
                         .format(operation, data1, data2))

    return nhdu
