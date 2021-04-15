"""

Main 3d rank filter module
--------------------------

The package is based on the 2d skimage.filters.rank filter module.

These filters compute the local histogram at each pixel, using a sliding window
similar to the method described in [1]_. A histogram is built using a moving
window in order to limit redundant computation. The moving window follows a
snake-like path:

...------------------------\
/--------------------------/
\--------------------------...

The local histogram is updated at each pixel as the structuring element window
moves by, i.e. only those pixels entering and leaving the structuring element
update the local histogram. The histogram size is 8-bit (256 bins) for 8-bit
images and 2- to 16-bit for 16-bit images depending on the maximum value of the
image.

The filter is applied up to the image border, the neighborhood used is
adjusted accordingly. The user may provide a mask image (same size as input
image) where non zero values are the part of the image participating in the
histogram computation. By default the entire image is filtered.

This implementation outperforms grey.dilation for large structuring elements.

Input image can be 8-bit or 16-bit, for 16-bit input images, the number of
histogram bins is determined from the maximum value present in the image.

Result image is 8-/16-bit or double with respect to the input image and the
rank filter operation.


References
----------

.. [1] Huang, T. ,Yang, G. ;  Tang, G.. "A fast two-dimensional
       median filtering algorithm", IEEE Transactions on Acoustics, Speech and
       Signal Processing, Feb 1979. Volume: 27 , Issue: 1, Page(s): 13 - 18.

"""
__note__ = "Code adpated to 3D images from skimage.filters.rank by Christoph Kirst, The Rockefeller University, New York City, 2017."


import functools
import numpy as np

from warnings import warn

from scipy import ndimage as ndi

from skimage import img_as_ubyte


import pyximport;
pyximport.install(setup_args={"include_dirs":np.get_include()}, reload_support=True)

from . import generic_cy


__all__ = ['autolevel', 'bottomhat', 'equalize', 'gradient', 'maximum', 'mean',
           'geometric_mean', 'subtract_mean', 'median', 'minimum', 'modal',
           'enhance_contrast', 'pop', 'threshold', 'tophat', 'noise_filter',
           'entropy', 'otsu']


def _handle_input(image, selem, out, mask, out_dtype=None, pixel_size=1):

    assert len(image.shape) == 3
    
    if image.dtype not in (np.uint8, np.uint16):
        image = img_as_ubyte(image)

    selem = np.ascontiguousarray(img_as_ubyte(selem > 0))
    
    #TODO: fix fortran vs c order here
    image = np.ascontiguousarray(image)

    #if mask is None:
    #    mask = np.ones(image.shape, dtype=np.uint8)
    #else:
    #    mask = img_as_ubyte(mask)
    #    mask = np.ascontiguousarray(mask)
    mask = img_as_ubyte(np.zeros((1,1,1), dtype = bool));

    if image is out:
        raise NotImplementedError("Cannot perform rank operation in place.")

    if out is None:
        if out_dtype is None:
            out_dtype = image.dtype
        out = np.empty(image.shape+(pixel_size,), dtype=out_dtype)
    else:
        if len(out.shape) == 3:
            out = out.reshape(out.shape+(pixel_size,))

    #is_8bit = image.dtype in (np.uint8, np.int8)

    #if is_8bit:
    max_bin = 255
    #else:
    #    max_bin = max(4, image.max())

    bitdepth = int(np.log2(max_bin))
    if bitdepth > 10:
         warn("Bitdepth of %d may result in bad rank filter performance due to large number of bins." % bitdepth);

    return image, selem, out, mask, max_bin


def _apply_scalar_per_pixel(func, image, selem, out, mask, shift_x, shift_y, shift_z,
                            out_dtype=None):

    image, selem, out, mask, max_bin = _handle_input(image, selem, out, mask,
                                                     out_dtype)

    func(image, selem, shift_x=shift_x, shift_y=shift_y, shift_z=shift_z, mask=mask,
         out=out, max_bin=max_bin)

    return out.reshape(out.shape[:3])


def _apply_vector_per_pixel(func, image, selem, out, mask, shift_x, shift_y, shift_z,
                            out_dtype=None, pixel_size=1):

    image, selem, out, mask, max_bin = _handle_input(image, selem, out, mask,
                                                     out_dtype,
                                                     pixel_size=pixel_size)

    func(image, selem, shift_x=shift_x, shift_y=shift_y, shift_z=shift_z, mask=mask,
         out=out, max_bin=max_bin)

    return out


def _default_selem(func):
    """Decorator to add a default structuring element to morphology functions.

    Parameters
    ----------
    func : function
        A morphology function such as erosion, dilation, opening, closing,
        white_tophat, or black_tophat.

    Returns
    -------
    func_out : function
        The function, using a default structuring element of same dimension
        as the input image with connectivity 1.
    """
    @functools.wraps(func)
    def func_out(image, selem=None, *args, **kwargs):
        if selem is None:
            selem = ndi.generate_binary_structure(image.ndim, image.ndim)
        return func(image, selem=selem, *args, **kwargs)

    return func_out


def autolevel(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Auto-level image using local histogram.

    This filter locally stretches the histogram of greyvalues to cover the
    entire range of values from "white" to "black".

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.

    """

    return _apply_scalar_per_pixel(generic_cy._autolevel, image, selem,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


def bottomhat(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Local bottom-hat of an image.

    This filter computes the morphological closing of the image and then
    subtracts the result from the original image.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : 3-D array
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.
    """

    return _apply_scalar_per_pixel(generic_cy._bottomhat, image, selem,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


def equalize(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Equalize image using local histogram.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.
    """

    return _apply_scalar_per_pixel(generic_cy._equalize, image, selem,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


def gradient(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Return local gradient of an image (i.e. local maximum - local minimum).

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.
    """

    return _apply_scalar_per_pixel(generic_cy._gradient, image, selem,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


def maximum(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Return local maximum of an image.

    Parameters
    ----------
    image :3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.

    See also
    --------
    skimage.morphology.dilation

    Notes
    -----
    The lower algorithm complexity makes `skimage.filters.rank.maximum`
    more efficient for larger images and structuring elements.
    """

    return _apply_scalar_per_pixel(generic_cy._maximum, image, selem,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


def mean(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Return local mean of an image.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.
    """

    return _apply_scalar_per_pixel(generic_cy._mean, image, selem, out=out,
                                   mask=mask, shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


def geometric_mean(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Return local geometric mean of an image.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.

    References
    ----------
    .. [1] Gonzalez, R. C. and Wood, R. E. "Digital Image Processing (3rd Edition)."
           Prentice-Hall Inc, 2006.
    """

    return _apply_scalar_per_pixel(generic_cy._geometric_mean, image, selem, out=out,
                                   mask=mask, shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


def subtract_mean(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Return image subtracted from its local mean.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.
    """

    return _apply_scalar_per_pixel(generic_cy._subtract_mean, image, selem,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


@_default_selem
def median(image, selem=None, out=None, mask=None,  shift_x=False, shift_y=False, shift_z=False):
    """Return local median of an image.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array, optional
        The neighborhood expressed as a 3-D array of 1's and 0's. If None, a
        full square of size 3 is used.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.
    """

    return _apply_scalar_per_pixel(generic_cy._median, image, selem,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


def minimum(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Return local minimum of an image.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.

    See also
    --------
    skimage.morphology.erosion

    Notes
    -----
    The lower algorithm complexity makes `skimage.filters.rank.minimum` more
    efficient for larger images and structuring elements.
    """

    return _apply_scalar_per_pixel(generic_cy._minimum, image, selem,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


def modal(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Return local mode of an image.

    The mode is the value that appears most often in the local histogram.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.
    """

    return _apply_scalar_per_pixel(generic_cy._modal, image, selem,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


def enhance_contrast(image, selem, out=None, mask=None, shift_x=False,
                     shift_y=False, shift_z=False):
    """Enhance contrast of an image.

    This replaces each pixel by the local maximum if the pixel greyvalue is
    closer to the local maximum than the local minimum. Otherwise it is
    replaced by the local minimum.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        The result of the local enhance_contrast.
    """

    return _apply_scalar_per_pixel(generic_cy._enhance_contrast, image, selem,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


def pop(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Return the local number (population) of pixels.

    The number of pixels is defined as the number of pixels which are included
    in the structuring element and the mask.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.
    """

    return _apply_scalar_per_pixel(generic_cy._pop, image, selem, out=out,
                                   mask=mask, shift_x=shift_x,
                                   shift_y=shift_y, shift_z=shift_z)


def sum(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Return the local sum of pixels.

    Note that the sum may overflow depending on the data type of the input
    array.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.
    """

    return _apply_scalar_per_pixel(generic_cy._sum, image, selem, out=out,
                                   mask=mask, shift_x=shift_x,
                                   shift_y=shift_y, shift_z=shift_z)


def threshold(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Local threshold of an image.

    The resulting binary mask is True if the greyvalue of the center pixel is
    greater than the local mean.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.
    """

    return _apply_scalar_per_pixel(generic_cy._threshold, image, selem,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


def tophat(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Local top-hat of an image.

    This filter computes the morphological opening of the image and then
    subtracts the result from the original image.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.
    """

    return _apply_scalar_per_pixel(generic_cy._tophat, image, selem,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


def noise_filter(image, selem, out=None, mask=None, shift_x=False,
                 shift_y=False, shift_z=False):
    """Noise feature.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    References
    ----------
    .. [1] N. Hashimoto et al. Referenceless image quality evaluation
                     for whole slide imaging. J Pathol Inform 2012;3:9.

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.
    """

    # ensure that the central pixel in the structuring element is empty
    centre_r = int(selem.shape[0] / 2) + shift_y
    centre_c = int(selem.shape[1] / 2) + shift_x
    # make a local copy
    selem_cpy = selem.copy()
    selem_cpy[centre_r, centre_c] = 0

    return _apply_scalar_per_pixel(generic_cy._noise_filter, image, selem_cpy,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z)


def entropy(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Local entropy.

    The entropy is computed using base 2 logarithm i.e. the filter returns the
    minimum number of bits needed to encode the local greylevel
    distribution.

    Parameters
    ----------
    image : 3-D array (uint8, uint16)
        Input image.
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : 3-D array (same dtype as input)
        If None, a new array is allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : ndarray (double)
        Output image.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Entropy_(information_theory)
    """

    return _apply_scalar_per_pixel(generic_cy._entropy, image, selem,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z,
                                   out_dtype=np.double)


def otsu(image, selem, out=None, mask=None, shift_x=False, shift_y=False, shift_z=False):
    """Local Otsu's threshold value for each pixel.

    Parameters
    ----------
    image : ndarray
        Image array (uint8 array).
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : ndarray
        If None, a new array will be allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).

    Returns
    -------
    out : 3-D array (same dtype as input image)
        Output image.

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Otsu's_method
    """

    return _apply_scalar_per_pixel(generic_cy._otsu, image, selem, out=out,
                                   mask=mask, shift_x=shift_x,
                                   shift_y=shift_y, shift_z=shift_z)


def windowed_histogram(image, selem, out=None, mask=None,
                       shift_x=False, shift_y=False, shift_z=False, n_bins=None):
    """Normalized sliding window histogram

    Parameters
    ----------
    image : ndarray
        Image array (uint8 array).
    selem : 3-D array
        The neighborhood expressed as a 3-D array of 1's and 0's.
    out : ndarray
        If None, a new array will be allocated.
    mask : ndarray
        Mask array that defines (>0) area of the image included in the local
        neighborhood. If None, the complete image is used (default).
    shift_x, shift_y, shift_z : int
        Offset added to the structuring element center point. Shift is bounded
        to the structuring element sizes (center must be inside the given
        structuring element).
    n_bins : int or None
        The number of histogram bins. Will default to ``image.max() + 1``
        if None is passed.

    Returns
    -------
    out : 3-D array with float dtype of dimensions (H,W,N), where (H,W) are
        the dimensions of the input image and N is n_bins or
        ``image.max() + 1`` if no value is provided as a parameter.
        Effectively, each pixel is a N-D feature vector that is the histogram.
        The sum of the elements in the feature vector will be 1, unless no
        pixels in the window were covered by both selem and mask, in which
        case all elements will be 0.

    """

    if n_bins is None:
        n_bins = int(image.max()) + 1

    return _apply_vector_per_pixel(generic_cy._windowed_hist, image, selem,
                                   out=out, mask=mask,
                                   shift_x=shift_x, shift_y=shift_y, shift_z=shift_z,
                                   out_dtype=np.double,
                                   pixel_size=n_bins)
