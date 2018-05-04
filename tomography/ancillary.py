
import os
import re
import numpy as np


def listfiles(path, pattern):
    """
    list files in a directory whose names match a regular expression

    Parameters
    ----------
    path: str
        the directory to be searched
    pattern: str
        the regular expression search pattern
    Returns
    -------
    list
        a list of absolute file names
    Example
    -------
    >>> listfiles('/path/to/some/data', 'file[0-9].tif')
    ['/path/to/some/data/file1.tif', '/path/to/some/data/file2.tif', '/path/to/some/data/file3.tif']
    """
    return sorted([os.path.join(path, x) for x in os.listdir(path) if re.search(pattern, x)])


def normalize(slice):
    """
    normalize a 1D array by its minimum and maximum values:

    .. math::
        y = \\frac{x-min(x)}{max(x)-min(x)}

    Parameters
    ----------
    slice: numpy.ndarray
        a 1d input array to be normalized

    Returns
    -------
    ndarray
        the normalized array
    """
    max = np.amax(slice)
    min = np.amin(slice)
    return np.divide((slice - min), (max - min))


def cbfi(slice, nTrack, height):
    """
    computation of capon beam forming inversion for a single pixel.
    This function is used internally by the core function :func:`~tomography.functions.capon_beam_forming_inversion`.

    Parameters
    ----------
    slice: numpy.ndarray
        an array containing the covariance matrix and wave number for a single pixel
    nTrack: int
        the number of original SLC files
    height: int
        the maximum inversion height

    Returns
    -------
    numpy.ndarray
        the tomographic result for one pixel
    """
    kz = np.transpose([slice[(nTrack ** 2):]])
    r0 = slice[0:(nTrack ** 2)].reshape((nTrack, nTrack))

    z_vector = np.matrix(np.arange(-height, height + 1, 1))

    # define the loading factor
    load_fac = (1 / 25.) * np.identity(nTrack)

    # define the steering matrix
    a_steer = np.exp(np.dot(complex(0, 1) * kz, z_vector))

    # define the numerator of the filter habf
    hnum = np.dot(np.linalg.inv(r0 + load_fac), a_steer)

    # define the denominator of the filter habf
    hden0 = np.diag(np.dot(np.conjugate(np.transpose(a_steer)), hnum))

    # replicate the diagonal by number of Tracks
    hden = np.array([hden0] * nTrack, dtype=np.complex64)

    h_abf = np.divide(hnum, hden)

    return np.diag(np.dot(np.dot(np.conjugate(np.transpose(h_abf)), r0), h_abf))
