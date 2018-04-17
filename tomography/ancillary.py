
import os
import re
import numpy as np
import matplotlib.pyplot as plt


def plot_profile_horizontal(slice, height, xlab, ylab, title, outname=None):
    """
    plotting of pixel profiles

    :param slice: the profile to be plotted
    :param height: the maximum inversion height (int)
    :param xlab: the x-axis label
    :param ylab: the y-axis label
    :param title: the plot title
    :param outname:  the name of the written file (png format); if not defined the plots are only displayed (str)
    :return: None
    """
    h = slice.shape[1]
    my_yticks = np.arange(-height, height, 10)
    my_tickpos = np.arange(0, h, 10)
    plt.yticks(my_tickpos, my_yticks)
    plt.imshow(np.rot90(slice, 1), origin='lower', cmap='jet')
    plt.xlabel(xlab, fontsize=8)
    plt.ylabel(ylab, fontsize=8)
    plt.title(title, fontsize=8)
    plt.tick_params(axis='both', which='major', labelsize=4)
    if outname:
        plt.savefig(outname, dpi=1000)
        plt.close()
    else:
        plt.show()


def plot_image(slice, xlab='range', ylab='azimuth', title='', outname=None):
    plt.imshow(slice, origin='upper', cmap='jet')
    plt.xlabel(xlab, fontsize=12)
    plt.ylabel(ylab, fontsize=12)
    plt.title(title, fontsize=12)
    plt.colorbar()
    if outname:
        plt.savefig(outname, dpi=1000)
        plt.close()
    else:
        plt.show()


def listfiles(path, pattern):
    """
    list files in a directory whose names match a regular expression
    :param path: the directory to be searched
    :param pattern: the regular expression
    :return: a list of strings with the absolute file names
    """
    return sorted([os.path.join(path, x) for x in os.listdir(path) if re.search(pattern, x)])


def sample(midr, mida, inrange=1):
    """
    generate range and azimuth indices for drawing a set of reflectivity profiles

    if inrange is 1 only the central pixel is returned
    example:
    >>> sample(30, 17, 3)
    [(29, 16), (29, 17), (29, 18), (30, 16), (30, 17), (30, 18), (31, 16), (31, 17), (31, 18)]

    :param midr: the range pixel coordinate
    :param mida: the azimuth pixel coordinate
    :param inrange: the the range of pixels to select around the central pixel
    :return: a list of tuples containing range and azimuth coordinates
    """
    if inrange % 2 == 0:
        raise RuntimeError('parameter inrange must be odd')
    rg_range = range(midr - inrange // 2, midr + inrange // 2 + 1)
    az_range = range(mida - inrange // 2, mida + inrange // 2 + 1)
    list_pixels = [(rg, az) for rg in rg_range for az in az_range]
    return list_pixels


def normalize(slice):
    """

    Parameters
    ----------
    slice: ndarray
        a 1d input array to be normalized
    Returns
    out: ndarray
        the normalized array
    -------

    """
    max = np.amax(slice)
    min = np.amin(slice)
    return np.divide((slice - min), (max - min))


def cbfi(slice, nTrack, height):
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
