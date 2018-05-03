import os
import pickle
import numpy as np
from osgeo import gdal
from scipy import ndimage
import matplotlib.pyplot as plt
from .ancillary import listfiles, plot_profile_horizontal, normalize, cbfi


def read_data(input, outname, overwrite=False):
    """
    read the raw input data to numpy arrays and write the results

    Parameters
    ----------
    input: str or list
        a single image file name or a list of multiple files
    outname: str
        the name of the file to be written. Default is False.
    overwrite: bool
        overwrite an existing file? Otherwise it is read from file and returned. Default is False.
    Returns
    -------
    out: numpy.ndarray
        a 2D (one file) or 3D (multiple files) array
    """

    if overwrite or not os.path.isfile(outname):
        if len(input) == 0:
            raise RuntimeError('img_list is empty')

        imgfile = gdal.Open(input[0])
        rows = imgfile.RasterYSize
        cols = imgfile.RasterXSize
        dtype = gdal.GetDataTypeName(imgfile.GetRasterBand(1).DataType)
        imgfile = None

        if dtype == 'CFloat32':
            offset = 0
        elif dtype == 'Float32':
            offset = 1
        else:
            raise RuntimeError('data type must be either "CFloat32" or "Float32"')

        img_stack = np.empty((rows, cols, len(input) + offset), np.complex64)
        for ii, img_path in enumerate(input):
            # Read the image file into the array
            imgfile = gdal.Open(img_path)
            img_stack[:, :, ii + offset] = imgfile.ReadAsArray()
            imgfile = None

        # save the variable
        with open(outname, 'wb') as f:
            pickle.dump(img_stack, f, 2)
    else:
        with open(outname, 'rb') as f:
            img_stack = pickle.load(f)
    return img_stack


def topo_phase_removal(img_stack, dem_stack, outname, overwrite=False):
    """
    Removal of Topographical Phase

    :param img_stack: the image stack (numpy.ndarray)
    :param dem_stack: (numpy.ndarray)
    :param outname: the name of the file to be written (str)
    :param overwrite: overwrite an existing file? Otherwise it is read from file and returned. (bool)
    :return: a matrix containing the normalized stack (numpy.ndarray)
    """
    if overwrite or not os.path.isfile(outname):

        normalized_stack = img_stack * np.exp(complex(0, 1) * dem_stack, dtype=np.complex64)

        # save the variable
        with open(outname, 'wb') as f:
            pickle.dump(normalized_stack, f, 2)
    else:
        with open(outname, 'rb') as f:
            normalized_stack = pickle.load(f)
    return normalized_stack


def calculate_covariance_matrix(img_stack, outname, kernelsize=10, overwrite=False):
    """
    calculate covariance matrix

    :param img_stack: the image stack (numpy.ndarray)
    :param outname: the name of the file to be written (str)
    :param kernelsize: the boxcar smoothing dimension (int)
    :param overwrite: overwrite an existing file? Otherwise it is read from file and returned. (bool)
    :return: the covariance matrix (numpy.ndarray)
    """

    if overwrite or not os.path.isfile(outname):
        rows, cols, nTrack = img_stack.shape

        kernel = np.ones((kernelsize, kernelsize))
        weight = kernel / np.sum(kernel)

        cov_matrix = np.empty((rows, cols, nTrack, nTrack), np.complex64)
        smooth = np.empty((rows, cols, nTrack), np.float32)

        for i in range(0, nTrack):
            mm = np.absolute(img_stack[:, :, i]) ** 2
            smooth[:, :, i] = ndimage.convolve(mm, weight, mode='reflect')

        for i in range(0, nTrack):
            for j in range(0, nTrack):
                # multiply S1 and complex conjugate of S2
                ms = np.multiply(img_stack[:, :, i], np.conjugate(img_stack[:, :, j]))

                ms_real_smooth = ndimage.convolve(ms.real, weight, mode='reflect')
                ms_imag_smooth = complex(0, 1) * ndimage.convolve(ms.imag, weight, mode='reflect')
                ms_smooth = ms_real_smooth + ms_imag_smooth

                cov_matrix[:, :, i, j] = ms_smooth / (smooth[:, :, i] * smooth[:, :, j]) ** 0.5

        # save the variable
        with open(outname, 'wb') as f:
            pickle.dump(cov_matrix, f, 2)
    else:
        with open(outname, 'rb') as f:
            cov_matrix = pickle.load(f)
    return cov_matrix


def capon_beam_forming_inversion(covmatrix, kz_array, outname, height=70, overwrite=False):
    """

    :param covmatrix:(numpy.ndarray)
    :param kz_array:(numpy.ndarray)
    :param height: the maximum inversion height (int)
    :param outname: the name of the file to be written (str)
    :param overwrite: overwrite an existing file? Otherwise it is read from file and returned. (bool)
    :return: the computed matrix (numpy.ndarray)
    """
    if overwrite or not os.path.isfile(outname):

        rows, cols, nTrack = covmatrix.shape[0:3]

        # reshape and stack the input arrays for easier vectorization
        stack = np.empty((rows, cols, nTrack + nTrack ** 2), np.complex64)
        stack[:, :, 0:(nTrack ** 2)] = covmatrix.reshape((rows, cols, nTrack ** 2))
        stack[:, :, (nTrack ** 2):] = kz_array

        # perform the actual computations
        capon_bf = np.apply_along_axis(cbfi, axis=2, arr=stack, nTrack=nTrack, height=height)

        capon_bf = np.real(capon_bf).astype(np.float32)
        with open(outname, 'wb') as f:
            pickle.dump(capon_bf, f, 2)
    else:
        with open(outname, 'rb') as f:
            capon_bf = pickle.load(f)
    return capon_bf


def plot_tomo_slices(capon_bf, height=70, outpath=None):
    """
    plot slices parallel to azimutha and range as well as horizontal images

    :param capon_bf:
    :param height: the maximum inversion height (int)
    :param outpath: the path to write the plots to (png format); if not defined the plots are only displayed (str)
    :return: None
    """
    rows = capon_bf.shape[0]
    cols = capon_bf.shape[1]
    h = capon_bf.shape[2]

    if outpath:
        subdir = os.path.join(outpath, 'tomo_slices')
        if not os.path.exists(subdir):
            os.makedirs(subdir)

    caponnorm = np.apply_along_axis(normalize, 2, capon_bf)

    ylab = 'height [m]'

    # slices parallel to range
    xlab = 'range pixel'
    for ii_az in range(0, rows, 100):
        title = 'range slice at range line ' + str(ii_az)
        if outpath:
            outname = os.path.join(subdir, 'range_slice_pix_{}.png'.format(ii_az))
        else:
            outname = None
        plot_profile_horizontal(caponnorm[ii_az, :, :], height, xlab, ylab, title, outname)

    # slices parallel to azimuth
    xlab = 'azimuth pixel'
    for ii_rg in range(0, cols, 100):
        title = 'azimuth slice at range line ' + str(ii_rg)
        if outpath:
            outname = os.path.join(subdir, 'azimuth_slice_pix_{}.png'.format(ii_rg))
        else:
            outname = None
        plot_profile_horizontal(caponnorm[:, ii_rg, :], height, xlab, ylab, title, outname)

    # horizontal slices
    for ii_h in range(0, h, 10):
        caponplot_h = caponnorm[:, :, ii_h]
        plt.imshow(caponplot_h, origin='upper', cmap='jet')
        plt.xlabel('range', fontsize=12)
        plt.ylabel('azimuth', fontsize=12)
        plt.title('horizontal slice at height ' + str(height - ii_h) + ' m', fontsize=12)
        if outpath:
            print('writing a png file for the horizontal slice at height {}'.format(height - ii_h))
            fname = os.path.join(subdir, 'horizontal_{}.png'.format(height - ii_h))
            plt.savefig(fname)
            plt.close()
        else:
            plt.show()


def plot_profiles(capon_bf, pixels, height=70, outpath=None, plot_separate=False):
    """
    plot vertical backscatter profiles

    :param capon_bf: the capon beam forming inversion result
    :param pixels: a list of pixel range-azimuth coordinate tuples
    :param height: the maximum inversion height (int)
    :param outpath: the path to write the plots to (png format);
        if not defined the plots are only displayed
    :param plot_separate: plot each ref. profiles separately?
    :return: None
    """

    for rgpix, azpix in pixels:
        x = np.flipud(capon_bf[azpix, rgpix, :])
        y = np.arange(-height, height + 1, 1)
        plt.plot(x, y, label='pixel at rg_{0}_az_{1}'.format(rgpix, azpix))
        plt.ylim(-height, height)
        plt.xlabel('reflectivity')
        plt.ylabel('height [m]')
        plt.title('reflectivity profiles')
        plt.legend(loc=0, prop={'size': 7}, markerscale=1)
        if plot_separate:
            if outpath:
                basename = 'reflectivity_profile_rg_{0}_az_{1}.png'.format(rgpix, azpix)
                fname = os.path.join(outpath, basename)
                plt.savefig(fname, dpi=1000)
            else:
                plt.show()
    if not plot_separate:
        if outpath:
            existing = listfiles(outpath, 'reflectivity_profile_stack_[0-9]{2}.png')
            fname = os.path.join(outpath + 'reflectivity_profile_stack_{:02}.png'.format(len(existing) + 1))
            plt.savefig(fname, dpi=1000)
        else:
            plt.show()
