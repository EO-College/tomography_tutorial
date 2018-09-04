import os
import pickle
import shutil
import numpy as np
from osgeo import gdal
from scipy import ndimage
import subprocess as sp
from .ancillary import cbfi


def start(notebook):
    """
    Create a custom copy of the jupyter notebook with a name defined buy the user and start it.
    The notebook is only copied from the package if it does not yet exist.
    Jupyter notebook files have the extension '.ipynb'.
    If the defined notebook does not contain this extension it is appended automatically.

    Parameters
    ----------
    directory: str
        the name of the custom notebook

    Returns
    -------

    """
    source_basename = 'tutorial.ipynb'
    target_dir = os.path.dirname(notebook)
    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)
    if not notebook.endswith('.ipynb'):
        notebook += '.ipynb'
    if not os.path.isfile(notebook):
        source = os.path.join(os.path.dirname(os.path.realpath(__file__)), source_basename)
        # copy the tutorial notebook from the directory of the installed package to the user directory
        shutil.copyfile(source, notebook)
    sp.check_call(['jupyter', 'notebook', notebook], cwd=target_dir)


def read_data(input, outname, overwrite=False):
    """
    read the raw input data into a numpy array and write the results

    Parameters
    ----------
    input: str or list
        a single image file name or a list of multiple files
    outname: str
        the name of the file to be written.
    overwrite: bool
        overwrite an existing file? Otherwise it is read from file and returned

    Returns
    -------
    numpy.ndarray
        an array in 2D (one file) or 3D (multiple files)
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
    Removal of Topographical Phase.
    If the target file already exists and ``overwrite=False`` this function acts as a simple file reader.

    Parameters
    ----------
    img_stack: numpy.ndarray
        the SLC image stack
    dem_stack: numpy.ndarray
        the image stack containing flat earth and topographic phase
    outname: str
        the name of the file to be written
    overwrite: bool
        overwrite an existing file? Otherwise it is read from file and returned

    Returns
    -------
    numpy.ndarray
        the normalized SLC stack
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
    compute the covariance matrix.
    If the target file already exists and ``overwrite=False`` this function acts as a simple file reader.

    Parameters
    ----------
    img_stack: numpy.ndarray
        the normalized SLC image stack
    outname: str
       the name of the file to be written
    kernelsize: int
        the boxcar smoothing dimension
    overwrite: bool
        overwrite an existing file? Otherwise it is read from file and returned

    Returns
    -------
    numpy.ndarray
        the covariance matrix
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
    perform the capon beam forming inversion to create the final tomographic result.
    If the target file already exists and ``overwrite=False`` this function acts as a simple file reader.

    Parameters
    ----------
    covmatrix: numpy.ndarray
        the covariance matrix
    kz_array: numpy.ndarray
        the wave number stack
    outname: str
        the name of the file to be written
    height: int
        the maximum inversion height
    overwrite: bool
        overwrite an existing file? Otherwise it is read from file and returned

    Returns
    -------
    numpy.ndarray
        the tomographic array
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
