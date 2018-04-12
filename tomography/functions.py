import os
import pickle
import numpy as np
from osgeo import gdal
from scipy import ndimage
import matplotlib.pyplot as plt
from .ancillary import listfiles, plot_profile_horizontal


def read_data(img_list, outname, overwrite=False):
    """
    read the input data

    :param img_list: a list of image files
    :param outname: the name of the file to be written
    :param overwrite: overwrite an existing file? Otherwise it is read from file and returned.
    :return: an array containing the stacked images (numpy.ndarray)
    """

    if overwrite or not os.path.isfile(outname):
        if len(img_list) == 0:
            raise RuntimeError('img_list is empty')
        fil = img_list[0]
        imgfile = gdal.Open(fil)
        rows = imgfile.RasterYSize
        cols = imgfile.RasterXSize
        band = imgfile.GetRasterBand(1)
        imgfile = None
        dtype = gdal.GetDataTypeName(band.DataType)

        if dtype == 'CFloat32':
            offset = 0
        elif dtype == 'Float32':
            offset = 1
        else:
            raise RuntimeError('data type must be either "CFloat32" or "Float32"')

        img_stack = np.empty((rows, cols, len(img_list) + offset), np.complex64)
        for ii, img_path in enumerate(img_list):
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
        rows = img_stack.shape[0]
        cols = img_stack.shape[1]
        nTrack = img_stack.shape[2]

        j_complex = complex(0, 1)

        normalized_stack = np.empty((rows, cols, nTrack), np.complex64)
        normalized_stack[:, :, 0] = img_stack[:, :, 0]

        for ii in range(1, nTrack):
            print('removing the flat-earth and the topographical phase of the slave : ' + str(ii))
            normalized_stack[:, :, ii] = img_stack[:, :, ii] * np.exp(j_complex * dem_stack[:, :, ii])

        # save the variable
        with open(outname, 'wb') as f:
            pickle.dump(normalized_stack, f, 2)
    else:
        with open(outname, 'rb') as f:
            normalized_stack = pickle.load(f)
    return normalized_stack


def calculate_covariance_matrix(img_stack, outname, kernelsize=10, overwrite=False):
    """
    calculate coherence &  covariance matrix

    :param img_stack: the image stack (numpy.ndarray)
    :param outname: the name of the file to be written (str)
    :param kernelsize: the boxcar smoothing dimension (int)
    :param overwrite: overwrite an existing file? Otherwise it is read from file and returned. (bool)
    :return: the covariance matrix (numpy.ndarray)
    """
    if overwrite or not os.path.isfile(outname):
        rows = img_stack.shape[0]
        cols = img_stack.shape[1]
        nTrack = img_stack.shape[2]

        j_complex = complex(0, 1)

        kernel = np.ones((kernelsize, kernelsize))
        cov_matrix = np.empty((rows, cols, nTrack, nTrack), np.complex64)
        weight = kernel / np.sum(kernel)

        for ii in range(0, nTrack):
            for jj in range(0, nTrack):
                img1 = img_stack[:, :, ii]
                img2 = img_stack[:, :, jj]
                ms = np.multiply(img1, np.conjugate(img2))  # multiply S1 and complex conjugate of S2
                mm = np.absolute(img1) ** 2
                ss = np.absolute(img2) ** 2
                ms_smooth = ndimage.convolve(ms.real, weight, mode='reflect') + \
                            j_complex * ndimage.convolve(ms.imag, weight, mode='reflect')
                mm_smooth = ndimage.convolve(mm, weight, mode='reflect')
                ss_smooth = ndimage.convolve(ss, weight, mode='reflect')

                coherence = np.divide(ms_smooth, (np.sqrt(np.multiply(mm_smooth, ss_smooth))))

                cov_matrix[:, :, ii, jj] = coherence

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
        z_vector = np.matrix(np.arange(-height, height + 1, 1))
        h = z_vector.shape[1]
        rows = covmatrix.shape[0]
        cols = covmatrix.shape[1]
        nTrack = covmatrix.shape[2]
        capon_bf = np.empty((rows, cols, h), np.complex64)

        j_complex = complex(0, 1)

        # define the loading factor
        load_fac = (1 / 25.) * np.identity(nTrack)

        for ii_az in range(0, rows):
            for ii_rg in range(0, cols):
                kz = np.transpose([kz_array[ii_az, ii_rg, :]])
                r0 = covmatrix[ii_az, ii_rg, :, :]

                # define the steering matrix
                a_steer = np.exp(np.dot(j_complex * kz, z_vector))

                # define the numerator of the filter habf
                hnum = np.dot(np.linalg.inv(r0 + load_fac), a_steer)

                # define the denominator of the filter habf
                hden0 = np.diag(np.dot(np.conjugate(np.transpose(a_steer)), hnum))

                # replicate the diagonal by number of Tracks
                hden = np.zeros((nTrack, h), dtype=np.complex64)
                for i in range(nTrack):
                    hden[i] = hden0

                h_abf = np.divide(hnum, hden)

                capon_bf[ii_az, ii_rg, :] = np.diag(np.dot(np.dot(np.conjugate(np.transpose(h_abf)), r0), h_abf))

        capon_bf = np.real(capon_bf)
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

    # normalize the values for each pixel
    caponnorm = np.empty((rows, cols, h), np.float_)

    for ii_az in range(0, rows):
        for ii_rg in range(0, cols):
            capon = capon_bf[ii_az, ii_rg, :]
            caponmax = np.amax(capon)
            caponmin = np.amin(capon)
            caponnorm[ii_az, ii_rg, :] = np.divide((capon - caponmin), (caponmax - caponmin))

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
            fname = os.path.join(outpath + 'reflectivity_profile_stack_{:02}.png'.format(len(existing)+1))
            plt.savefig(fname, dpi=1000)
        else:
            plt.show()
