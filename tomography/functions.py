
import os
import pickle
import numpy as np
from osgeo import gdal
from scipy import ndimage
import matplotlib.pyplot as plt

# *****************************************************************************
# =============================================================================
#     ==================== Part 1: read the input data ====================
# =============================================================================

nTrack = 0
height = 0
outpath = '?'
j_complex = '?'


def read_data(img_list, fout):

    print('reading the data')
    fil = img_list[0]
    imgfile = gdal.Open(fil)
    rows = imgfile.RasterYSize
    cols = imgfile.RasterXSize
    band = imgfile.GetRasterBand(1)
    dtype = gdal.GetDataTypeName(band.DataType)
    print(rows, cols, dtype)

    if dtype == 'CFloat32':
        img_stack = np.empty((rows, cols, nTrack), np.complex64)
        for ii, img_path in enumerate(img_list):
            imgfile = gdal.Open(img_path)
            # Read the image file
            img = imgfile.ReadAsArray()
            # store each images in the list
            img_stack[:, :, ii] = img

    elif dtype == 'Float32':
        img_stack = np.empty((rows, cols, nTrack), np.float_)
        img_stack[:, :, 0] = 0
        for ii, img_path in enumerate(img_list):
            imgfile = gdal.Open(img_path)
            # Read the image file
            img = imgfile.ReadAsArray()
            # store each images in the list
            img_stack[:, :, ii+1] = img

    # save the variable
    f = open(os.path.join(outpath + fout), 'wb')
    pickle.dump(img_stack, f, 2)

    return img_stack

# =============================================================================
#     ============== Part 2: Removal of Topographical Phase ==============
# =============================================================================


def topo_phase_removal(img_stack, dem_stack):

    rows = img_stack.shape[0]
    cols = img_stack.shape[1]

    normalized_stack = np.empty((rows, cols, nTrack), np.complex64)
    normalized_stack[:, :, 0] = img_stack[:, :, 0]

    for ii in range(1, nTrack):
        print('removing the flat-earth and the topographical phase of the slave : ' + str(ii))
        normalized_stack[:, :, ii] = img_stack[:, :, ii] * np.exp(j_complex * dem_stack[:, :, ii])

    # save the variable
    f = open(os.path.join(outpath + 'normalized_stack'), 'wb')
    pickle.dump(normalized_stack, f, 2)

    return normalized_stack

# =============================================================================
#     ========= Part 3: Calculate Coherence &  Covariance Matrix =========
# =============================================================================


def calculate_covariance_matrix(img_stack, kernelsize):

    rows = img_stack.shape[0]
    cols = img_stack.shape[1]

    kernel = np.ones((kernelsize, kernelsize))
    cov_matrix = np.empty((rows, cols, nTrack,  nTrack), np.complex64)
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
    f = open(os.path.join(outpath, 'cov_matrix'), 'wb')
    pickle.dump(cov_matrix, f, 2)

    return cov_matrix

# =============================================================================
#     ============== Part 4: CAPON BEAM FORMING INVERSION ==============
# =============================================================================


def capon_beam_forming_inversion(covmatrix, kz_array, z_vector):

    h = z_vector.shape[1]
    rows = covmatrix.shape[0]
    cols = covmatrix.shape[1]
    print(rows, cols, h)
    capon_bf = np.empty((rows, cols, h), np.complex64)

    # define the loading factor
    load_fac = (1/25.)*np.identity(nTrack)

    for ii_az in range(0, rows):
        for ii_rg in range(0, cols):
            kz = np.transpose([kz_array[ii_az, ii_rg, :]])
            r0 = covmatrix[ii_az, ii_rg, :, :]

            # define the steering matrix
            a_steer = np.exp(np.dot(j_complex*kz, z_vector))

            # define the numerator of the filter habf
            hnum = np.dot(np.linalg.inv(r0+load_fac), a_steer)

            # define the denominator of the filter habf
            hden0 = np.diag(np.dot(np.conjugate(np.transpose(a_steer)), hnum))

            # replicate the diagonal by number of Tracks
            hden = np.zeros((nTrack, h), dtype=np.complex64)
            for i in range(nTrack):
                hden[i] = hden0

            h_abf = np.divide(hnum, hden)

            capon_bf[ii_az, ii_rg, :] = np.diag(np.dot(np.dot(np.conjugate(np.transpose(h_abf)), r0), h_abf))

    # save the variable capon_bf
    f = open(os.path.join(outpath + 'capon_bf'), 'wb')
    pickle.dump(np.real(capon_bf), f, 2)

    return np.real(capon_bf)

# =============================================================================
#     ============== Part 5: Plot the Tomography Slices ==============
# =============================================================================


def plot_tomo_slices(capon_bf):

    rows = capon_bf.shape[0]
    cols = capon_bf.shape[1]
    h = capon_bf.shape[2]
    print(rows, cols, h)

    if not os.path.exists(outpath + 'tomo_slices'):
        os.makedirs(outpath + 'tomo_slices')

    # normalize the values for each pixel
    caponnorm = np.empty((rows, cols, h), np.float_)

    for ii_az in range(0, rows):
        for ii_rg in range(0, cols):
            capon = capon_bf[ii_az, ii_rg, :]
            caponmax = np.amax(capon)
            caponmin = np.amin(capon)
            caponnorm[ii_az, ii_rg, :] = np.divide((capon-caponmin), (caponmax-caponmin))

    # slices paralel to range
    for ii_az in range(0, rows, 100):
        print('writing a png file for the range slice at azimuth ' + str(ii_az))
        fname = os.path.join(outpath + 'tomo_slices\\range_slice_pix_' + str(ii_az) + '.png')
        caponplot_rg = caponnorm[ii_az, :, :]
        my_yticks = np.arange(-height, height, 10)
        my_tickpos = np.arange(0, h, 10)
        plt.yticks(my_tickpos, my_yticks)
        plt.imshow(np.rot90(caponplot_rg, 1), origin='lower', cmap='jet')
        plt.xlabel('range pixel', fontsize=8)
        plt.ylabel('height [m]', fontsize=8)
        plt.title('range slice at azimuth line ' + str(ii_az), fontsize=8)
        plt.tick_params(axis='both', which='major', labelsize=4)
        # plt.show()
        plt.savefig(fname, dpi=1000)
        plt.close()

    # slices paralel to azimuth
    for ii_rg in range(0, cols, 100):
        print('writing a png file for the azimuth slice at range ' + str(ii_rg))
        fname = os.path.join(outpath + 'tomo_slices\\azimuth_slice_pix' + str(ii_rg) + '.png')
        caponplot_az = caponnorm[:, ii_rg, :]
        my_yticks = np.arange(-height, height, 10)
        my_tickpos = np.arange(0, h, 10)
        plt.yticks(my_tickpos, my_yticks)
        plt.imshow(np.rot90(caponplot_az, 1), origin='lower', cmap='jet')
        plt.xlabel('azimuth pixel', fontsize=8)
        plt.ylabel('height [m]', fontsize=8)
        plt.title('azimuth slice at range line ' + str(ii_rg), fontsize=8)
        plt.tick_params(axis='both', which='major', labelsize=4)
        # plt.show()
        plt.savefig(fname, dpi=1000)
        plt.close()

    # horizontal slices
    for ii_h in range(0, h, 10):
        print('writing a png file for the horizontal slice at height ' + str(height-ii_h))
        fname = os.path.join(outpath + 'tomo_slices\\horizontal_' + str(height-ii_h) + '.png')
        caponplot_h = caponnorm[:, :, ii_h]
        plt.imshow(caponplot_h, origin='upper', cmap='jet')
        plt.xlabel('range', fontsize=12)
        plt.ylabel('azimuth', fontsize=12)
        plt.title('horizontal slice at height ' + str(height-ii_h) + ' m', fontsize=12)
        plt.savefig(fname)
        plt.close()

# =============================================================================
#     ========= Part 6: Plot the Backscatterer Profile =========
# =============================================================================


def plot_profiles(capon_bf, pixels):

    npix = len(pixels)
    plt.figure()

    for ii_pix in range(0, npix):
        print('writing a png file for the profile at pixel ' + str(ii_pix))
        # to plot each ref. profiles separately, remove '#' from the below code
        # plt.figure()
        pix = pixels[ii_pix]
        rgpix = pix[0]
        azpix = pix[1]
        fname = os.path.join(outpath + 'reflectivity_profile_rg_' + str(rgpix) + '_az_' + str(azpix) + '.png')
        x = np.flipud(capon_bf[azpix, rgpix, :])
        y = np.arange(-height, height+1, 1)
        plt.plot(x, y, label='pixel at rg_' + str(rgpix) + '_az_' + str(azpix))
        plt.ylim(-height, height)
        plt.xlabel('reflectivity')
        plt.ylabel('height [m]')
        plt.title('reflectivity profiles')
        plt.legend(loc=0, prop={'size': 7}, markerscale=1)
        # plt.show()
        plt.savefig(fname, dpi=1000)