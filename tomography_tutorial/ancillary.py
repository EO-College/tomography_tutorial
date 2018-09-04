import os
import re
import numpy as np

from osgeo import gdal, gdal_array

from scipy.ndimage.measurements import label
from scipy.ndimage import generate_binary_structure, find_objects


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
    >>> listfiles('/path/to/somedata', 'file[0-9].tif')
    ['/path/to/somedata/file1.tif', '/path/to/somedata/file2.tif', '/path/to/somedata/file3.tif']
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
        the 1d input array to be normalized

    Returns
    -------
    numpy.ndarray
        the normalized array
    """
    max = np.amax(slice)
    min = np.amin(slice)
    return np.divide((slice - min), (max - min))


def cbfi(slice, nTrack, height):
    """
    computation of capon beam forming inversion for a single pixel.
    This function is used internally by the core function
    :func:`~tomography_tutorial.functions.capon_beam_forming_inversion`.

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


def lut_crop(lut_rg, lut_az,
             range_min=0, range_max=None,
             azimuth_min=0, azimuth_max=None):
    """
    compute indices for subsetting the range and azimuth lookup tables (LUTs).
    The returned slices describe the minimum LUT subset,
    which contains all radar coordinates within the range-azimuth subset.

    Parameters
    ----------
    lut_rg: numpy.ndarray
        the lookup table for range direction
    lut_az: numpy.ndarray
        the lookup table for azimuth direction
    range_min: int
        first range pixel
    range_max: int
        last range pixel
    azimuth_min: int
        first azimuth pixel
    azimuth_max: int
        last azimuth pixel

    Returns
    -------
    tuple of slices
        the pixel indices for subsetting: (ymin:ymax, xmin:xmax)
    """

    def get_indices(mask):
        """
        helper function to get the array slices containing all rows and columns
        where a binary mask is one/True
        Parameters
        ----------
        mask: numpy.ndarray
            the mask to be subsetted

        Returns
        -------
        tuple of slices
            two slices for subsetting numpy arrays in the row (azimuth) and column (range) dimension
        """
        # get the azimuth array indices for subsetting the LUT files
        az_index = np.where(np.any(mask, axis=1))
        azmin = np.min(az_index)
        azmax = np.max(az_index) + 1

        # get the range array indices for subsetting the LUT files
        rg_index = np.where(np.any(mask, axis=0))
        rgmin = np.min(rg_index)
        rgmax = np.max(rg_index) + 1
        return slice(azmin, azmax), slice(rgmin, rgmax)

    # create a LUT mask containing all radar coordinates of the defined subset
    lut_mask = (float(range_min) <= lut_rg) & \
               (lut_rg <= float(range_max)) & \
               (float(azimuth_min) <= lut_az) & \
               (lut_az <= float(azimuth_max))

    # subset the original mask to the computed subset coordinates
    indices = get_indices(lut_mask)
    lut_mask_sub = lut_mask[indices]

    # refine the mask using an image segmentation approach
    # create a new mask by selecting only the largest segment of the original mask
    s = generate_binary_structure(2, 2)
    labeled_array, num_features = label(lut_mask_sub, structure=s)
    obj_slices = find_objects(labeled_array)
    # identify the largest object
    obj_sizes = [(x.stop - x.start) * (y.stop - y.start) for x, y in obj_slices]
    obj_max = obj_sizes.index(max(obj_sizes))

    # create a new mask subset containing only the largest object
    lut_mask_sub = labeled_array == (obj_max + 1)

    # create a new mask with original extent and place the values of the mask subset
    lut_mask = np.zeros(lut_mask.shape)
    lut_mask[indices] = lut_mask_sub

    # return the position indices of the mask object on the main mask
    return get_indices(lut_mask)


def geowrite(data, outname, reference, indices, nodata=-99):
    """
    write an array to a file using an already geocoded file as reference.
    The output format is either GeoTiff (for 2D arrays) or ENVI (for 3D arrays).

    Parameters
    ----------
    data: numpy.ndarray
        the array to write to the file; must be either 2D or 3D
    outname: str
        the file name
    reference: `gdal.Dataset <https://gdal.org/python/osgeo.gdal.Dataset-class.html>`_
        the geocoded reference dataset
    indices: tuple of slices
        the slices which define the subset of data in reference,
        i.e. where is the data located within the reference pixel dimensions;
        see :func:`lut_crop`
    nodata: int
        the nodata value to write to the file

    Returns
    -------

    """

    geo_keys = ['xmin', 'xres', 'rotation_x', 'ymax', 'rotation_y', 'yres']
    geo = dict(zip(geo_keys, reference.GetGeoTransform()))

    row_f, row_l = indices[0].start, indices[0].stop
    col_f, col_l = indices[1].start, indices[1].stop

    if data.ndim == 2:
        nrow, ncol = data.shape
        nbands = 1
    elif data.ndim == 3:
        nrow, ncol, nbands = data.shape
    else:
        raise RuntimeError("parameter 'data' must be an array with either two or three dimensions")

    if (row_l - row_f, col_l - col_f) != (nrow, ncol):
        raise IndexError('mismatch of data dimensions and subset indices')

    geo['xmin'] += col_f * geo['xres']
    geo['ymax'] -= row_f * abs(geo['yres'])
    geotransform = [geo[x] for x in geo_keys]

    driver = gdal.GetDriverByName('GTiff' if nbands == 1 else 'ENVI')
    dtype = gdal_array.NumericTypeCodeToGDALTypeCode(data.dtype)
    outDataset = driver.Create(outname, ncol, nrow, nbands, dtype)
    driver = None
    outDataset.SetMetadata(reference.GetMetadata())
    outDataset.SetGeoTransform(geotransform)
    outDataset.SetProjection(reference.GetProjection())

    if nbands == 1:
        outband = outDataset.GetRasterBand(1)
        outband.SetNoDataValue(nodata)
        outband.WriteArray(data)
        outband.FlushCache()
    else:
        for i in range(nbands):
            outband = outDataset.GetRasterBand(i+1)
            outband.SetNoDataValue(nodata)
            outband.WriteArray(data[:, :, i])
            outband.FlushCache()

    ref_data = None
    outband = None
    outDataset = None


def geocode(data, lut_rg_name, lut_az_name, outname=None,
            range_min=0, range_max=None, azimuth_min=0, azimuth_max=None):
    """
    Geocode a radar image using lookup tables.
    The LUTs are expected to be georeferenced and contain range and azimuth radar coordinates
    for a specific image data set which is linked to these LUTs.
    If parameter `data` is a subset of this data set,
    the pixel coordinates of this subset need to be defined.

    Parameters
    ----------
    data: numpy.ndarray
        the image data in radar coordinates
    lut_rg_name: str
        the name of the range coordinates lookup table file
    lut_az_name: str
        the name of the azimuth coordinates lookup table file
    outname: str or None
        the name of the file to write;
        if None, the geocoded array is returned and no file is written.
        See function :func:`geowrite` for details on how the file is written.
    range_min: int
        the minimum range coordinate
    range_max: int
        the maximum range coordinate
    azimuth_min: int
        the minimum azimuth coordinate
    azimuth_max: int
        the maximum azimuth coordinate

    Returns
    -------

    Example
    -------
    >>> from osgeo import gdal
    >>> from tomography_tutorial.ancillary import geocode
    >>> image_name = 'path/to/somedata/image.tif'
    >>> lut_rg_name = 'path/to/somedata/lut_rg.tif'
    >>> lut_az_name = 'path/to/somedata/lut_az.tif'
    >>> outname = 'path/to/somedata/image_sub_geo.tif'
    >>> image_ras = gdal.Open(image_name)
    >>> image_mat = image_ras.ReadAsArray()
    >>> image_ras = None
    >>> image_mat_sub = image_mat[0:100, 200:400]
    >>> geocode(image_mat_sub, outname, lut_rg_name, lut_az_name, \
range_min=200, range_max=400, azimuth_min=0, azimuth_max=100)
    """
    if data.ndim == 2:
        nazimuth, nrange = data.shape
        nbands = 1
    elif data.ndim == 3:
        nazimuth, nrange, nbands = data.shape
    else:
        raise RuntimeError("parameter 'data' must be an array with either two or three dimensions")

    if not range_max:
        range_max = nrange

    if not azimuth_max:
        azimuth_max = nazimuth

    imgfile = gdal.Open(lut_rg_name)
    lut_rg = imgfile.ReadAsArray()

    imgfile = gdal.Open(lut_az_name)
    lut_az = imgfile.ReadAsArray()

    # compute indices for subsetting the LUTs to only contain the subset of selected radar coordinates
    indices = lut_crop(lut_rg, lut_az, range_min, range_max, azimuth_min, azimuth_max)

    # crop the LUTs to the selected subset
    lut_rg_sub = lut_rg[indices]
    lut_az_sub = lut_az[indices]

    # recompute the LUT radar coordinates
    lut_rg_sub = lut_rg_sub - range_min
    lut_az_sub = lut_az_sub - azimuth_min

    ########################################################################################################
    # actual geocoding

    # create a mask containing only radar coordinates of the selected subset
    mask = ((lut_rg_sub < 0) | (lut_rg_sub > nrange)) | ((lut_az_sub < 0) | (lut_az_sub > nazimuth))

    # combine the two LUTs by computing indices for a 1-D row-major flattened array
    lut_combi = lut_az_sub.astype(int) * nrange + lut_rg_sub.astype(int)

    # set all indices out of bounds to zero
    lut_combi[(lut_combi < 0) | ((nrange * nazimuth) < lut_combi)] = 0

    if nbands == 1:
        # create the geo-coded array by flattening the radar image to 1D and indexing it by the combined LUT
        # the resulting array has the same dimensions as the LUT, which is still 2D
        mat_geo = data.flatten()[lut_combi]

        # mask out all areas outside the selected subset
        mat_geo[mask] = np.nan
    else:
        mat_geo = np.empty(lut_combi.shape + (nbands,), data.dtype)
        for i in range(nbands):
            sub = data[:, :, i].flatten()[lut_combi]
            sub[mask] = np.nan
            mat_geo[:, :, i] = sub
    ########################################################################################################
    # write the result to disk
    if outname:
        geowrite(mat_geo, outname, imgfile, indices)
    imgfile = None
    if not outname:
        return mat_geo
