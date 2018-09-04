import os
from osgeo import gdal, osr
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display
from ipywidgets import IntSlider, IntRangeSlider, Checkbox, Button, Layout, HBox, VBox, interactive_output
from mpl_toolkits.axes_grid1 import make_axes_locatable


class DataViewer(object):
    """
    functionality for displaying the input data (SLC, topographic phase and wave number)

    Parameters
    ----------
    slc_list: list of str
        the names of the SLC input files
    phase_list: list of str
        the names of the topographic phase input files
    kz_list: list of str
        the names of the Kappa-Zeta wave number input files
    slc_stack: numpy.ndarray
        the SLC images
    phase_stack: numpy.ndarray
        the topographic phase images
    kz_stack: numpy.ndarray
        the wave number images
    """

    def __init__(self, slc_list, phase_list, kz_list, slc_stack, phase_stack, kz_stack):
        self.slc_list = slc_list
        self.phase_list = phase_list
        self.kz_list = kz_list

        self.slc_stack = slc_stack
        self.phase_stack = phase_stack
        self.kz_stack = kz_stack

        # define some options for display of the widget box
        self.layout = Layout(
            display='flex',
            flex_flow='row',
            border='solid 2px',
            align_items='stretch',
            width='94.5%'
        )

        # define a slider for changing a plotted image
        self.slider = IntSlider(min=1, max=len(self.slc_list) - 1, step=1, continuous_update=False,
                                description='image number',
                                style={'description_width': 'initial'},
                                layout=self.layout)

        display(self.slider)

        self.fig = plt.figure(num='visualization of input data')
        # display of SLC amplitude
        self.ax1 = self.fig.add_subplot(131)
        # display of topographical phase
        self.ax2 = self.fig.add_subplot(132)
        # display of wave number
        self.ax3 = self.fig.add_subplot(133)

        self.ax1.set_xlabel('range', fontsize=12)
        self.ax1.set_ylabel('azimuth', fontsize=12)

        # format the cursor value displays
        self.ax1.format_coord = lambda x, y: 'range={0}, azimuth={1}, amplitude='.format(int(x), int(y))
        self.ax2.format_coord = lambda x, y: 'range={0}, azimuth={1}, phase='.format(int(x), int(y))
        self.ax3.format_coord = lambda x, y: 'range={0}, azimuth={1}, wave number='.format(int(x), int(y))

        # enable interaction with the slider
        out = interactive_output(self.__onslide, {'h': self.slider})

        plt.tight_layout()

    def __onslide(self, h):
        """
        a function to respond to slider value changes

        Parameters
        ----------
        h: int
            the slider value

        Returns
        -------
        """
        slc_name = os.path.basename(self.slc_list[h])
        phase_name = os.path.basename(self.phase_list[h - 1])
        kz_name = os.path.basename(self.kz_list[h - 1])
        self.ax1.set_title('SLC intensity: {}'.format(slc_name), fontsize=12)
        self.ax2.set_title('phase: {}'.format(phase_name), fontsize=12)
        self.ax3.set_title('wavenumber: {}'.format(kz_name), fontsize=12)
        # logarithmic scaling of SLC amplitude
        amp_log = 10 * np.log10(np.absolute(self.slc_stack[:, :, h]))
        # computation of image percentiles for linear image stretching
        p02, p98 = np.percentile(amp_log, (2, 98))
        self.ax1.imshow(amp_log, cmap='gray', vmin=p02, vmax=p98)
        self.ax2.imshow(np.absolute(self.phase_stack[:, :, h]))
        self.ax3.imshow(np.absolute(self.kz_stack[:, :, h]))
        self._set_colorbar(self.ax1, 'amplitude [$dB$]')
        self._set_colorbar(self.ax2, 'phase [$rad$]')
        self._set_colorbar(self.ax3, 'wave number [$m\cdot rad^{-1}$]')

    def _set_colorbar(self, axis, label):
        if len(axis.images) > 1:
            axis.images[0].colorbar.remove()
            del axis.images[0]

        divider = make_axes_locatable(axis)
        cax = divider.append_axes('right', size='5%', pad=0.05)

        cbar = self.fig.colorbar(axis.images[0], cax=cax)
        cbar.ax.set_ylabel(label, fontsize=12)


class Tomographyplot(object):
    """
    functionality for creating the main tomography_tutorial analysis plot

    Parameters
    ----------
    capon_bf_abs: numpy.ndarray
        the absolute result of the capon beam forming inversion
    caponnorm: numpy.ndarray
        the normalized version of `capon_bf_abs`;
        see function :func:`~tomography_tutorial.ancillary.normalize`.
    """

    def __init__(self, capon_bf_abs, caponnorm):
        if not caponnorm.shape == capon_bf_abs.shape:
            raise RuntimeError('mismatch of input arrays')

        self.height = capon_bf_abs.shape[2] // 2
        self.capon_bf_abs = capon_bf_abs
        self.caponnorm = caponnorm
        #############################################################################################
        # widget box setup

        # define a slider for changing the horizontal slice image
        self.slider = IntSlider(min=-self.height,
                                max=self.height,
                                step=10,
                                continuous_update=False,
                                description='inversion height',
                                style={'description_width': '140px'},
                                layout={'width': '400px'})

        self.rangeslider = IntRangeSlider(value=[-self.height, self.height],
                                          min=-self.height,
                                          max=self.height,
                                          step=1,
                                          continuous_update=False,
                                          description='inversion height range',
                                          style={'description_width': '140px'},
                                          layout={'width': '400px'})

        # an internal variable to keep track of the set inversion height range
        self.heightrange = (-self.height, self.height)

        # a simple checkbox to enable/disable stacking of vertical profiles into one plot
        self.checkbox = Checkbox(value=True, description='stack vertical profiles', indent=False)

        # a button to clear the vertical profile plot
        self.clearbutton = Button(description='clear vertical plot')
        self.clearbutton.on_click(lambda x: self.__init_vertical_plot())

        layout = Layout(
            justify_content='space-around',
            border='solid 2px',
            width='88%'
        )

        form = HBox(children=[VBox([self.slider, self.rangeslider]),
                              VBox([self.checkbox, self.clearbutton])],
                    layout=layout)

        display(form)
        #############################################################################################
        # main plot setup

        # set up the subplot layout
        self.fig = plt.figure()
        # the horizontal slice plot
        self.ax1 = self.fig.add_subplot(221)
        # the vertical profile plot
        self.ax2 = self.fig.add_subplot(222)
        # the range slice plot
        self.ax3 = self.fig.add_subplot(413)
        # the azimuth slice plot
        self.ax4 = self.fig.add_subplot(414)
        plt.subplots_adjust(left=0.1, right=0.2, top=0.3, bottom=0.2)

        # set up the plots for range and azimuth slices
        self.ax3.set_xlim(0, capon_bf_abs.shape[1])
        self.ax4.set_xlim(0, capon_bf_abs.shape[0])
        self.ax3.set_ylim(0, self.height * 2)
        self.ax4.set_ylim(0, self.height * 2)

        # set up the vertical profile plot
        self.__init_vertical_plot()

        # add a cross-hair to the horizontal slice plot
        self.lhor = self.ax1.axhline(0, linewidth=1, color='r')
        self.lver = self.ax1.axvline(0, linewidth=1, color='r')

        # make the figure respond to mouse clicks by executing method __onclick
        self.cid1 = self.fig.canvas.mpl_connect('button_press_event', self.__onclick)

        # enable interaction with the sliders
        out1 = interactive_output(self.__onslide, {'h': self.slider})
        out2 = interactive_output(self.__onslide_range, {'h': self.rangeslider})
        #############################################################################################
        # general formatting

        # format the cursor value displays
        self.ax1.format_coord = lambda x, y: \
            'range={0}, azimuth={1}, reflectivity='.format(int(x), int(y))
        self.ax2.format_coord = lambda x, y: \
            'reflectivity={0:.3f}, height={1}'.format(x, int(y))
        self.ax3.format_coord = lambda x, y: \
            'range={0}, height={1}, reflectivity='.format(int(x), int(y - self.height))
        self.ax4.format_coord = lambda x, y: \
            'range={0}, height={1}, reflectivity='.format(int(x), int(y - self.height))

        # arrange the subplots to make best use of space
        plt.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.1)
        #############################################################################################

    def __reset_crosshair(self, range, azimuth):
        """
        redraw the cross-hair on the horizontal slice plot

        Parameters
        ----------
        range: int
            the range image coordinate
        azimuth: int
            the azimuth image coordinate

        Returns
        -------
        """
        self.lhor.set_ydata(azimuth)
        self.lver.set_xdata(range)
        plt.draw()

    def __init_vertical_plot(self):
        """
        set up the vertical profile plot

        Returns
        -------
        """
        # clear the plot if lines have already been drawn on it
        if len(self.ax2.lines) > 0:
            self.ax2.cla()
        # set up the vertical profile plot
        self.ax2.set_ylabel('height [m]', fontsize=12)
        self.ax2.set_xlabel('reflectivity', fontsize=12)
        self.ax2.set_title('vertical point profiles', fontsize=12)
        self.ax2.set_ylim(self.heightrange[0], self.heightrange[1])

    def __rename_sliceplot_ticklabels(self):
        """
        rename the ticks of the slice plots from pixel coordinates (0:nbands) to
        inversion height range (-height:height)

        Returns
        -------

        """
        ticks = self.ax3.get_yticks()
        labels_new = [str(int(x - self.height)) for x in ticks]
        self.ax3.set_yticklabels(labels_new)
        self.ax4.set_yticklabels(labels_new)

    def __onslide(self, h):
        """
        a function to respond to slider value changes by redrawing the horizontal slice plot

        Parameters
        ----------
        h: int
            the slider value
        Returns
        -------
        """
        p1 = self.ax1.imshow(self.caponnorm[:, :, self.height - h], origin='upper', cmap='jet')
        self.ax1.set_xlabel('range', fontsize=12)
        self.ax1.set_ylabel('azimuth', fontsize=12)
        self.ax1.set_title('horizontal slice at height {} m'.format(h), fontsize=12)
        # remove the previous plot and its color bar
        if len(self.ax1.images) > 1:
            self.ax1.images[0].colorbar.remove()
            del self.ax1.images[0]
        cbar = self.fig.colorbar(p1, ax=self.ax1)
        cbar.ax.set_ylabel('reflectivity', fontsize=12)
        plt.show()

    def __onslide_range(self, h):
        """
        respond to changes on the range slider.
        This changes the displayed y-axis range of the vertical profile and
        slice plots and renames the tick labels
        Parameters
        ----------
        h

        Returns
        -------

        """
        self.heightrange = h
        self.ax2.set_ylim(self.heightrange[0], self.heightrange[1])
        self.ax3.set_ylim(self.heightrange[0] + self.height, self.heightrange[1] + self.height)
        self.ax4.set_ylim(self.heightrange[0] + self.height, self.heightrange[1] + self.height)
        self.__rename_sliceplot_ticklabels()

    def __onclick(self, event):
        """
        respond to mouse clicks in the plot.
        This function responds to clicks on the first (horizontal slice) plot and
        updates the vertical profile and slice plots

        Parameters
        ----------
        event: matplotlib.backend_bases.MouseEvent
            the click event object containing image coordinates

        """
        # only do something if the first plot has been clicked on
        if event.inaxes == self.ax1:

            # retrieve the click coordinates
            rg = int(event.xdata)
            az = int(event.ydata)

            # redraw the cross-hair
            self.__reset_crosshair(rg, az)

            subset_vertical = self.capon_bf_abs[az, rg, :]
            subset_range = self.caponnorm[az, :, :]
            subset_azimuth = self.caponnorm[:, rg, :]

            # redraw/clear the vertical profile plot in case stacking is disabled
            if not self.checkbox.value:
                self.__init_vertical_plot()

            # plot the vertical profile
            label = 'rg: {0:03}; az: {1:03}'.format(rg, az)
            self.ax2.plot(np.flipud(subset_vertical), range(-self.height, self.height + 1), label=label)
            self.ax2_legend = self.ax2.legend(loc=0, prop={'size': 7}, markerscale=1)

            # plot the range slice
            self.ax3.imshow(np.rot90(subset_range, 1), origin='lower', cmap='jet', aspect='auto')
            self.ax3.set_title('range slice at azimuth line {}'.format(az), fontsize=12)

            # plot the azimuth slice
            self.ax4.imshow(np.rot90(subset_azimuth, 1), origin='lower', cmap='jet', aspect='auto')
            self.ax4.set_title('azimuth slice at range line {}'.format(rg), fontsize=12)

            self.__rename_sliceplot_ticklabels()


class GeoViewer(object):
    """
    plotting utility for displaying a geocoded image stack file.

    On moving the slider, the band at the slider position is read from the file and displayed.

    Parameters
    ----------
    filename: str
        the name of the file to display
    cmap: str
        the color map for displaying the image. See :func:`matplotlib.pyplot.imshow`.
    band_indices: list
        a list of indices for renaming the individual bands in `filename` such that one can
        scroll trough the range of inversion heights, e.g. -70:70, instead of the raw band
        indices, e.g. 1:140.
        The number of unique elements must of same length as the number of bands in `filename`.
    """

    def __init__(self, filename, cmap='jet', band_indices=None):
        self.filename = filename
        ras = gdal.Open(filename)
        self.rows = ras.RasterYSize
        self.cols = ras.RasterXSize
        self.bands = ras.RasterCount
        geo = ras.GetGeoTransform()
        srs = osr.SpatialReference(wkt=ras.GetProjection())
        ras = None
        srs.AutoIdentifyEPSG()
        epsg = int(srs.GetAuthorityCode(None))

        xmin = geo[0]
        ymax = geo[3]
        xres = geo[1]
        yres = abs(geo[5])

        xmax = xmin + xres * self.cols
        ymin = ymax - yres * self.rows

        self.extent = (xmin, xmax, ymin, ymax)

        # define some options for display of the widget box
        self.layout = Layout(
            display='flex',
            flex_flow='row',
            border='solid 2px',
            align_items='stretch',
            width='88%'
        )

        self.colormap = cmap

        if band_indices is not None:
            if len(list(set(band_indices))) != self.bands:
                raise RuntimeError('length mismatch of unique provided band indices ({0}) '
                                   'and image bands ({1})'.format(len(band_indices), self.bands))
            else:
                self.indices = sorted(band_indices)
        else:
            self.indices = range(1, self.bands + 1)

        # define a slider for changing a plotted image
        self.slider = IntSlider(min=min(self.indices),
                                max=max(self.indices),
                                step=1,
                                continuous_update=False,
                                value=self.indices[len(self.indices)//2],
                                description='band index',
                                style={'description_width': 'initial'},
                                layout=self.layout)

        display(self.slider)

        self.fig = plt.figure()
        self.ax = plt.gca()
        self.ax.get_xaxis().get_major_formatter().set_useOffset(False)
        self.ax.get_yaxis().get_major_formatter().set_useOffset(False)

        self.ax.format_coord = lambda x, y: \
            'easting={0:.2f}, northing={1:.2f}, reflectivity='.format(x, y)

        # enable interaction with the slider
        out = interactive_output(self.__onslide, {'h': self.slider})

    def __onslide(self, h):
        mat = self.__read_band(self.indices.index(h) + 1)
        self.ax.imshow(mat, extent=self.extent, cmap=self.colormap)

    def __read_band(self, band):
        ras = gdal.Open(self.filename)
        mat = ras.GetRasterBand(band).ReadAsArray()
        ras = None
        return mat
