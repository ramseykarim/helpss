import numpy as np
import matplotlib.pyplot as plt
from astropy.io.fits import getdata as fits_getdata
from astropy.wcs import WCS
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit

import utils_regrid as rgu
import path_config as cfg
import utils_planck as plu


class GNILCModel:
    band_dictionary_template = {'F545': None, 'F857': None}

    def __init__(self, *target_args, target_bandpass='PACS160um', extension=0):
        """
        Object representing the GNILC dust model and packaging its capability
        to predict flux in various instrument bands.
        This object also contains all the statistical tools needed to calculate
        the target's most common offset from the GNILC-predicted flux, under a
        mask where the GNILC model can be trusted. This  mask is calculated here
        as well by predicting Planck HFI fluxes and comparing back to the observed
        HFI flux used to generate the GNILC model.

        EXAMPLE USE (basic):
        >>> model = GNILCModel(pacs_filename)
        >>> model.get_offset()
        this will print the offset and show the histogram diagnostic plot

        :param target_args: either a filename (str) or a tuple(array, header)
            describing the target whose offset is sought.
            If filename, note the FITS extension argument (default 0)
        :param target_bandpass: stub describing target filter profile.
            Must be one of the (Herschel) filters supported by this code
        :param extension: FITS extension of the data; relevant only if
            target_args gives a string file path
        """
        if len(target_args) == 1 and type(target_args[0]) == str:
            # String file path
            target_filename = target_args[0]
            target_data, target_head = fits_getdata(target_filename, ext=extension, header=True)
        elif len(target_args) == 2:
            # Data, header array
            target_data, target_head = target_args
        else:
            msg = f"Unsupported argument: {target_args}."
            msg += " Accepted formats are: str(filename) or tuple(data, header)."
            raise TypeError(msg)
        self.target_data = target_data
        self.target_wcs = WCS(target_head)
        self.target_bandpass_stub = target_bandpass
        self.projector = rgu.HEALPix2FITS(target_data, target_head, pixel_scale_arcsec=300)
        self.T = self.load_component('Temperature')
        self.beta = self.load_component('Spectral-Index')
        self.tau = self.load_component('Opacity')
        # Declare variables to be set within functions
        self.masks = GNILCModel.band_dictionary_template.copy()
        self.predicted_target_flux = None
        self.mask = None
        self.difference = None
        self.stats = {}
        # Basic operations with reusable results
        self.accumulate_masks()
        self.difference_to_target()
        self.offset_statistics()

    @plu.shape_of_component
    def load_component(self, stub):
        """
        Load a GNILC model component and project it to intermediate grid.
        Decorator adjusts the shape of the array.
        :param stub: GNILC component stub
        :return: GNILC component map at intermediate pixel scale (pixel_scale_arcsec)
        """
        component_filename = cfg.PlanckConfig.component_filename(stub)
        self.projector.healpix_to_intermediate(rgu.open_healpix(component_filename))
        return self.projector.pop_intermediate()

    @plu.shape_of_frequency
    def load_bandpass(self, stub):
        """
        Load the frequency and transmission arrays for a filter profile.
        Decorator adjusts the shape of the array.
        :param stub: Planck HFI or Herschel band stub
        :return: tuple(frequency (Hz) array, transmission array)
        """
        return cfg.get_bandpass_data(stub)

    def predict_flux_in_band(self, band_stub):
        """
        Predict the flux in the specified band.
        Integrate flux at intermediate pixel scale for efficiency.
        Return map gridded to target resolution.
        :param band_stub: Planck HFI or Herschel band stub
        :return: flux map (2D array) in MJy/sr matching target grid
        """
        frequency, transmission = self.load_bandpass(band_stub)
        band_center = cfg.get_bandpass_center(band_stub)
        image = plu.calculate_gnilc_flux(band_center, frequency, transmission,
                                         self.T, self.tau, self.beta)
        return self.projector.intermediate_to_target(intermediate=image)

    def generate_planck_ratio_mask(self, planck_band_stub):
        """
        For a single Planck HFI band, load and remap the observed flux,
        and use the GNILC model to predict the flux in that band.
        Take the ratio of observed to predicted.
        Generate a mask based on constraining the observed-to-predicted
        ratio to within 10% (0.9 < ratio < 1.1).
        Save this mask.
        The logic here is that wherever the Planck GNILC model agrees to
        within 10% with the measurements used to  make this model, it must
        be an appropriate model for predicting flux. Thus, in these regions,
        it should also predict the target (probably PACS 160 micron) flux
        fairly well.
        :param planck_band_stub: HFI band stub
        """
        planck_map_fn = cfg.PlanckConfig.light_map_filename(planck_band_stub)
        observed_planck_flux = self.projector.project(rgu.open_healpix(planck_map_fn))
        observed_planck_flux = cfg.PlanckConfig.unit_conversion(planck_band_stub,
                                                                observed_planck_flux)
        predicted_planck_flux = self.predict_flux_in_band(planck_band_stub)
        ratio = (observed_planck_flux / predicted_planck_flux)
        # Constrain ratio to be within 10% (and not Nan)
        self.masks[planck_band_stub] = ~np.isnan(ratio) & (ratio > 0.9) & (ratio < 1.1)

    def accumulate_masks(self):
        """
        Generate the Planck observed-to-predicted ratio masks for all
        bands listed in the band_dictionary_template.
        Take a bitwise AND of the masks to create the final mask.
        This final mask is used with the difference image to generate
        the offset statistics.
        """
        for band_stub in self.masks:
            self.generate_planck_ratio_mask(band_stub)
        self.mask = np.all(np.array(list(self.masks.values())), axis=0)

    def difference_to_target(self):
        """
        Predict the flux in the target band using the GNILC model.
        Save that map. Convolve the observed target data up to the GNILC
        resolution and create & save a difference image.
        """
        # Assume this map to be at GNILC resolution
        self.predicted_target_flux = self.predict_flux_in_band(self.target_bandpass_stub)
        # Convolve Herschel map to Planck GNILC resolution
        h_beam_size = cfg.HerschelConfig.beam_size(self.target_bandpass_stub)
        p_beam_size = cfg.gnilc_resolution
        convolution_beam_size = np.sqrt(p_beam_size ** 2 - h_beam_size ** 2)
        convolution_beam = rgu.prepare_convolution(self.target_wcs,
                                                   convolution_beam_size,
                                                   self.target_data.shape)
        target_data_convolved = rgu.convolve_properly(self.target_data, convolution_beam)
        # [predicted] - [observed] difference map
        self.difference = self.predicted_target_flux - target_data_convolved

    def offset_statistics(self, n_bins=128, x_range=(25, 80)):
        """
        Make masked difference image histogram and use it to calculate
        a few basic statistics. Save all this to the self.stats dictionary.
        :param n_bins: number of histogram bins. Default=128
        :param x_range: histogram limits. Default=(25, 80).
            This range is acceptable for the Perseus test region, but may
            need to be tweaked for future regions.
        """
        # Run numpy histogram function with default bins and limits
        # Rerun this with new bins/limits if these don't work
        d_hist, d_edges = np.histogram(self.difference[self.mask].ravel(),
                                       bins=n_bins, range=x_range)
        # Set up plot-friendly arrays
        hist_x, hist_y = map(lambda x: np.array([x[0], x[1]]).T.flatten(),
                             ((d_edges[:-1], d_edges[1:]), (d_hist, d_hist)))
        # Get centers of bins using edges
        bin_centers = (d_edges[:-1] + d_edges[1:]) / 2
        # Find maximum value (number of times mode occurred)
        peak_val = np.max(d_hist)
        # Get crude mode using this peak value
        mode = bin_centers[d_hist == peak_val][0]
        # Save all these findings to a dictionary for later use
        self.stats.update(dict(d_hist_edges=(d_hist, d_edges),
                               hist_xy=(hist_x, hist_y),
                               bin_centers=bin_centers,
                               peak_val=peak_val, peak_mode=mode))

    def calculate_offset(self):
        """
        Fit a Gaussian to the difference image histogram to refine the
        mode calculation and return that mode.
        Also saves the fitted parameters to the self.stats dictionary.
        :return: the fit Gaussian mean, representing the mode of the histogram
        """
        # Set up the Gaussian fit to the histogram
        peak_location, peak_val = self.stats['peak_mode'], self.stats['peak_val']
        # Default the standard deviation to 10; should be of that order
        p0 = [peak_location, 10, peak_val]
        x_to_fit = self.stats['bin_centers']
        y_to_fit = self.stats['d_hist_edges'][0]
        # We will try to fit above half-max to avoid bias from a noisy/uneven floor
        # Find full-width half-max locations
        try:
            # Use spline to find roots of function sunk by half-max
            spline = UnivariateSpline(x_to_fit, y_to_fit - peak_val/2, s=0)
            r1, r2 = spline.roots()
            mask_to_fit = (x_to_fit > r1) & (x_to_fit < r2)
        except ValueError:
            # Approximate roots by limiting function to above half max
            # Depending on function shape, this may not be ideal. Assumes Gaussian.
            print("spline fit estimate of histogram FWHM failed; approximating")
            mask_to_fit = y_to_fit > peak_val/2
        # noinspection SpellCheckingInspection,PyTypeChecker
        popt, pcov = curve_fit(gaussian, x_to_fit[mask_to_fit],
                               y_to_fit[mask_to_fit], p0=p0)
        # Save fitted Gaussian parameters to the stats dictionary
        self.stats['gauss_fit'] = popt
        # Return the mean of the Gaussian (more accurate mode of the distribution)
        mode = popt[0]
        return mode

    def diagnostic_difference_histogram(self):
        """
        Plot the difference image histogram with overlaid Gaussian fit.
        """
        mode, sigma, amplitude = self.stats['gauss_fit']
        plt.figure("Difference Histogram", figsize=(8, 5))
        # Plot the histogram itself
        plt.plot(*self.stats['hist_xy'], '-', color=(.1, .5, .1),
                 linewidth=3, label="$F_{GNILC} - F_{obs}$")
        curve_x = self.stats['hist_xy'][0]
        # Trim Gaussian fit x extent to 3 sigma to make it look better
        trim = 3 * sigma
        lo_trim, hi_trim = mode - trim, mode + trim
        # Plot the Gaussian fit
        curve_x = curve_x[(curve_x > lo_trim) & (curve_x < hi_trim)]
        curve_y = gaussian(curve_x, mode, sigma, amplitude)
        eq_string = r"$Ae^{-(x - \mu)^{2}/2\sigma^{2}}$"
        plt.plot(curve_x, curve_y, '-', color='k', alpha=.7,
                 label="{:s}->({:s}: {:.1f}, {:s}: {:.1f}, A: {:.1f})".format(
                     eq_string, r"$\mu$", mode,
                     r"$\sigma$", sigma,
                     amplitude))
        plt.plot([mode], [self.stats['peak_val']], 'x', markersize=12,
                 color='k', alpha=0.7, label="Offset: {:.2f} MJy/sr".format(mode))
        plt.title("Difference histogram: Predicted $-$ Observed")
        plt.xlabel("Pixel difference between predicted/observed (MJy/sr)")
        plt.ylabel("Histogram count")
        plt.legend()

    def diagnostic_flux_map(self):
        """
        Plot the predicted and observed flux.
        The observed flux is shown at native, not GNILC, resolution.
        """
        fig = plt.figure("Predicted & Observed Flux", figsize=(14, 7))
        # First plot predicted data
        ax = plt.subplot(121)
        # Find appropriate color limits
        fq, tq = visual_min_max(self.predicted_target_flux)
        last_plot = plt.imshow(self.predicted_target_flux, origin='lower',
                               vmin=fq, vmax=tq)
        fig.colorbar(last_plot, ax=ax)
        plt.title("GNILC Model Predicted Flux (GNILC resolution)")
        # Now repeat for observed data
        ax = plt.subplot(122)
        fq, tq = visual_min_max(self.target_data)
        last_plot = plt.imshow(self.target_data, origin='lower',
                               vmin=fq, vmax=tq)
        fig.colorbar(last_plot, ax=ax)
        plt.title("Observed {:s} flux (native resolution)".format(
            self.target_bandpass_stub))

    def diagnostic_mask(self):
        """
        Plot the mask used to calculate difference statistics
        """
        plt.figure("GNILC Trust Mask")
        plt.imshow(self.mask.astype(int), origin='lower',
                   vmin=0, vmax=1, cmap='Greys_r')
        plt.colorbar()
        plt.title("Mask for Difference Image; Included if True (1)")

    def get_offset(self, no_diagnostic=False, full_diagnostic=False):
        """
        Calculate and return the offset.
        Plot the histogram diagnostic unless no_diagnostic is True.
        Plot the flux and mask diagnostics if full_diagnostic is True
            (and no_diagnostic is False)
        :param no_diagnostic: True if you do NOT want ANYTHING to plot
        :param full_diagnostic: True if you want THREE plots
        :return: offset calibration needed by target, in MJy/sr
        """
        offset = self.calculate_offset()
        print("="*25)
        print("="*25)
        print("= OFFSET: {:6.2f} MJy/sr =".format(offset))
        print("="*25)
        print("="*25)
        if not no_diagnostic:
            self.diagnostic_difference_histogram()
            if full_diagnostic:
                self.diagnostic_flux_map()
                self.diagnostic_mask()
            plt.show()
        return offset


def visual_min_max(array):
    """
    Get locations of some visually appropriate min/max for an n-D array.
    Uses locations of 10th and 90th data percentiles
    :param array: n-dimensional array of values
    :return: tuple of the values near (10%, 90%) in sorted order
    """
    # Clean, flatten, and sort array
    array_1d = np.array(sorted(array[~np.isnan(array)].flatten()))  # FIXME need numpy sort
    # Index for quartiles. Values won't be exact, but will be close
    # enough if the array is large.
    first_quartile = array_1d[array_1d.size // 10]
    third_quartile = array_1d[9 * array_1d.size // 10]
    return first_quartile, third_quartile


def gaussian(x, mu, sigma, amplitude):
    coefficient = amplitude / (np.sqrt(2 * np.pi) * sigma)
    exponent = -((x - mu) ** 2 / (2 * sigma * sigma))
    return coefficient * np.exp(exponent)
