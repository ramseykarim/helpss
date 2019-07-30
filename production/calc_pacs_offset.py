import numpy as np
import matplotlib.pyplot as plt
from astropy.io.fits import getdata as fits_getdata
from astropy.wcs import WCS
from scipy.optimize import curve_fit

import utils_regrid as rgu
import path_config as cfg
import utils_planck as plu


class GNILCModel:

    def __init__(self, target_filename, target_bandpass='PACS160um', extension=0):
        target_data, target_head = fits_getdata(target_filename, ext=extension, header=True)
        band_dictionary_template = {'F545': None, 'F857': None}
        self.target_data = target_data
        self.target_wcs = WCS(target_head)
        self.target_bandpass_stub = target_bandpass
        self.projector = rgu.HEALPix2FITS(target_data, target_head, pixel_scale_arcsec=300)
        self.T = self.load_component('Temperature')
        self.beta = self.load_component('Spectral-Index')
        self.tau = self.load_component('Opacity')
        self.observed_planck_fluxes = band_dictionary_template.copy()
        self.predicted_planck_fluxes = band_dictionary_template.copy()
        self.ratios = band_dictionary_template.copy()
        self.masks = band_dictionary_template.copy()
        self.mask = None
        self.difference = None
        # self.accumulate_masks()
        # self.difference_to_target()

    @plu.shape_of_component
    def load_component(self, stub):
        component_filename = cfg.PlanckConfig.component_filename(stub)
        self.projector.healpix_to_intermediate(rgu.open_healpix(component_filename))
        return self.projector.pop_intermediate()

    @plu.shape_of_frequency
    def load_bandpass(self, stub):
        return cfg.get_bandpass_data(stub)

    def predict_flux_in_band(self, band_stub):
        frequency, transmission = self.load_bandpass(band_stub)
        band_center = cfg.get_bandpass_center(band_stub)
        image = plu.calculate_gnilc_flux(band_center, frequency, transmission,
                                         self.T, self.tau, self.beta)
        return self.projector.intermediate_to_target(intermediate=image)

    def generate_planck_ratio_mask(self, planck_band_stub):
        planck_map_fn = cfg.PlanckConfig.light_map_filename(planck_band_stub)
        observed_planck_flux = self.projector.project(rgu.open_healpix(planck_map_fn))
        observed_planck_flux = cfg.PlanckConfig.unit_conversion(planck_band_stub,
            observed_planck_flux)
        predicted_planck_flux = self.predict_flux_in_band(planck_band_stub)
        ratio = (observed_planck_flux / predicted_planck_flux)
        # Save these maps as diagnostics
        self.observed_planck_fluxes[planck_band_stub] = observed_planck_flux
        self.predicted_planck_fluxes[planck_band_stub] = predicted_planck_flux
        self.ratios[planck_band_stub] = ratio
        # Constrain ratio to be within 10%
        self.masks[planck_band_stub] = (ratio > 0.9) & (ratio < 1.1)

    def accumulate_masks(self):
        for band_stub in self.masks:
            self.generate_planck_ratio_mask(band_stub)
        self.mask = np.all(np.array(list(self.masks.values())), axis=0)

    def difference_to_target(self):
        # Assume this map to be at GNILC resolution
        predicted_target_flux = self.predict_flux_in_band(self.target_bandpass_stub)
        # Convolve Herschel map to Planck GNILC resolution
        h_beam_size = cfg.HerschelConfig.beam_size(self.target_bandpass_stub)
        p_beam_size = cfg.gnilc_resolution
        convolution_beam_size = np.sqrt(p_beam_size**2 - h_beam_size**2)
        convolution_beam = rgu.prepare_convolution(self.target_wcs,
                                                   convolution_beam_size,
                                                   self.target_data.shape)
        target_data_convolved = rgu.convolve_properly(self.target_data, convolution_beam)
        # [predicted] - [observed] difference map
        self.difference = predicted_target_flux - target_data_convolved

    def calculate_offset(self):
        dhist, dedges = np.histogram(self.difference[self.mask].ravel(),
                                     bins=128, range=(25, 80))
        bin_centers = (dedges[:-1] + dedges[1:]) / 2
        peak_val = np.max(dhist)
        mode = bin_centers[dhist == peak_val][0]
        p0 = [mode, 10, peak_val]
        popt, pcov = curve_fit(gaussian, bin_centers, dhist, p0=p0)
        mode, sigma, A = popt
        return mode


def gaussian(x, mu, sigma, amplitude):
    coefficient = amplitude / (np.sqrt(2 * np.pi) * sigma)
    exponent = -((x - mu) ** 2 / (2 * sigma * sigma))
    return coefficient * np.exp(exponent)


def calculate_pacs_offset(pacs_filename):
    model = GNILCModel(pacs_filename)
    return model.calculate_offset()

