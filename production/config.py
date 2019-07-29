from scipy.constants import c
import numpy as np
from astropy.io.fits import open as fits_open
"""
Configuration module for file paths and related items

The items in this file should be changed by the user (probably Lee) when appropriate.
The items you will most likely change:

gnilc_directory
rimo_directory
herschel_bandpass_directory
"""

"""
***************************************************************************
***************************************************************************
GNILC COMPONENT MAPS

This is the path to the GNILC component maps.
Place them all in the same directory & put the path to the directory here.
***************************************************************************
***************************************************************************
"""
# *************************************************************************
gnilc_directory = ""
# *************************************************************************

"""
***************************************************************************
***************************************************************************
PLANCK HFI RIMO

Name/path to the Planck HFI Reduced Instrument Model (RIMO)
List the full path to the file here
***************************************************************************
***************************************************************************
"""
# *************************************************************************
rimo_directory = ""
# *************************************************************************

"""
***************************************************************************
***************************************************************************
HERSCHEL BANDPASS PROFILES

These are names/paths to the Herschel filter profiles.
At present time, these are the same ones Kevin is using;
    I copied them out of the manticore source code.
Place them all in the same directory & put the path to the directory here.
***************************************************************************
***************************************************************************
"""
# *************************************************************************
herschel_bandpass_directory = ""
# *************************************************************************

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Should not be any need to change anything below this (but who knows)
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# Bandpass stubs for each telescope (Planck HFI, Herschel PACS/SPIRE)
# These are the only filters this code supports
hfi_bandpass_stubs = ("F100", "F143", "F217", "F353", "F545", "F857")
herschel_bandpass_stubs = ("PACS160um", "SPIRE250um", "SPIRE350um", "SPIRE500um")
# FITS extension for each HFI band in the HFI RIMO
hfi_bandpass_indices = (3, 4, 5, 6, 7, 8)
# Bandpass centers for Herschel (wavelength, Angstroms)
herschel_bandpass_centers = (1539451.3, 2471245.1, 3467180.4, 4961067.7)


def if_bandpass_valid(bandpass_stub_list):
    """
    This is honestly just me having fun with decorators.
    :param bandpass_stub_list: List of valid bandpass "stub" names.
    :return: Decorator that limits its function to arguments within
        this bandpass_stub_list and throws an error if not.
    """
    def if_bandpass_valid_decorator(func_to_decorate):
        # Nested decorator so the original decorator can take arguments!
        def decorated_function(bandpass_stub):
            # Some function of a given bandpass, with the stub as the arg
            if bandpass_stub in bandpass_stub_list:
                return func_to_decorate(bandpass_stub)
            else:
                basic_message = "is not a valid bandpass stub."
                extra_info = f"Supported options are: {bandpass_stub_list}"
                raise RuntimeError(f"{bandpass_stub} {basic_message} {extra_info}")
        return decorated_function
    return if_bandpass_valid_decorator


def angstroms_to_hz(x_angstroms):
    """
    This is exactly what it looks like
    :param x_angstroms: wavelength in Angstroms (in vacuum)
    :return: frequency in Hz
    """
    return c / (x_angstroms * 1e-10)


def limit_planck_frequency(arr, f_arr):
    """
    Helper function for bandpass_profile functions
    Trims a frequency-related array using the frequency array
    :param arr: the array to trim
    :param f_arr: the frequency array
    :return: the trimmed array
    """
    mask = (f_arr > 1e10) & (f_arr < 2e12)
    return arr[mask]


class PlanckConfig:
    """
    Class to hold "static variables" describing locations of Planck-related files
    Doesn't need to be instanced; everything should be a class variable
    """

    # Unlikely that you need to change these
    gnilc_temp = gnilc_directory + "COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.00.fits"
    gnilc_tau = gnilc_directory + "COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.00.fits"
    gnilc_beta = gnilc_directory + "COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.00.fits"

    hfi_rimo = rimo_directory + "HFI_RIMO_R3.00.fits"

    @staticmethod
    @if_bandpass_valid(hfi_bandpass_stubs)
    def bandpass_profile(bandpass_stub):
        """
        Load bandpass profile from the HFI RIMO
        :param bandpass_stub: HFI band stub
        :return: arrays: frequency (Hz), weight (arbitrary)
        """
        with fits_open(PlanckConfig.hfi_rimo) as hdul:
            i = hfi_bandpass_indices[hfi_bandpass_stubs.index(bandpass_stub)]
            assert bandpass_stub == hdul[i].header['EXTNAME'][-4]
            frequency_hz = hdul[i].data['WAVENUMBER'] * c * 1e2
            weight = hdul[i].data['TRANSMISSION']
        weight = limit_planck_frequency(weight, frequency_hz)
        frequency_hz = limit_planck_frequency(frequency_hz, frequency_hz)
        return frequency_hz, weight

    @staticmethod
    @if_bandpass_valid(hfi_bandpass_stubs)
    def bandpass_center(bandpass_stub):
        """
        Calculate bandpass center for HFI bands
        :param bandpass_stub: HFI band stub
        :return: band center (Hz)
        """
        return float(bandpass_stub[1:]) * 1e9


class HerschelConfig:
    """
    Class to hold "static variables" describing locations of Herschel-related files
    Doesn't need to be instanced; everything should be a class variable
    """

    @staticmethod
    @if_bandpass_valid(herschel_bandpass_stubs)
    def bandpass_profile(bandpass_stub):
        """
        Load bandpass profile from saved versions of the Herschel profiles
        :param bandpass_stub: Herschel band stub
        :return: frequency (Hz), weight (arbitrary)
        """
        bandpass_filename = f"{herschel_bandpass_directory}{bandpass_stub}_fromManticore.dat"
        bandpass_data = np.loadtxt(bandpass_filename)
        frequency_hz, weight = bandpass_data[:, 0], bandpass_data[:, 1]
        return frequency_hz, weight

    @staticmethod
    @if_bandpass_valid(herschel_bandpass_stubs)
    def bandpass_center(bandpass_stub):
        """
        Calculate bandpass center for Herschel bands
        :param bandpass_stub: Herschel band stub
        :return: band center (Hz)
        """
        return angstroms_to_hz(herschel_bandpass_centers[herschel_bandpass_stubs.index(bandpass_stub)])

