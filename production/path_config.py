from scipy.constants import c
import numpy as np
from astropy.io.fits import open as fits_open
"""
Configuration module for file paths and related items, some file I/O

The items in this file should be changed by the user (probably Lee) when appropriate.
The items you will most likely change (at top of this file):

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
rimo_directory = "/home/ramsey/Documents/Research/Filaments/HFI_stuff/"
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
herschel_bandpass_directory = "/home/ramsey/Documents/Research/Filaments/Herschel_bands/"
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


def limit_planck_frequency(f_arr):
    """
    Helper function for bandpass_profile functions
    Returns a mask for trimming a frequency-related array
        based on the frequency array of the same shape
    Limits are generous and therefore hardcoded
    :param f_arr: the frequency array
    :return: the trimmed array
    """
    return (f_arr > 1e10) & (f_arr < 2e12)


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
            # Use extension index i to reference RIMO info
            i = hfi_bandpass_indices[hfi_bandpass_stubs.index(bandpass_stub)]
            # Sanity check!
            try:
                assert bandpass_stub == hdul[i].header['EXTNAME'][-4:]
            except AssertionError:
                msg = "Filter name mismatch. "
                msg += "Failed assertion between: "
                msg += f"{bandpass_stub} // {hdul[i].header['EXTNAME'][-4:]}"
                raise RuntimeError(msg)
            # Wavenumber delivered in cm-1; convert to frequency in Hz
            frequency_hz = hdul[i].data['WAVENUMBER'] * c * 1e2
            weight = hdul[i].data['TRANSMISSION']
        # Trim high and low frequencies
        f_mask = limit_planck_frequency(frequency_hz)
        weight = weight[f_mask]
        frequency_hz = frequency_hz[f_mask]
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


"""
Convenience functions exposed to the outside world.
"""


def get_bandpass_data(stub):
    """
    Get photometry weighting function for either Herschel or Planck HFI.
    Returns tuple(frequency array in Hz, weight array).
    Relies on this config module to have accurate filenames
        to all Herschel filter profiles as well as the Planck HFI RIMO.
    :param: stub: short string name indicating the bandpass filter
    :returns: tuple(array, array) of frequencies (Hz) and filter transmission.
        Normalization of the transmission curve is arbitrary
    """
    # Check if Planck ('F') or Herschel
    instrument = PlanckConfig if stub[0] == 'F' else HerschelConfig
    return instrument.bandpass_profile(stub)


def get_bandpass_center(stub):
    """
    Get effective band centers in Hz for either Herschel or Planck HFI
    Returns frequency in Hz
    Relies on this config module to have accurate filenames
        to all Herschel filter profiles as well as the Planck HFI RIMO
    :param stub: short string name indicating the bandpass filter
    :return: float frequency in Hz
    """
    # Check if Planck ('F') or Herschel
    instrument = PlanckConfig if stub[0] == 'F' else HerschelConfig
    return instrument.bandpass_center(stub)
