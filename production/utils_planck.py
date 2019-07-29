import numpy as np
from scipy.constants import k, h, c
from astropy.io import fits


def sizeof_fmt(num, suffix='B'):
    """
    Convert number of bytes to string expression with proper binary prefix.
    Taken directly from StackOverflow
    :param num: size of some object in number of bytes
    :param suffix: symbol for basic unit. Default='B' for Bytes
    :return: formatted, human-readable string describing the size
    """
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)


planck_freq_lim = (1e10, 2e12)  # Hz


def get_bandpass_data(stub):
    """
    Get photometry weighting function for either Herschel or Planck HFI.
    Returns tuple(frequency array in Hz, weight array).
    Relies on a bandpass_files dictionary to have accurate filenames
        to all Herschel filter profiles as well as the Planck HFI RIMO.
    :param: stub: short string name indicating the bandpass filter
    :returns: tuple(array, array) of frequencies (Hz) and filter transmission.
        Normalization of the transmission curve is arbitrary
    """
    # Check if Planck ('F') or Herschel
    if stub[0] == "F":
        # Planck; use RIMO
        with fits.open(p_RIMO) as hdul:
            i = bandpass_files[stub]
            # Ensure this is the right band
            assert stub == hdul[i].header['EXTNAME'][-4:]
            frequency_hz = hdul[i].data['WAVENUMBER'] * c * 1e2
            weight = hdul[i].data['TRANSMISSION']
        # Limit the frequency range of the outputs to avoid large/small floats
        weight = limit_planck_frequency(weight, frequency_hz)
        # fixme why didn't this work? vvvv
        # limit_planck_frequency(frequency_hz, frequency_hz)
        frequency_hz = frequency_hz[(frequency_hz > 1e10) & (frequency_hz < 2e12)]
    else:
        # Herschel; use Kevin's filter profiles
        fn = bandpass_files[stub]
        bp_data = np.loadtxt(fn)
        frequency_hz, weight = bp_data[:, 0], bp_data[:, 1]
    return frequency_hz, weight
