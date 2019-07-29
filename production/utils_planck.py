import numpy as np
from scipy.constants import k, h, c
from astropy.io import fits


def sizeof_fmt(num, suffix='B'):
    # Taken directly from stackoverflow
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)


pfreq_lims = (1e10, 2e12) # Hz
def limit_planck_frequency(arr, f_arr):
    mask = (f_arr > 1e10) & (f_arr < 2e12)
    return arr[mask]


def get_bandpass_data(stub):
    """
    Get photometry weighting function for either Herschel or Planck HFI.
    Returns tuple(frequency array in Hz, weight array).
    Relies on a bandpass_files dictionary to have accurate filenames
        to all Herschel filter profiles as well as the Planck HFI RIMO.
    :param: stub: short string name indicating the bandpass filter
    :returns: tuple(array, array) of frequencies (SI) and filter transmission.
        Normalization of the transmission curve is arbitrary
    """
    # Check if Planck ('F') or Herschel
    if stub[0] == "F":
        # Planck; use RIMO
        with fits.open(p_RIMO) as hdul:
            i = bandpass_files[stub]
            # Ensure this is the right band
            assert stub == hdul[i].header['EXTNAME'][-4:]
            frequency_SI = hdul[i].data['WAVENUMBER'] * c * 1e2
            weight = hdul[i].data['TRANSMISSION']
        # Limit the frequency range of the outputs to avoid large/small floats
        weight = limit_planck_frequency(weight, frequency_SI)
        #fixme why didn't this work? vvvv
        frequency_SI = frequency_SI[(frequency_SI > 1e10) & (frequency_SI < 2e12)] #limit_planck_frequency(frequency_SI, frequency_SI)
    else:
        # Herschel; use Kevin's filter profiles
        fn = bandpass_files[stub]
        bp_data = np.loadtxt(fn)
        frequency_SI, weight = bp_data[:, 0], bp_data[:, 1]
    return frequency_SI, weight
