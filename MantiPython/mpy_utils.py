import numpy as np
from astropy.io import fits
import scipy.constants as cst
from collections import deque

# Constants

# Kevin's values
cst.c = 2.99792458e8
H2mass = 1.6737237e-24 * 2.75 # g


# scipy value
# H2mass = 2.75 * cst.m_p * 1e3 * 1.001  # g (mass of H2 molecule in grams)

BINS = 32

meters_micron = 1e6
f_hz_meters = lambda hz: cst.c / hz
# so for NU_hz = frequency in hertz
# LAMBDA_micron = f_hz_meters(NU_hz) * meters_micron
# and NU_hz = f_hz_meters( LAMBDA_micron / meters_micron )
f_hz_micron = lambda hz: f_hz_meters(hz) * meters_micron

H_WL = [160, 250, 350, 500]
H_stubs = {
    160: "PACS160um",
    250: "SPIRE250um",
    350: "SPIRE350um",
    500: "SPIRE500um",
}

MANTICORE_KEYS = ["Tc", "dTc", "Nc", "dNc",
    "Th", "dTh", "Nh", "dNh",
    "obs160", "obs250", "obs350", "obs500",
    "err160", "err250", "err350", "err500",
    "mod160", "mod250", "mod350", "mod500"
]
# --------
# indices for all these are:
MANTICORE_INDICES_3p = [
    1, 2, 3, 4,
    5, 6, 7, 8,
    15, 17, 19, 21,
    16, 18, 20, 22,
    11, 12, 13, 14,
]
MANTICORE_INDICES_2p = [
    1, 2, 3, 4,
    None, None, None, None,
    11, 13, 15, 17,
    12, 14, 16, 18,
    7, 8, 9, 10,
]

PIXELS_OF_INTEREST = (
    (1, 'B1', (551-1, 307-1), 'dodgerblue'),
    (1, 'NGC 1333', (590-1, 520-1), 'midnightblue'),
    (1, 'L1455', (353-1, 554-1), 'forestgreen'),
    (0, 'L1455', (465-1, 513-1), 'olive'), # dNh larger than Nh
    (0, 'L1455', (425-1, 401-1), 'darkmagenta'),
    (0, 'B1', (587-1, 264-1), 'salmon'), # dNc ~ 0.5*Nc
    (0, 'B1', (587-1, 260-1), 'firebrick'), # Tc ~ 4 K
    (0, 'B1', (587, 260), 'deeppink'), # Tc ~ 6 K

)

def get_obs(info_dict):
    return [info_dict[x] for x in ("obs160", "obs250", "obs350", "obs500")]

def get_err(info_dict):
    return [info_dict[x] for x in ("err160", "err250", "err350", "err500")]

def get_mod(info_dict):
    return [info_dict[x] for x in ("mod160", "mod250", "mod350", "mod500")]


def get_manticore_info(source, *args):
    # i is the ROW, which is Y in FITS
    # zero indexed, so subtract 1 from FITS coordinates
    return_dict = dict()

    if isinstance(source, str):
        hdul = fits.open(source)
    else:
        hdul = source

    if len(args) == 2:
        mask = (args[0], args[1])
    elif isinstance(args[0], tuple):
        mask = tuple(zip(*args[0]))
    elif isinstance(args[0], np.ndarray):
        mask = args[0].astype(bool)
    else:
        msg = ["Unknown call signature ({:d} args):".format(len(args))]
        msg += [str(a) for a in args]
        raise RuntimeError(" ".join(msg))

    indices = MANTICORE_INDICES_3p if len(hdul) > 20 else MANTICORE_INDICES_2p

    for k, idx in zip(MANTICORE_KEYS, indices):
        if idx is not None:
            return_dict[k] = hdul[idx].data[mask]
    # but the last four (mod) need to be added to the obs
    for n in range(4):
        # obs is 8-11, so n+8
        # mod is 16-19, so n+16
        return_dict[MANTICORE_KEYS[n+16]] += return_dict[MANTICORE_KEYS[n+8]]
    if isinstance(source, str):
        hdul.close()
    return return_dict


prep_arr = lambda a, b: np.array([a, b]).T.flatten()
def histogram(x, x_lim=None):
    if x_lim is None:
        x_lim = (np.min(x), np.max(x))
    dhist, dedges = np.histogram(x.ravel(), bins=BINS, range=x_lim)
    histx, histy = prep_arr(dedges[:-1], dedges[1:]), prep_arr(dhist, dhist)
    return histx, histy
