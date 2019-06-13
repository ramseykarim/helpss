import numpy as np
from astropy.io import fits
import scipy.constants as cst
from collections import deque
import emcee
import corner
import pickle
import matplotlib.pyplot as plt

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
    "mod160", "mod250", "mod350", "mod500",
    "chi_sq"
]
# --------
# indices for all these are:
MANTICORE_INDICES_3p = [
    1, 2, 3, 4,
    5, 6, 7, 8,
    15, 17, 19, 21,
    16, 18, 20, 22,
    11, 12, 13, 14,
    9,
]
MANTICORE_INDICES_2p = [
    1, 2, 3, 4,
    None, None, None, None,
    11, 13, 15, 17,
    12, 14, 16, 18,
    7, 8, 9, 10,
    5,
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

PIXEL_CHAIN = (
    (554-1, 271-1),
    (588-1, 305-1)
)

def gen_CHAIN_dict(source):
    start, end = PIXEL_CHAIN
    coords = []
    current_coord = start
    for i in range(35):
        coords.append(current_coord)
        current_coord = tuple(x+1 for x in current_coord)
    return get_manticore_info(source, tuple(coords))


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

def emcee_3p(index, info_dict,
        dust=None, instrument=None, goodnessoffit=None,
        niter=800, burn=400, nwalkers=60):
    ndim = 3
    p_labels = ('Tc', 'Nh', 'Nc')
    nominal = [info_dict[x][index] for x in p_labels]
    Th = info_dict['Th'][index]
    for i in (1, 2):
        nominal[i] = np.log10(nominal[i])
    if dust is None:
        raise RuntimeError("dust!")
    if not isinstance(dust, list):
        dust = dust()
    if instrument is None:
        raise RuntimeError("instrument!")
    if not isinstance(instrument, list):
        instrument = instrument()
    if goodnessoffit is None:
        raise RuntimeError("goodnessoffit!")
    obs = [x[index] for x in get_obs(info_dict)]
    err = [x[index] for x in get_err(info_dict)]
    dof = 1.
    arguments = (dust, obs, err, instrument, Th, dof)
    Tlo, Thi, Nlo, Nhi = 0., Th, 17., 26.
    Tcenter, Ncenter = 10, 21
    def lnposterior(x):
        T, N1, N2 = x
        if T<Tlo or T>Th or N1<Nlo or N1>Nhi or N2<Nlo or N2>Nhi:
            return -np.inf
        else:
            return -1*goodnessoffit(x, *arguments)
    p0 = np.concatenate([
        np.random.normal(scale=3, size=(nwalkers, 1)) + Tcenter,
        np.random.normal(scale=1.5, size=(nwalkers, 2)) + Ncenter
    ], axis=1)
    badTmask = ((p0[:, 0]<Tlo)|(p0[:, 0]>Thi))
    p0[badTmask, 0] = np.random.normal(scale=.5, size=p0[badTmask, 0].shape) + Tcenter
    badNmask = ((p0[:, 1]<Nlo)|(p0[:, 1]>Nhi))
    p0[badNmask, 1] = np.random.normal(scale=.3, size=p0[badNmask, 1].shape) + Ncenter
    badNmask = ((p0[:, 2]<Nlo)|(p0[:, 2]>Nhi))
    p0[badNmask, 2] = np.random.normal(scale=.3, size=p0[badNmask, 2].shape) + Ncenter
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnposterior)
    sampler.run_mcmc(p0, niter+burn)
    samples = sampler.chain[:, burn:, :].reshape((-1, ndim))
    fig = corner.corner(samples, labels=p_labels, truths=nominal,
        range=[(0, 16), (19.5, 21.6), (17., 23.5)],)
    fig.set_size_inches((10, 10))
    plt.title("pixel #{:02d}, n={:d}".format(index, niter*nwalkers))
    plt.savefig("./corner1_{:02d}.pdf".format(index))
    # with open("./emcee_imgs/samples1_{:02d}.pkl".format(index), 'wb') as pfl:
    #     pickle.dump(samples, pfl)

def grid_3d(index, info_dict,
    dust=None, instrument=None, goodnessoffit=None,
    Tcgrid=None, Nhgrid=None, Ncgrid=None,
    empty_grid=None):
    p_labels = ('Tc', 'Nh', 'Nc')
    nominal = [info_dict[x][index] for x in p_labels]
    Th = info_dict['Th'][index]
    for i in (1, 2):
        nominal[i] = np.log10(nominal[i])
    if dust is None:
        raise RuntimeError("dust!")
    if not isinstance(dust, list):
        dust = dust()
    if instrument is None:
        raise RuntimeError("instrument!")
    if not isinstance(instrument, list):
        instrument = instrument()
    if goodnessoffit is None:
        raise RuntimeError("goodnessoffit!")
    obs = [x[index] for x in get_obs(info_dict)]
    err = [x[index] for x in get_err(info_dict)]
    dof = 1.
    arguments = (dust, obs, err, instrument, Th, dof)
    own_grid = (empty_grid is None)
    if own_grid:
        Tclim, Nclim = (0, 16), (19, 22.5)
        Nhlim = (20, 21.5)
        Tcrange = np.arange(*Tclim, 0.1)
        Nhrange = np.arange(*Nhlim, 0.05)
        Ncrange = np.arange(*Nclim, 0.05)
        Tcgrid, Nhgrid, Ncgrid = np.meshgrid(Tcrange, Nhrange, Ncrange, indexing='ij')
        empty_grid = np.empty(Tcgrid.size)
    for i, pvec in enumerate(zip(Tcgrid.ravel(), Nhgrid.ravel(), Ncgrid.ravel())):
        gof = goodnessoffit(pvec, *arguments)
        empty_grid[i] = gof
    empty_grid = np.log10(empty_grid.reshape(Tcgrid.shape))
    with open("./emcee_imgs/grid1_{:02d}.pkl".format(index), 'wb') as pfl:
        pickle.dump(empty_grid, pfl)
    return empty_grid

def render_grid(index, Tclim=None, Nhlim=None, Nclim=None):
    if Tclim is None:
        Tclim, Nclim = (0, 16), (19, 22.5)
        Nhlim = (20, 21.5)
    return
