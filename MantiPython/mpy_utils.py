import numpy as np
from astropy.io import fits
import scipy.constants as cst
from collections import deque
import emcee
import corner
import pickle
import matplotlib.pyplot as plt
import sys

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

LIMS_hidef_00 = (
    (11.645157051086425, 13.245157051086426,), # Tc
    (20.694322967529295, 20.8943229675293,), # Nh
    (20.65826873779297, 21.058268737792968,), # Nc
    0.01, 0.0025) # dT, dN (arange)

LIMS_grid1 = ((0, 16.), (20, 21.5), (19, 22.5), 0.1, 0.05) # Tc,Nh,Nc,dT,dN (arange)
LIMS_grid2 = ((4, 16.01), (20, 21.21), (19, 22.41), 0.05, 0.025) # same ^

def genranges(lims, differentials):
    # lims: Tc, Nh, Nc
    dT, dN = differentials
    return [np.arange(*l, d) for l, d in zip(lims, (dT, dN, dN))]

def gengrids(ranges):
    # output of genranges
    return np.meshgrid(*ranges, indexing='ij')

P_LABELS = ('Tc', 'Nh', 'Nc')
PE_LABELS = ('dTc', 'dNh', 'dNc')

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
    empty_grid=None, fname_override=None):
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
        dT, dN = 0.1, 0.05
        Tcrange, Nhrange, Ncrange = genranges((Tclim, Nhlim, Nclim), (dT, dN))
        Tcgrid, Nhgrid, Ncgrid = gengrids(Tcrange, Nhrange, Ncrange)
        empty_grid = np.empty(Tcgrid.size)
    for i, pvec in enumerate(zip(Tcgrid.ravel(), Nhgrid.ravel(), Ncgrid.ravel())):
        gof = goodnessoffit(pvec, *arguments)
        empty_grid[i] = gof
    empty_grid = np.log10(empty_grid.reshape(Tcgrid.shape))
    if fname_override is not None:
        fname = fname_override
    else:
        fname = "./emcee_imgs/grid1_{:02d}.pkl".format(index)
    with open(fname, 'wb') as pfl:
        pickle.dump(empty_grid, pfl)
    return empty_grid


def render_grid(index, info_dict, fname=None, savename=None,
    gofgrid=None, grids=None, ranges=None,
    fig_size=(1200, 1050), Tscale=4,
    more_contours=False, point_size=0.2, focalpoint_nominal=False,
    mlab=None):
    # mayavi rendering of some contours over chi squared surface
    if mlab is None:
        raise RuntimeError("Please pass the mlab module as kwarg 'mlab'")
    if grids is None:
        raise RuntimeError("Can't just leave grids blank these days")
    Tcgrid, Nhgrid, Ncgrid = grids
    if ranges is None:
        raise RuntimeError("Can't just leave ranges blank these days")
    Tcrange, Nhrange, Ncrange = ranges
    if (fname is None) and (gofgrid is None):
        # gofgrid, if given, is assumed to be multiplied by -1 and un-logged already
        raise RuntimeError("Give filename of grid pickle or the actual grid")
    nominal = [info_dict[x][index] for x in P_LABELS]
    for x in (1, 2):
        nominal[x] = np.log10(nominal[x])
    chi_sq = info_dict['chi_sq'][index]
    if gofgrid is None:
        with open(fname, 'rb') as pfl:
            gofgrid = -1*(10**pickle.load(pfl))
    fig = mlab.figure(figure="main", size=fig_size)
    src = mlab.pipeline.scalar_field(Tcgrid/Tscale, Nhgrid, Ncgrid, gofgrid)
    if more_contours:
        ####### Grey Xs=100 contour
        mlab.pipeline.iso_surface(src, contours=[-100,],
            opacity=0.1, vmin=-101, vmax=-100, colormap='gist_yarg')
        ####### Blue Xs~5 contours
        mlab.pipeline.iso_surface(src, contours=[-5, -3],
            colormap='cool', opacity=0.15, vmin=-8, vmax=-2)
    ####### Red/Yellow Xs~1 contours
    mlab.pipeline.iso_surface(src, contours=[-1.5, -1, -.5],
        colormap='hot', opacity=0.25, vmin=-2, vmax=-.3)
    ####### Axes
    mlab.axes(ranges=sum(([x.min(), x.max()] for x in (Tcrange, Nhrange, Ncrange)), []),
        extent=sum(([x.min(), x.max()] for x in (Tcrange/Tscale, Nhrange, Ncrange)), []),
        nb_labels=5, xlabel="Tc", ylabel='Nh', zlabel="Nc")
    ####### Title
    mlab.title("pt({:02d}) [Tc: {:04.1f}, Nh: {:4.1f}, Nc: {:4.1f}], ChiSq: {:6.2f}".format(
        index, *nominal, chi_sq),
        size=0.25, height=.9)
    nominal[0] /= Tscale # since the grid is rescaled
    ####### Manticore solution point
    mlab.points3d(*([x] for x in nominal),
        colormap='flag', mode='axes', scale_factor=point_size, line_width=4)
    for x in ("x", "y", "z"):
        pts = mlab.points3d(*([x] for x in nominal),
            colormap='flag', mode='axes', scale_factor=point_size, line_width=4)
        eval("pts.glyph.glyph_source._trfm.transform.rotate_{:s}(45)".format(x))
    ####### Favorable camera angle
    if focalpoint_nominal:
        focalpoint = nominal
    else:
        focalpoint = [10./Tscale, 20.75, 20.75]
    mlab.view(azimuth=45., elevation=92., distance=9.,
        focalpoint=focalpoint)
    if savename is not None:
        mlab.savefig(savename, figure=fig, magnification=1)
        mlab.clf()
    else:
        mlab.show()
    return
