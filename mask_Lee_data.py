plotting_remotely = False
import numpy as np
import matplotlib
if plotting_remotely:
    matplotlib.use('Agg')
SAVE_NAME = "/home/rkarim/Downloads/Figure_X_current.png"
import matplotlib.pyplot as plt
import glob, sys
from astropy.io import fits
from astropy.wcs import WCS
import compare_images as cimg
from math import ceil
from scipy.signal import convolve2d
from pipeline.calc_offset import flquantiles
import single_temperature_varying as ip

"""
Testing masks on Lee's data (small frames)
Creation date: Oct 12, 2019
"""
__author__ = "Ramsey Karim"


def show_plot():
	if plotting_remotely:
		plt.savefig(SAVE_NAME)
	else:
		plt.show()

# lee_dir = "/n/sgraraid/filaments/lgm/manticore-test/L723/"
# comp_stub = "TEST1-"
# comp1_dir = comp_stub + "1comp/"
# comp2_dir = comp_stub + "2comp*/"
# fits_stub = "*.fits"
# comp1_fn = glob.glob(lee_dir+comp1_dir+fits_stub)[0]
# comp2_fns = glob.glob(lee_dir+comp2_dir+fits_stub)
#
# nominal_2p_fn = comp1_fn
# nominal_3p_fn = next(x for x in comp2_fns if 'grid' in x)

nominal_2p_fn = "../full-1.5-L723-pow-1000-0.1-1.80.fits"
nominal_2p_fn_210 = "../full-1.5-L723-pow-1000-0.1-2.10.fits"
nominal_3p_fn = "../full-1.5-L723-pow-1000-0.1-1.80,pow-1000-0.1-2.10-Th16.40.fits"

BINS = 64
def histogram(x, x_lim=None):
    # yet another histogram routine
    # x should be flat
    if x_lim is None:
        x_lim = (np.min(x), np.max(x))
    dhist, dedges = np.histogram(x, bins=BINS, range=x_lim)
    prep_arr = lambda a, b: np.array([a, b]).T.flatten()
    histx, histy = prep_arr(dedges[:-1], dedges[1:]), prep_arr(dhist, dhist)
    bin_centers = (dedges[:-1]+dedges[1:])/2
    return {"hist": dhist, "edges": dedges, "bin_centers": bin_centers, "histxy": (histx, histy)}

def apply_nanmask(img, true_if_nan):
    img[true_if_nan] = np.nan

def flatten(img):
    return img[~np.isnan(img)].flatten()

def log10(x):
    return np.log10(x + 1e18)

def masking_gridsamp():
    olkw = dict(origin='lower')
    dTkw = dict(**olkw, vmin=0, vmax=2)
    Tkw = dict(**olkw, vmin=9, vmax=14)
    Nkw = dict(**olkw, vmin=1e21, vmax=3e21)
    with fits.open(nominal_2p_fn) as hdul:
        T = hdul[1].data
        N = hdul[3].data
        w = WCS(hdul[1].header)
    frames_to_get = {"Nh": 3, "Tc": 5, "dTc": 6, 'Nc': 7, "band160": 20}
    with fits.open(nominal_3p_fn) as hdul:
        for k in frames_to_get:
            frames_to_get[k] = hdul[frames_to_get[k]].data
    frames_to_get.update({'T': T, 'N': N})
    with fits.open(nominal_2p_fn_210) as hdul:
        f2g = {'N2.1':3, 'T2.1':1}
        for k in f2g:
            frames_to_get[k] = hdul[f2g[k]].data
    pacs_nanmask = np.isnan(frames_to_get['band160'])
    frames_flat = {}
    for x in frames_to_get:
        apply_nanmask(frames_to_get[x], pacs_nanmask)
        img_flat = flatten(frames_to_get[x])
        if 'N' in x:
            img_flat = np.log10(img_flat + 1e18)
        frames_flat[x] = img_flat

    def Nh_is_too_high_in_dense_regions():
        mask = (frames_to_get['Nh'] > 10**21.06) & (frames_to_get['N'] > 21.65)
        fig, axes = plt.subplots(ncols=2, sharex=True, sharey=True)
        axes[0].imshow(mask, **olkw)
        axes[1].imshow(frames_to_get['Nc'], **Nkw)

    def Nh_is_too_low_in_thin_regions():
        mask = (frames_to_get['Nh'] < 10**20.4)
        fig, axes = plt.subplots(ncols=2, sharex=True, sharey=True)
        axes[0].imshow(mask, **olkw)
        axes[1].imshow(frames_to_get['Nh'], **Nkw)

    def Nc_is_too_high_in_thin_regions():
        mask = (frames_to_get['Nc'] > frames_to_get['N'] * 10**0.2) & (frames_to_get['N'] < 10**21.43)
        fig, axes = plt.subplots(ncols=2, sharex=True, sharey=True)
        axes[0].imshow(mask, **olkw)
        axes[1].imshow(frames_to_get['Nc'], **Nkw)

    def N_scatter_plot():
        axScatter = plt.subplot(111)
        plt.scatter(frames_flat['N'], frames_flat['Nh'], marker='.', label='Nh')
        plt.scatter(frames_flat['N'], frames_flat['Nc'], marker='.', label='Nc')
        plt.scatter(frames_flat['N'], frames_flat['N2.1'], marker='.', label='N2.1')
        plt.legend()
        lolim = np.max([plt.xlim()[0], plt.ylim()[0]])
        hilim = np.min([plt.xlim()[1], plt.ylim()[1]])
        plt.plot([lolim, hilim], [lolim, hilim], '--', label='N=N')
        plt.plot([lolim, hilim], [lolim+0.2, hilim+0.2], '--', label='N=N+0.2')

    def T_scatter_plot():
        axScatter = plt.subplot(111)
        plt.scatter(frames_flat['T2.1'], frames_flat['T'], marker='.', label='T1.8', alpha=0.7)
        plt.scatter(frames_flat['T2.1'], frames_flat['Tc'], marker='.', label='Tc', alpha=0.7)
        plt.plot([np.min(frames_flat['T2.1']), np.max(frames_flat['T2.1'])], [16.4]*2, '--', label='Th')
        lolim = np.max([plt.xlim()[0], plt.ylim()[0]])
        hilim = np.min([plt.xlim()[1], plt.ylim()[1]])
        plt.plot([lolim, hilim], [lolim, hilim], '--', label='T=T')
        plt.legend()


    # Nc_is_too_high_in_thin_regions()
    T_scatter_plot()
    show_plot()

def inpainting():
    olkw = dict(origin='lower')
    plotT_kwargs = dict(origin='lower', vmin=14, vmax=18)
    plotN_kwargs = dict(origin='lower', vmin=1e21, vmax=2.5e21)
    with fits.open(nominal_2p_fn) as hdul:
        T = hdul[1].data
        N = hdul[3].data
        w = WCS(hdul[1].header)
    Torig, Norig = T.copy(), N.copy()
    method = 'manual'
    T, N, conv_kernel = ip.prepare_TN_maps(Torig, Norig, w,
        conv_sigma_mult='noconv', sigma_mult=2, method=method)
    if method == 'scipy' or method == 'manual':
        print("KERNEL SHAPE", conv_kernel.shape)
        T = np.pad(T, tuple([x//2]*2 for x in conv_kernel.shape), mode='constant', constant_values=np.nan)
        N = np.pad(N, tuple([x//2]*2 for x in conv_kernel.shape), mode='constant', constant_values=np.nan)

    validmask = ~(np.isnan(T) & np.isnan(N))
    ipmask = (N > 2.5e21)
    source_mask = validmask & ~ipmask
    has_borders = ip.boolean_edges(source_mask, validmask, valid_borders=True)
    # This cleans out lone pixels in a sea of NaNs
    ipmask[np.where(~has_borders)] = False # not a problem with this dataset
    final_img = T.copy()
    final_img[np.where(ipmask | ~validmask)] = 0
    plotting = False
    if plotting:
        fig = plt.figure(figsize=(16, 7))
        axes = tuple(plt.subplot(x) for x in (131, 132, 133))
        axes[0].imshow(Torig, **plotT_kwargs)
        axes[1].imshow(final_img, **plotT_kwargs)

    final_img = ip.paint(ipmask, validmask, final_img, conv_kernel, method=method)
    final_img[~validmask] = np.nan
    if method == 'scipy' or method == 'manual':
        print("ORIGINAL SHAPE", Torig.shape)
        print("PADDED SHAPE", final_img.shape)
        final_img = final_img[(conv_kernel.shape[0]//2):(-conv_kernel.shape[0]//2 + 1), (conv_kernel.shape[1]//2):(-conv_kernel.shape[1]//2 + 1)]
        print("FINAL SHAPE", final_img.shape)
    print("DONE")
    final_img = ip.prepare_TN_maps(final_img, N, w, conv_sigma_mult=(1./np.sqrt(2)))[0]
    print("NANMEDIAN", np.nanmedian(final_img))

    fitssavename = lee_dir + 'inpained_T.fits'

    with fits.open(nominal_2p_fn) as hdul:
        global_ext = fits.PrimaryHDU(header=hdul[0].header)
        T_header = hdul[1].header
        N_data, N_header = hdul[3].data, hdul[3].header
        Xs_data, Xs_header = hdul[5].data, hdul[5].header
        T_header['HISTORY'] = "Regions of high N(H2) removed, inpainted"
        T_ext = fits.ImageHDU(final_img, header=T_header)
        N_ext = fits.ImageHDU(N_data, header=N_header)
        Xs_ext = fits.ImageHDU(Xs_data, header=Xs_header)
        fits.HDUList([global_ext, T_ext, N_ext, Xs_ext]).writeto(
            fitssavename
        )

    if plotting:
        axes[2].imshow(final_img, **plotT_kwargs)
        show_plot()

def masking():
    with fits.open(nominal_2p_fn) as hdul:
        T2 = hdul[1].data
        N2 = hdul[3].data
        # pacs_mask True where VALID
        pacs_mask = ~np.isnan(hdul[12].data)

    with fits.open(nominal_3p_fn) as hdul:
        Nh = hdul[3].data
        Tc = hdul[5].data
        Nc = hdul[7].data

    def imglim(img):
        lims = flquantiles(img[~np.isnan(img)].ravel(), 10)
        return {v: l for v, l in zip(('vmin', 'vmax'), lims)}

    pltimg = Tc
    maskN = (N2 > 2.5e21) & (Nc > 0)
    lims = imglim(pltimg[pacs_mask&maskN])
    pltimg_c = pltimg.copy()
    pltimg_c[~(maskN&pacs_mask)] = np.nan
    """
    It seems like in single_temperature_varying.py,
    I used a cutoff of N2 = 1.5e21.
    This does not work here, since most of the frame (in L723)
    is above 1.5e21 to begin with.

    How about N2 > 2e21? This sort of works.
    It still lets in a lot of the stripey noise.

    I found that Nc == 0 for a lot of the stripey noise, so & Nc>0 helps.
    This still lets in some Tc~5K.

    How about examining the N2 map under these masks? Reveals N2>2.5e21
    could help.
    However, there is STILL Tc~5K. Maybe we need the Nh smoothing iteration.
    """

    plt.figure(figsize=(16, 9))
    plt.subplot(121)
    plt.imshow(pltimg, origin='lower', **lims)
    plt.colorbar()
    plt.subplot(122)
    plt.imshow(pltimg_c, origin='lower', **lims)
    plt.colorbar()
    plt.show()

if __name__ == "__main__":
    masking_gridsamp()
