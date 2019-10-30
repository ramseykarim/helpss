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
# nominal_3p_fn = comp2_fns[0]

nominal_2p_fn = "../full-1.5-L723-pow-1000-0.1-1.80.fits"
nominal_3p_fn = "../full-1.5-L723-pow-1000-0.1-1.80,pow-1000-0.1-2.10-Th16.40.fits"

def masking_gridsamp():
    olkw = dict(origin='lower')
    dTkw = dict(**olkw, vmin=0, vmax=2)
    Nkw = dict(**olkw, vmin=1e21, vmax=3e21)
    with fits.open(nominal_2p_fn) as hdul:
        T = hdul[1].data
        N = hdul[3].data
        w = WCS(hdul[1].header)
    frames_to_get = {"Tc": 5, "dTc": 6, "band160": 20}
    with fits.open(nominal_3p_fn) as hdul:
        for k in frames_to_get:
            frames_to_get[k] = hdul[frames_to_get[k]].data
    frames_to_get.update({'T': T, 'N': N})
    plt.figure(figsize=(14, 7))
    plt.subplot(121)
    plt.imshow(frames_to_get['N'], **Nkw)
    frames_to_get['N'][np.isnan(frames_to_get['band160'])] = np.nan
    plt.subplot(122)
    plt.imshow(frames_to_get['N'], **Nkw)
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
