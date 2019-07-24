import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import compare_images as cimg
from math import ceil
from scipy.signal import convolve2d


per1_dir = "../"
per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"

soln_5pcterr = "T4-absdiff-Per1J-3param-plus045-plus05.0pct-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-bcreinit-Th15.95-Nh5E19,2E22.fits"
soln_2p_5pcterr = "T4-absdiff-Per1J-plus045-plus05.0pct-pow-1000-0.1-1.80.fits"

manticore_nominal_2p = "T4-absdiff-Per1J-plus045-plus05.0pct-pow-1000-0.1-1.80-crop6.fits"
manticore_nominal_3p = "T4-absdiff-Per1J-3param-plus045-plus05.0pct-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-bcreinit-Th15.95-Nh5E19,2E22-crop6.fits"

def prepare_convolution(w, beam, n_sigma=3):
    # Given a WCS object and beam FWHMs in arcminutes,
    #  returns the Gaussian needed to convolve image by this kernel
    # Gaussian is returned in smaller array that includes contributions out to 5sigma
    # Find pixel scale, in arcminutes
    dtheta_dpix_i = np.sqrt(np.sum([(x2 - x1)**2 for (x1, x2) in zip(w.wcs_pix2world(0, 0, 0), w.wcs_pix2world(0, 1, 0))]))*60
    dtheta_dpix_j = np.sqrt(np.sum([(x2 - x1)**2 for (x1, x2) in zip(w.wcs_pix2world(0, 0, 0), w.wcs_pix2world(1, 0, 0))]))*60
    dthetas = [dtheta_dpix_i, dtheta_dpix_j]
    dtheta_dpix_avg = (dtheta_dpix_i + dtheta_dpix_j)/2  # arcminutes
    sigma_arcmin = beam / 2.35 # FWHM to standard deviation
    ij_arrays = [None, None]
    for x in range(2):
        several_sigma_pixels = int(ceil(n_sigma * sigma_arcmin / dthetas[x])) # number of pixels in 5 sigma
        kernel_width = 2*several_sigma_pixels + 1 # several_sigma_pixels on each side, plus center for zero
        x_array = np.arange(kernel_width).astype(float) - several_sigma_pixels
        x_array *= dthetas[x]
        print(x_array)
        y_array = np.exp(-x_array*x_array/(2*sigma_arcmin*sigma_arcmin))
        y_array = y_array / np.trapz(y_array)
        ij_arrays[x] = y_array
    i, j = ij_arrays
    convolution_beam = i[:, np.newaxis] * j[np.newaxis, :]
    return convolution_beam


def convolve_properly(image, kernel):
    # Preserve NaNs
    # also mitigate edge effects / normalization from NaN correction
    # try to use scipy
    convkwargs = dict(mode='same', boundary='fill', fillvalue=0.)
    image = image.copy()
    nanmask = np.isnan(image)
    image[nanmask] = 0.
    result = convolve2d(image, kernel, **convkwargs)
    # now account for edge effects / normalization
    image[~nanmask] = 1.
    norm = convolve2d(image, kernel, **convkwargs)
    image[:] = 1.
    norm /= convolve2d(image, kernel, **convkwargs)
    result /= norm
    result[nanmask] = np.nan
    return result


def fitsopen(fn):
    return fits.open(per1_dir+fn)


def get_T_wcs(fn):
    with fitsopen(fn) as hdul:
        T = hdul[1].data
        w = WCS(hdul[1].header)
    return T, w


def test_compare_convolution_algs():
    T, w = get_T_wcs(soln_2p_5pcterr)
    new_beam = cimg.bandpass_beam_sizes['SPIRE500um'] / np.sqrt(2)
    # FT (homegrown)
    conv_beam_ft = cimg.prepare_convolution(w, 0, new_beam, T.shape)
    convolved_img_ft = cimg.convolve_properly(T, conv_beam_ft)
    # scipy (probably more trustworthy)
    conv_beam = prepare_convolution(w, new_beam)
    convolved_img = convolve_properly(T, conv_beam)
    fig, axes = plt.subplots()
    plt.imshow(convolved_img - convolved_img_ft, origin='lower', vmin=-0.3, vmax=0.3)
    plt.colorbar()
    plt.show()


def test_convolve_T():
    T, w = get_T_wcs(soln_2p_5pcterr)
    new_beam = cimg.bandpass_beam_sizes['SPIRE500um'] / np.sqrt(2)
    conv_beam = prepare_convolution(w, new_beam)
    convolved_img = convolve_properly(T, conv_beam)    
    fig, axes = plt.subplots(ncols=2, nrows=1, sharex=True, sharey=True)
    p = axes[1].imshow(convolved_img, origin='lower', vmin=13, vmax=19)
    fig.colorbar(p, ax=axes[1], orientation='horizontal')
    p = axes[0].imshow(T, origin='lower', vmin=13, vmax=19)
    fig.colorbar(p, ax=axes[0], orientation='horizontal')
    plt.show()


def test_convolve_TNXs():
    with fitsopen(soln_2p_5pcterr) as hdul:
        T = hdul[1].data
        N = hdul[3].data
        Xs = hdul[5].data
        w = WCS(hdul[1].header)
    new_beam = cimg.bandpass_beam_sizes['SPIRE500um'] / np.sqrt(2)
    conv_beam = prepare_convolution(w, new_beam)
    fig, axes = plt.subplots(ncols=3, nrows=2, sharex=True, sharey=True)
    imgs = (T, N, Xs)
    lims = ((13, 19), (0, 4e21), (0, 2))
    axkwarg = dict(orientation='horizontal')
    for i in range(len(imgs)):
        img, lim = imgs[i], lims[i]
        pkwargs = dict(vmin=lim[0], vmax=lim[1], origin='lower')
        convolved_img = convolve_properly(img, conv_beam)
        p = axes[0, i].imshow(img, **pkwargs)
        fig.colorbar(p, ax=axes[0, i], **axkwarg)
        p = axes[1, i].imshow(convolved_img, **pkwargs)
        fig.colorbar(p, ax=axes[1, i], **axkwarg)
    plt.show()


def test_filter_conv_N():
    with fitsopen(soln_2p_5pcterr) as hdul:
        T = hdul[1].data
        N = hdul[3].data
        w = WCS(hdul[1].header)
    new_beam = cimg.bandpass_beam_sizes['SPIRE500um'] / np.sqrt(2)
    conv_beam = prepare_convolution(w, new_beam)
    T_conv = convolve_properly(T, conv_beam)
    N_conv = convolve_properly(N, conv_beam)
    mask_N = (N_conv < 1.5e21)
    T_conv[~mask_N] = np.nan
    axkwarg = dict(orientation='vertical')
    pkwarg = dict(origin='lower')
    fig, axes = plt.subplots(ncols=2, nrows=2, sharex=True, sharey=True)
    # plot T
    lims = dict(vmin=13, vmax=19)
    p = axes[0, 0].imshow(T, **pkwarg, **lims)
    fig.colorbar(p, ax=axes[0, 0], **axkwarg)
    p = axes[1, 0].imshow(T_conv, **pkwarg, **lims)
    fig.colorbar(p, ax=axes[1, 0], **axkwarg)
    # plot N
    lims = dict(vmin=0, vmax=4e21)
    p = axes[0, 1].imshow(N, **pkwarg, **lims)
    fig.colorbar(p, ax=axes[0, 1], **axkwarg)
    p = axes[1, 1].imshow(N_conv, **pkwarg, **lims)
    fig.colorbar(p, ax=axes[1, 1], **axkwarg)
    plt.show()


def inpaint_mask():
    # returns (ipmask, validmask)
    # ipmask is true where data needs to be inpainted
    # validmask is true where data is valid (not nan)
    with fitsopen(soln_2p_5pcterr) as hdul:
        T = hdul[1].data
        N = hdul[3].data
        w = WCS(hdul[1].header)
    nanmask = np.isnan(T)
    new_beam = cimg.bandpass_beam_sizes['SPIRE500um'] / np.sqrt(2)
    conv_beam = prepare_convolution(w, new_beam)
    T_conv = convolve_properly(T, conv_beam)
    N_conv = convolve_properly(N, conv_beam)
    mask_N = (N_conv < 1.5e21)
    return ~mask_N, ~nanmask
    T_conv[~mask_N] = np.nan
    plt.imshow(T_conv, origin='lower', vmin=13, vmax=19)
    plt.show()
    # take the masked T_conv map and fill in the internal nans
    # don't fill in the nans outside (nanmask)
    pass


def boolean_edges(mask):
    # takes in a bool array
    # returns bool array of 1s from mask that border 0s

    mshape = mask.shape
    # cast to int so we can add
    mask = mask.astype(int)
    border_count = np.zeros(mshape, dtype=int)
    """
    [:, :]     [:, :] # self
    [:, 1:]    [:, :-1] # right
    [:, -1]    [:, 1:] # left
    [1:, :]    [:-1, :] # up
    [1:, 1:]   [:-1, :-1] # top right
    [1:, :-1]  [:-1, 1:] # top left
    [:-1, :]   [1:, :] # down
    [:-1, 1:]  [1:, :-1] # bottom right
    [:-1, :-1] [1:, 1:] # bottom left
    """
    options = ((":", ":"), ("1:", ":-1"), (":-1", "1:"))
    for i in range(3):
        row = options[i]
        for j in range(3):
            col = options[j]
            stmt_l = "[{:s}, {:s}]".format(row[0], col[0])
            stmt_r = "[{:s}, {:s}]".format(row[1], col[1])
            stmt = "dst{:s} = dst{:s} + src{:s}".format(stmt_l, stmt_l, stmt_r)
            print(stmt)
            exec(stmt, {}, {"dst": border_count, "src": mask})
    # at this point, can add 
    return border_count.astype(bool) ^ mask


if __name__ == "__main__":
    ipmask, validmask = inpaint_mask()
    notipmask = ~(validmask&ipmask)
    border = boolean_edges(notipmask)
    plt.imshow(border, origin='lower')
    # arr = np.arange(64).reshape(8, 8)
    # mask = np.ones(arr.shape)
    # mask[arr > 47] = 0
    # mask[arr < 16] = 0
    # mask[arr%8 < 1] = 0
    # mask[arr%8 > 6] = 0
    # mask = ~(mask.astype(bool))
    # count = boolean_edges(mask)
    # plt.imshow(count, origin='lower')
    plt.show()
