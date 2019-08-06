import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import compare_images as cimg
from math import ceil
from scipy.signal import convolve2d
import sys


per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
per1_dir = "../"

soln_5pcterr = "T4-absdiff-Per1J-3param-plus045-plus05.0pct-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-bcreinit-Th15.95-Nh5E19,2E22.fits"
soln_2p_5pcterr = "T4-absdiff-Per1J-plus045-plus05.0pct-pow-1000-0.1-1.80.fits"

manticore_nominal_2p = "T4-absdiff-Per1J-plus045-plus05.0pct-pow-1000-0.1-1.80-crop6.fits"
manticore_nominal_3p = "T4-absdiff-Per1J-3param-plus045-plus05.0pct-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-bcreinit-Th15.95-Nh5E19,2E22-crop6.fits"

def prepare_convolution(w, beam, n_sigma=3, method='scipy', data_shape=None):
    # Given a WCS object and beam FWHMs in arcminutes,
    #  returns the Gaussian needed to convolve image by this kernel
    # Gaussian is returned in smaller array that includes contributions out to n_sigma
    # Find pixel scale, in arcminutes
    # methods = 'scipy' or 'fft'
    dtheta_dpix_i = np.sqrt(np.sum([(x2 - x1)**2 for (x1, x2) in zip(w.wcs_pix2world(0, 0, 0), w.wcs_pix2world(0, 1, 0))]))*60
    dtheta_dpix_j = np.sqrt(np.sum([(x2 - x1)**2 for (x1, x2) in zip(w.wcs_pix2world(0, 0, 0), w.wcs_pix2world(1, 0, 0))]))*60
    dthetas = [dtheta_dpix_i, dtheta_dpix_j]
    dtheta_dpix_avg = (dtheta_dpix_i + dtheta_dpix_j)/2  # arcminutes
    sigma_arcmin = beam / 2.35 # FWHM to standard deviation
    ij_arrays = [None, None]
    for x in range(2):
        if method == 'scipy' or method == 'manual':
            several_sigma_pixels = int(ceil(n_sigma * sigma_arcmin / dthetas[x])) # number of pixels in 5 sigma
            kernel_width = 2*several_sigma_pixels + 1 # several_sigma_pixels on each side, plus center for zero
            x_array = np.arange(kernel_width).astype(float) - several_sigma_pixels
        elif method == 'fft':
            x_array = np.arange(data_shape[x]).astype(float) - data_shape[x]//2
        else:
            raise RuntimeError("Method {:s} is not implemented".format(method))
        x_array *= dthetas[x]
        y_array = np.exp(-x_array*x_array/(2*sigma_arcmin*sigma_arcmin))
        y_array = y_array / np.trapz(y_array)
        ij_arrays[x] = y_array
    i, j = ij_arrays
    convolution_beam = i[:, np.newaxis] * j[np.newaxis, :]
    return convolution_beam


def convolve_fft(image, kernel, **kwargs):
    # assumes no nans, just the convolution
    # kwargs are ignored;
    #   needs to match call signature of scipy.signal's convolve2d
    ft = np.fft.fft2(image)*np.fft.fft2(kernel)
    return np.real(np.fft.fftshift(np.fft.ifft2(ft)))

def convolve_properly(image, kernel, method='scipy'):
    # Preserve NaNs
    # also mitigate edge effects / normalization from NaN correction
    # offers both scipy and fft methods
    if method == 'scipy':
        convolve_routine = convolve2d
    elif method == 'fft':
        convolve_routine = convolve_fft
    else:
        raise RuntimeError("Method {:s} is not implemented".format(method))
    convkwargs = dict(mode='same', boundary='fill', fillvalue=0.)
    image = image.copy()
    nanmask = np.isnan(image)
    image[nanmask] = 0.
    result = convolve_routine(image, kernel, **convkwargs)
    # now account for edge effects / normalization
    image[~nanmask] = 1.
    norm = convolve_routine(image, kernel, **convkwargs)
    image[:] = 1.
    norm /= convolve_routine(image, kernel, **convkwargs)
    result /= norm
    result[nanmask] = np.nan
    return result

def convolve_from_source(coordinates_to_fill, source_mask, source_image, kernel,
    method='manual'):
    """
    Coordinates should be tuples(i coords array, j coords array)
    Assume that no coordinates will force kernel over the edge of the image
    Source mask should be float 1s and 0s, 1 where valid
    Source image should be 0 wherever it's invalid (no nans!)
    You can assume that the kernel side lengths are odd
    HEY there's two ways to do this
    1) manually loop thru points and find source contribution
    2) scipy convolve, normalize, and mask out the stuff that doesn't matter
    3) FFT convolve, normalize, and mask out the stuff that doesn't matter
    it's unclear which of these will actually take less time. we'll have to...
    # TODO: TRY BOTH AND COMPARE!
    """
    if method == 'manual':
        coordinates_to_fill = np.stack([*coordinates_to_fill], axis=1)
        kernel_extents = tuple(x // 2 for x in kernel.shape)
        values = []
        for ij in coordinates_to_fill:
            (i_lo, i_hi), (j_lo, j_hi) = tuple((i - kernel_extents[n], i + kernel_extents[n] + 1) for i, n in zip(ij, range(2)))
            predicted_value = np.sum(kernel * source_image[i_lo:i_hi, j_lo:j_hi]) / np.sum(kernel * source_mask[i_lo:i_hi, j_lo:j_hi])
            values.append(predicted_value)
        return np.array(values)
    else:
        return convolve_properly(source_image, kernel, method=method)[coordinates_to_fill] / convolve_properly(source_mask, kernel, method=method)[coordinates_to_fill]


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

def prepare_TN_maps(T, N, w, conv_sigma_mult=None, n_sigma=3, sigma_mult=1, method='manual'):
    """
    T is temperature map, N is column density map
    T, N are already convolved up a little bit
    w is WCS object from WCS(header)
    returns (T convolved, N convolved, convolution kernel)
    The convolution kernel can be used for inpainting
    """
    new_beam = cimg.bandpass_beam_sizes['SPIRE500um']
    if conv_sigma_mult != 'noconv':
        if conv_sigma_mult is None:
            conv_sigma_mult = 1./np.sqrt(2)
        conv_beam = prepare_convolution(w, new_beam*conv_sigma_mult, n_sigma=n_sigma, method='scipy')
        T = convolve_properly(T, conv_beam, method='scipy')
        N = convolve_properly(N, conv_beam, method='scipy')
    conv_beam = prepare_convolution(w, new_beam*sigma_mult, n_sigma=n_sigma, data_shape=T.shape, method=method)
    return T, N, conv_beam

def inpaint_mask(T, N):
    """
    T is temperature map, N is column density map
    T, N are already convolved up a little bit
    returns (ipmask, validmask)
    ipmask is true where data needs to be inpainted
    validmask is true where data is valid (not nan)
    """
    nanmask = np.isnan(T) & np.isnan(N)
    mask_N = (N > 1.5e21) # this is what we want to inpaint
    return mask_N, ~nanmask


def boolean_edges(mask, valid_mask, valid_borders=False):
    # takes in a bool array
    # returns bool array of 1s from mask that border 0s
    # will not allow 1s anywhere valid_mask is 0

    mshape = mask.shape

    if valid_borders:
        # return mask of pixels that border 8 or more valid_mask pixels
        mask = ~valid_mask
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
            exec(stmt, {}, {"dst": border_count, "src": mask})
    if valid_borders:
        return border_count < 8
    else:
        return border_count.astype(bool) & ~mask & valid_mask

def paint(paint_mask, valid_mask, source_image, kernel, method='scipy'):
    """
    paint_mask should be true where you want to paint
    valid_mask should be true at the initially valid source pixels
    source_image need not have the paint pixels removed; this is the
        image that the painting will be based on
    kernel should be small w.r.t. the source image. it need not be square,
        but each side length should be odd
    """
    work_in_progress = source_image.copy()
    source_mask = (valid_mask & ~paint_mask)
    while np.any(source_mask ^ valid_mask):
        # plt.imshow(source_mask ^ valid_mask, origin='lower')
        # plt.title("currently painting...")
        # plt.show()
        border_mask = boolean_edges(source_mask, valid_mask)
        border_coords = np.where(border_mask)
        border_values = convolve_from_source(border_coords,
            source_mask.astype(float), work_in_progress, kernel,
            method=method)
        work_in_progress[border_coords] = border_values
        source_mask[np.where(border_mask)] = True
    return work_in_progress


if __name__ == "__main__":
    with fitsopen(soln_2p_5pcterr) as hdul:
        T = hdul[1].data
        N = hdul[3].data
        w = WCS(hdul[1].header)
    # T = np.arange(50*50).astype(float).reshape(50, 50)
    # T = (T+8)*16/((50*50)+8)
    # T[20:30, 20:30] = 0
    # T[:4, :] = np.nan
    # T[-4:, :] = np.nan
    # T[:, :4] = np.nan
    # T[:, -4:] = np.nan
    # N = np.arange(50*50).astype(float).reshape(50, 50)
    # N[20:30, 20:30] = 2e22
    Torig, Norig = T.copy(), N.copy()
    method = 'manual'
    T, N, conv_kernel = prepare_TN_maps(Torig, Norig, w, n_sigma=3, conv_sigma_mult='noconv', sigma_mult=2, method=method)
    if method == 'scipy' or method == 'manual':
        print("KERNEL SHAPE", conv_kernel.shape)
        T = np.pad(T, tuple([x//2]*2 for x in conv_kernel.shape), mode='constant', constant_values=np.nan)
        N = np.pad(N, tuple([x//2]*2 for x in conv_kernel.shape), mode='constant', constant_values=np.nan)
    plot_it = True
    ipmask, validmask = inpaint_mask(T, N)
    source_mask = validmask & ~ipmask
    has_borders = boolean_edges(source_mask, validmask, valid_borders=True)
    ipmask[np.where(~has_borders)] = False
    # plt.figure()
    # plt.subplot(121)
    # plt.imshow(T, origin='lower')
    # plt.subplot(122)
    # plt.imshow(N, origin='lower')
    # fig, axes = plt.subplots(ncols=3, sharex=True, sharey=True)
    # axes[0].imshow(validmask, origin='lower')
    # axes[1].imshow(ipmask, origin='lower')
    # axes[2].imshow(ipmask, origin='lower')
    final_img = T
    final_img[np.where(ipmask | ~validmask)] = 0.
    final_img = paint(ipmask, validmask, final_img, conv_kernel, method=method)
    if method == 'scipy' or method == 'manual':
        print("ORIGINAL SHAPE", Torig.shape)
        final_img = final_img[(conv_kernel.shape[0]//2):(-conv_kernel.shape[0]//2 + 1), (conv_kernel.shape[1]//2):(-conv_kernel.shape[1]//2 + 1)]
        print("FINAL SHAPE", final_img.shape)
    print("DONE")
    final_img = prepare_TN_maps(final_img, N, w, conv_sigma_mult=3)[0]
    plt.figure(figsize=(14, 9))
    plt.imshow(final_img, origin='lower', vmin=14, vmax=18)
    plt.show()
