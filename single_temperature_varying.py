import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import datetime

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u

from math import ceil
from scipy.signal import convolve2d

# OpenCV, computer vision library
import cv2

# This is mainly used to convolve and retrieve bandpass-specific info
# The pipeline/scripts code can do all this and is more specific and better
# documented. Should eventually switch over to that
import compare_images as cimg


"""
Created: Unsure!! would have to look back at git commits???
Edited August 19, 2020 to add/compare to OpenCV inpainting!
"""


# per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/" # jup
per1_dir = "/n/sgraraid/filaments/Perseus/Herschel/results/" # Jup
if not os.path.isdir(per1_dir):
    per1_dir = "../" # laptop

# IDK what these were for, but they'll fail if opened
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
    dtheta_dpix_i = w.array_index_to_world(0, 0).separation(w.array_index_to_world(0, 1)).to('arcmin').to_value()
    dtheta_dpix_j = w.array_index_to_world(0, 0).separation(w.array_index_to_world(1, 0)).to('arcmin').to_value()
    dthetas = [dtheta_dpix_i, dtheta_dpix_j]
    dtheta_dpix_avg = (dtheta_dpix_i + dtheta_dpix_j)/2  # arcminutes
    sigma_arcmin = beam / 2.35 # FWHM to standard deviation
    ij_arrays = [None, None]
    for x in range(2):
        if method == 'scipy' or method == 'manual':
            several_sigma_pixels = int(ceil(n_sigma * sigma_arcmin / dthetas[x])) # number of pixels in n sigma
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

def inpaint_mask(T, N, Ncutoff=1.5e21, Tcutoff=None):
    """
    T is temperature map, N is column density map
    T, N are already convolved up a little bit
    Tcutoff is not applied if it is None; otherwise, T > Tcutoff is included
        in the mask. Designed to eliminate extremely hot regions, SF (NGC 1333)
    returns (ipmask, validmask)
    ipmask is true where data needs to be inpainted
    validmask is true where data is valid (not nan)
    """
    nanmask = np.isnan(T) & np.isnan(N)
    main_mask = (N > Ncutoff) # this is what we want to inpaint
    if Tcutoff is not None:
        # Include the T > Tcutoff mask in this
        mask_T = (T > Tcutoff)
        main_mask |= mask_T
    return main_mask, ~nanmask


def display_inpaint_mask(mask):
    plt.figure(figsize=(14, 9))
    plt.title("Mask")
    plt.imshow(mask, origin='lower')

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


"""
INPAINTING METHODS
"""

def inpaint_byhand():
    """
    September 1, 2020: I am going to use a T < 20 filter to get rid of
        things like NGC 1333
        I am also going to use a lower N threshold, like 1e21
        The hope is that T inpainted is always 14 < T < 20
    """
    # Setup
    plot_kwargs = dict(origin='lower', vmin=13, vmax=18)
    fn_old = f"{per1_dir}full-1.5.3-Per1-pow-750-0.05625-1.75.fits"
    # Open files
    with fits.open(fn_old) as hdul:
        T = hdul['T'].data
        N = hdul['N(H2)'].data
        w = WCS(hdul[1].header)
    # Save originals
    Torig, Norig = T.copy(), N.copy()
    method = 'manual' # I think this is just about setting up the inpainting kernel
    # Don't convolve images yet; set up ~2-3x beam inpaint kernel
    # sigma_mult = (14./36.) * 8 # INPAINT KERNEL in Herschel-beams
    sigma_mult = 10. / cimg.bandpass_beam_sizes['SPIRE500um']
    T, N, conv_kernel = prepare_TN_maps(Torig, Norig, w, n_sigma=3, conv_sigma_mult='noconv', sigma_mult=sigma_mult, method=method)
    if method == 'scipy' or method == 'manual':
        print("KERNEL SHAPE", conv_kernel.shape)
        T = np.pad(T, tuple([x//2]*2 for x in conv_kernel.shape), mode='constant', constant_values=np.nan)
        N = np.pad(N, tuple([x//2]*2 for x in conv_kernel.shape), mode='constant', constant_values=np.nan)
    plot_it = True
    # Ncutoff = 7e21 # this was what I used with DL3, but I don't recall what that mask looked like
    Ncutoff = 1e21
    Tcutoff = 20.
    ipmask, validmask = inpaint_mask(T, N, Ncutoff=Ncutoff, Tcutoff=Tcutoff)
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


    fig = plt.figure(figsize=(18, 9))
    fig.canvas.set_window_title("Inpainting by hand")

    plt.subplot(221)
    plt.imshow(Torig, **plot_kwargs)
    plt.colorbar()

    plt.subplot(223)
    img_to_show = final_img.copy()
    img_to_show[img_to_show == 0] = np.nan
    plt.imshow(img_to_show, **plot_kwargs)
    plt.colorbar()

    final_img = paint(ipmask, validmask, final_img, conv_kernel, method=method)
    if method == 'scipy' or method == 'manual':
        print("ORIGINAL SHAPE", Torig.shape)
        final_img = final_img[(conv_kernel.shape[0]//2):(-conv_kernel.shape[0]//2 + 1), (conv_kernel.shape[1]//2):(-conv_kernel.shape[1]//2 + 1)]
        print("FINAL SHAPE", final_img.shape)
    print("DONE")
    final_img, nothing1, nothing2 = prepare_TN_maps(final_img, N, w, conv_sigma_mult=(2./np.sqrt(2)))
    final_img[final_img == 0] = np.nan

    final_img[final_img < 5] = Torig[final_img < 5]

    plt.subplot(122)
    plt.imshow(final_img, **plot_kwargs)
    plt.colorbar()


    # save the image in a manticore-ready format
    if False:
        with fits.open(fn_old) as hdul:
            Thdu = hdul[1]
            Xshdu = hdul[5]
            Thdu = hdul[0]
            Thdu.header['DATE'] = (datetime.datetime.now(datetime.timezone.utc).astimezone().isoformat(), "File creation date")
            Thdu.header['CREATOR'] = ("Ramsey: {}".format(str(__file__)), "FITS file creator")
            Thdu.header['HISTORY'] = "Ramsey inpainted this on Jan 13 2020."
            Thdu.header['HISTORY'] = "Cutoff: N={:4.0E}. Value from talks with LGM.".format(Ncutoff)
            Thdu.data = final_img
            Thdu.header['HISTORY'] = "Inpainted"
            Xshdu.header['HISTORY'] = "not changed, straight from 1-component manticore"
            hdulnew = fits.HDUList([Thdu, Thdu, Xshdu])
            hdulnew.writeto("/n/sgraraid/filaments/Perseus/Herschel/results/single-DL3-vary.fits", overwrite=True)
    elif True:
        with fits.open(fn_old) as hdul:
            Thdu = hdul[1]
        Thdu.header['DATE'] = (datetime.datetime.now(datetime.timezone.utc).astimezone().isoformat(), "File creation date")
        Thdu.header['EXTNAME'] = Thdu.header['EXTNAME']
        Thdu.header['BUNIT'] = Thdu.header['BUNIT']
        Thdu.header['CREATOR'] = ("Ramsey: {}".format(str(__file__)), "FITS file creator")
        Thdu.header['HISTORY'] = "Ramsey inpainted this on Sept 1 2020."
        Thdu.header['HISTORY'] = f"N Cutoff: N={Ncutoff:4.0E}. Value from talks with LGM."
        Thdu.header['HISTORY'] = f"T Cutoff: T={Tcutoff:.1f}. Prevents inpainting influenced by NGC 1333."
        Thdu.header['HISTORY'] = "Inpainting kernel had FWHM of 10 arcmin."
        Thdu.header['HISTORY'] = "Post-process convolution used 2/sqrt(2) 500um beams."
        Thdu.header['HISTORY'] = "Pixels < 5 K were replaced with the original temperature solution,"
        Thdu.header['HISTORY'] = " however, low-T edge effects still exist due to post-process convolution"
        Thdu.data = final_img
        Thdu.writeto(os.path.join(os.path.dirname(os.path.abspath(fn_old)), "single_T_inpainted_10arcmin.fits"), overwrite=True)

    plt.subplots_adjust(top=0.953,
        bottom=0.088,
        left=0.0,
        right=0.984,
        hspace=0.159,
        wspace=0.086)

    plt.show()


def inpaint_openCV(method=cv2.INPAINT_TELEA):
    """
    For descriptions of methods
    https://opencv-python-tutroals.readthedocs.io/en/latest/py_tutorials/py_photo/py_inpainting/py_inpainting.html
    For example code:
    https://github.com/opencv/opencv/blob/master/samples/python/inpaint.py
    """
    # Methods:
    # Alexandru Telea 2004 method, "Fast Marching"
    # Navier-Stokes

    # Now generic stuff
    plot_kwargs = dict(origin='lower', vmin=13, vmax=18)
    fn_old = f"{per1_dir}full-1.5.3-Per1-pow-750-0.05625-1.75.fits"
    with fits.open(fn_old) as hdul:
        T = hdul['T'].data
        N = hdul['N(H2)'].data
        w = WCS(hdul[1].header)
    Torig, Norig = T.copy(), N.copy()
    Ncutoff = 2.5e21
    ipmask, validmask = inpaint_mask(T, N, Ncutoff=Ncutoff)
    source_mask = validmask & ~ipmask
    has_borders = boolean_edges(source_mask, validmask, valid_borders=True)
    ipmask[np.where(~has_borders)] = False

    fig = plt.figure("OpenCV inpainting", figsize=(18, 9))
    plt.subplot(221)
    plt.imshow(Torig, **plot_kwargs)
    plt.colorbar()

    plt.subplot(223)
    img_to_show = Torig.copy()
    img_to_show[ipmask] = np.nan
    plt.imshow(img_to_show, **plot_kwargs)
    plt.colorbar()

    T = T.astype(np.float32) # openCV only works with floats as 32 bits
    ipmask_8 = ipmask.astype(np.uint8)
    inpaint_radius = 25 # pixels
    T_inpainted = cv2.inpaint(T, ipmask_8, inpaint_radius, method)

    plt.subplot(122)
    plt.imshow(T_inpainted, **plot_kwargs)
    plt.colorbar()

    plt.subplots_adjust(top=0.953,
        bottom=0.088,
        left=0.0,
        right=0.984,
        hspace=0.159,
        wspace=0.086)

    plt.show()


if __name__ == "__main__":
    inpaint_byhand()
