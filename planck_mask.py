import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
import sys
from functools import reduce
from compare_images import prepare_convolution, convolve_properly, plot_compare, gaussian

XsMAX, XsMIN = 1.e0, 1.e-3
XsMOD = 2
XsMASK = lambda x: (x < XsMAX) & (x > XsMIN)
BINS = 128
diffLIM = (18, 68)
small_diffLIM = (-10, 10)
fluxLIM = (0, 100)
XsPLOTLIM = (-5, 3.1)

CURRENT_DIR = "/n/sgraraid/filaments/data/TEST4/helpss_scratch_work/"
SAVE_NAME = "Figure_X_current.png"

def show_plot():
    if plotting_remotely:
        plt.savefig(SAVE_NAME)
    else:
        plt.show()

bandpass_beam_sizes = {
    "PACS160um": np.sqrt(11.64 * 15.65)/60,
    "SPIRE250um": 18.4/60,
    "SPIRE350um": 25.2/60,
    "SPIRE500um": 36.7/60,
    "F100": 9.66,
    "F143": 7.27,
    "F217": 5.01,
    "F353": 4.86,
    "F545": 4.84,
    "F857": 4.63,
}

bands = ["PACS160um", "SPIRE250um", "SPIRE350um", "SPIRE500um"]

low_color_limits = {
    "PACS160um": (-40, -10),
    "SPIRE250um": (0, 30),
    "SPIRE350um": (0, 20),
    "SPIRE500um": (0, 10),
}


#def gaussian(x, *args):
#    mu, sigma, A = args
#    return A*np.exp(-((x - mu)**2)/(sigma**2))


def get_chisq_and_diffManticore(single_temp_fn, band_stub):
    with fits.open(single_temp_fn) as hdul:
        chisq = hdul[5].data * XsMOD
        diffManticore = hdul[6 + bands.index(band_stub)].data
    return chisq, diffManticore


def get_data(fn, head=False):
    with fits.open(fn) as hdul:
        data = hdul[0].data
        if head:
            header = hdul[0].header
            header['CUNIT1'] = "deg"
            header['CUNIT2'] = "deg"
            data = (data, header)
    return data


def gen_hist_and_stats(*args, x_lim=None, log=False, band=None, setting=0):
    if len(args) == 2:
        img, mask = args
    elif len(args) == 1:
        img = args[0]
        mask = np.full(img.shape, True)
    else:
        raise InputError("Too many inputs (%d) to gen_hist_and_stats" % (len(args)))
    if band is None:
        band = "PACS160um"
    if x_lim is None:
        if band == "PACS160um":
            x_lim = diffLIM
        else:
            x_lim = small_diffLIM
    dhist, dedges = np.histogram(img[mask].ravel(), bins=BINS, range=x_lim)
    prep_arr = lambda a, b: np.array([a, b]).T.flatten()
    histx, histy = prep_arr(dedges[:-1], dedges[1:]), prep_arr(dhist, dhist)
    bin_centers = (dedges[:-1]+dedges[1:])/2
    peak_val = np.max(dhist)
    mode = bin_centers[dhist == peak_val][0]
    error = (None, None, None)
    if setting < 3:
        factor = 0.5
    elif setting == 3:
        factor = 0.75
    try:
        spline = UnivariateSpline(bin_centers, dhist - peak_val*factor, s=0)
        r1, r2 = spline.roots()
        # this is interesting to investigate...
        # print(max(np.linspace(r1, r2, 200), key=spline))
        # should probably try gauss-fit to hist>FWHM
        fwhm = np.abs(r1 - r2)
        sigma = fwhm/2.355
        if setting < 0:
            p0 = [mode, sigma, peak_val]
            popt, pcov = curve_fit(gaussian, bin_centers, dhist, p0=p0)
            mode, sigma, A = popt
            error = np.sqrt(np.diag(pcov))
            if setting < -2:
                histy = gaussian(bin_centers, *popt)
                histx = bin_centers
    except ValueError:
        fwhm, sigma = np.nan, np.nan
    if log:
        histy = np.log10(histy)
    if setting <= 0:
        if setting < -1:
            return (histx, histy), (mode, sigma, error)
        else:
            return (histx, histy), (mode, sigma)
    elif setting == 1:
        return mode, sigma
    elif setting > 1:
        return mode, sigma, r2


def get_referenced_mask(reference_image, setting=1):
    mask = ~np.isnan(reference_image)
    assert setting != 0
    if setting <= 3:
        stats_ref = gen_hist_and_stats(reference_image, mask,
            x_lim=fluxLIM, log=False, setting=setting)
        if setting == 1:
            ref_value = stats_ref[0]
        elif setting > 1:
            ref_value = stats_ref[2]
        mask &= (reference_image < ref_value)
    elif setting > 3:
        tolerance = 0.05
        if setting == 4:
            # Within 10%
            tolerance = 0.1
        elif setting == 5:
            # Within 20%
            tolerance = 0.2
        # Setting == 6 is now implicitly within 5%
        mask &= (reference_image < (1 + tolerance)) & (reference_image > (1 - tolerance))
    return mask

def get_referenced_mask_positive(reference_image, setting=1):
    mask = np.isnan(reference_image)
    assert setting != 0
    if setting <= 3:
        stats_ref = gen_hist_and_stats(reference_image, ~mask,
            x_lim=fluxLIM, log=False, setting=setting)
        if setting == 1:
            ref_value = stats_ref[0]
        elif setting > 1:
            ref_value = stats_ref[2]
        mask |= (reference_image < ref_value)
    elif setting > 3:
        tolerance = 0.05
        if setting == 4:
            # Within 10%
            tolerance = 0.1
        elif setting == 5:
            # Within 20%
            tolerance = 0.2
        # Setting == 6 is now implicitly within 5%
        mask |= (reference_image < (1 + tolerance)) & (reference_image > (1 - tolerance))
    return ~mask


def quickrun_flux_offset():
    if len(sys.argv) > 1:
        band_stub = sys.argv[1]
        field_stub = sys.argv[2]
    else:
        band_stub = "PACS160um"
        field_stub = "Per1"

    ref_stub = "SPIRE350um"
    planck_stubs = ['353', '545', '857']

    # Set paths, filename generators
    per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
    gen_herschel_fn = lambda b: "%s%s-image-remapped.fits" % (per1_dir, b)
    gen_planck_fab_fn = lambda b: "%s%s-image-remapped-PLANCKfabricated.fits" % (per1_dir, b)
    planck_ratio_namegen = lambda b: "%sPlanck_ObsToModel_ratio_%s_convolved.fits" % (CURRENT_DIR, b)

    # Generate filenames
    ref_fn = gen_herschel_fn(ref_stub)
    herschel_fn = gen_herschel_fn(band_stub)
    singleT_manticore_fn = "%sT4-absdiff-%sJ-4bandLErr.fits" % (per1_dir, field_stub)
    planck_fn = gen_planck_fab_fn(band_stub)

    # Load stuff in
    ref_img = get_data(ref_fn)
    herschel_img, herschel_head = get_data(herschel_fn, head=True)
    chisq, diffManticore = get_chisq_and_diffManticore(singleT_manticore_fn, band_stub)
    planck_img = get_data(planck_fn)

    # Get the offset mask the old way, with the manticore chi squared filter
    manticore_mask = XsMASK(chisq) & (~np.isnan(diffManticore))
    # Reference SPIRE for 3 masks: mode, 1/2 max, and 3/4 max
    spire_masks = (get_referenced_mask(ref_img, setting=x) for x in range(1, 4))
    # Get Planck ratio mask, within 10% and 20%
    planck_masks = (reduce(lambda x, y: x&y, (get_referenced_mask(get_data(x), setting=s) for x in (planck_ratio_namegen(b) for b in planck_stubs))) for s in (5, 4, 6))

    # Convolve the Herschel map down to Planck resolution
    h_beam, p_beam = bandpass_beam_sizes[band_stub], 5  # arcminutes
    conv_beam = prepare_convolution(WCS(herschel_head), h_beam, p_beam, herschel_img.shape)
    herschel_img = convolve_properly(herschel_img, conv_beam)

    # Subtract Herschel image from Planck image
    diffPlanck = planck_img - herschel_img

    # Get stats and make plots
    plt.figure(figsize=(18, 9))
    ax = plt.subplot(111)

    horiz_offset = 0.7
    # Get stats for manticore + chisq mask and make plot
    hist_manticore, stats_manticore = gen_hist_and_stats(diffManticore, manticore_mask, band=band_stub)
    ax.plot(*hist_manticore, '-', color='r', label='manticore diff, manticore chisq mask', alpha=0.5)
    ax.text(horiz_offset, 0.85, "Mode: %.1f\nSTD: %.2f" % (stats_manticore[0], stats_manticore[1]),
        bbox=dict(fill=False, edgecolor='r', linestyle='-'),
        horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    # Get stats for Planck + chisq mask and make plot
    hist_p_chisq, stats_p_chisq = gen_hist_and_stats(diffPlanck, manticore_mask, band=band_stub)
    ax.plot(*hist_p_chisq, '-', color='b', label='Planck - Herschel, manticore chisq mask', alpha=0.5)
    ax.text(horiz_offset, 0.75, "Mode: %.1f\nSTD: %.2f" % (stats_p_chisq[0], stats_p_chisq[1]),
        bbox=dict(fill=False, edgecolor='b', linestyle='-'),
        horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    # Get stats for SPIRE masks and make plots
    for mask, linestyle, mask_type in zip(spire_masks, ['-', '-.', ':'], ['mode', '1/2 max', '3/4 max']):
        hist_spiremask, stats_spiremask = gen_hist_and_stats(diffPlanck, mask, band=band_stub)
        ax.plot(*hist_spiremask, linestyle, color='teal', label='Planck - Herschel, %s %s mask' % (ref_stub, mask_type), alpha=0.5)
        ax.text(horiz_offset, 0.65, "Mode: %.1f\nSTD: %.2f" % (stats_spiremask[0], stats_spiremask[1]),
            bbox=dict(fill=False, edgecolor='teal', linestyle=linestyle),
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        horiz_offset += 0.08

    horiz_offset = 0.7
    # Get stats for Planck ratio masks and make plots
    for mask, linestyle, mask_type in zip(planck_masks, ['-', '-.', ':'], ['20%', '10%', '5%']):
        hist_planckmask, stats_planckmask = gen_hist_and_stats(diffPlanck, mask, band=band_stub)
        ax.plot(*hist_planckmask, linestyle, color='magenta', label="Planck - Herschel, Planck ratio within %s" % mask_type, alpha=0.5)
        ax.text(horiz_offset, 0.55, "Mode: %.1f\nSTD: %.2f" % (stats_planckmask[0], stats_planckmask[1]),
            bbox=dict(fill=False, edgecolor='magenta', linestyle=linestyle),
            horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        horiz_offset += 0.08


    ax.set_xlabel("Positive difference from original PACS160um flux")
    ax.set_ylabel("Pixel count")
    ax.set_title("%s %s flux difference, manticore and Planck methods, with chi squared, SPIRE flux, and Planck ratio masks"
        % (field_stub, band_stub))
    ax.legend()
    plt.tight_layout()
    # plt.savefig("Figure_8_%s_%s.png" % (field_stub, band_stub))
    show_plot()


def quickrun_masked_region():
    if len(sys.argv) > 1:
        band_stub = sys.argv[1]
        field_stub = sys.argv[2]
    else:
        band_stub = "PACS160um"
        field_stub = "Per1"

    band_stub = sys.argv[1]
    field_stub = "Per1"

    # Set paths, filename generators
    per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
    gen_herschel_fn = lambda b: "%s%s-image-remapped.fits" % (per1_dir, b)
    gen_planck_fab_fn = lambda b: "%s%s-image-remapped-PLANCKfabricated.fits" % (per1_dir, b)

    # Generate filenames
    herschel_fn = gen_herschel_fn(band_stub)
    singleT_manticore_fn = "%sT4-absdiff-%sJ-4bandLErr.fits" % (per1_dir, field_stub)
    planck_fn = gen_planck_fab_fn(band_stub)

    # Load stuff in
    herschel_img, herschel_head = get_data(herschel_fn, head=True)
    chisq, diffManticore = get_chisq_and_diffManticore(singleT_manticore_fn, band_stub)
    planck_img = get_data(planck_fn)

    # Get the chi squared mask
    manticore_mask = XsMASK(chisq) & (~np.isnan(diffManticore))

    # Reference SPIRE350 for a mask
    spire350_mask = get_referenced_mask(get_data(gen_herschel_fn("SPIRE350um")))
    # Get SPIRE500 reference mask as well
    spire500_mask = get_referenced_mask(get_data(gen_herschel_fn("SPIRE500um")))
    # Get SPIRE250 reference mask as well
    spire250_mask = get_referenced_mask(get_data(gen_herschel_fn("SPIRE250um")))

    plt.figure(figsize=(18, 9))
    vmin, vmax = low_color_limits[band_stub]

    plt.subplot(141)
    herschel_img_mm = herschel_img.copy()
    herschel_img_mm[~manticore_mask] = np.nan
    plt.imshow(herschel_img_mm, origin='lower', vmin=vmin, vmax=vmax)
    plt.title("(Herschel %s) Manticore mask" % band_stub)
    plt.colorbar()

    plt.subplot(142)
    herschel_img_2m = herschel_img.copy()
    herschel_img_2m[~spire250_mask] = np.nan
    plt.imshow(herschel_img_2m, origin='lower', vmin=vmin, vmax=vmax)
    plt.title("SPIRE250um<MODE mask")
    plt.colorbar()

    plt.subplot(143)
    herschel_img_3m = herschel_img.copy()
    herschel_img_3m[~spire350_mask] = np.nan
    plt.imshow(herschel_img_3m, origin='lower', vmin=vmin, vmax=vmax)
    plt.title("SPIRE350um<MODE mask")
    plt.colorbar()

    plt.subplot(144)
    herschel_img_5m = herschel_img.copy()
    herschel_img_5m[~spire500_mask] = np.nan
    plt.imshow(herschel_img_5m, origin='lower', vmin=vmin, vmax=vmax)
    plt.title("SPIRE500um<MODE mask")
    plt.colorbar()

    plt.tight_layout()
    # plt.savefig("Figure_9_%s_%s.png" % (field_stub, band_stub))
    show_plot()


def quickrun_investigate_larger_mask():
    np.warnings.filterwarnings('ignore')

    field_stub = "Per1"
    ref_stub = "SPIRE350um"

    # Set paths, filename generators
    per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
    gen_herschel_fn = lambda b: "%s%s-image-remapped.fits" % (per1_dir, b)

    # Generate filenames
    ref_fn = gen_herschel_fn(ref_stub)

    # Load stuff in
    ref_img, ref_head = get_data(ref_fn, head=True)

    img = ref_img
    dhist, dedges = np.histogram(img.ravel(), bins=BINS, range=(0, 40))
    prep_arr = lambda a, b: np.array([a, b]).T.flatten()
    histx, histy = prep_arr(dedges[:-1], dedges[1:]), prep_arr(dhist, dhist)
    bin_centers = (dedges[:-1]+dedges[1:])/2

    mode_ref = bin_centers[dhist==np.max(dhist)][0]
    count_notNaN = img[~np.isnan(img)].size
    n_hi = img[img > mode_ref].size
    n_lo = img[img < mode_ref].size
    n_hi /= count_notNaN
    n_lo /= count_notNaN

    partial_height = np.max(dhist)/2
    hm_stub = "Half Max"
    try:
        spline = UnivariateSpline(bin_centers, dhist - partial_height, s=0)
        r1, r2 = spline.roots()
        fwhm = np.abs(r1 - r2)
        sigma = fwhm/2.355
    except ValueError:
        r1, r2 = np.nan, np.nan
        fwhm, sigma = np.nan, np.nan

    hi_hm_ref = r2
    n2_hi = img[img > hi_hm_ref].size / count_notNaN
    n2_lo = img[img < hi_hm_ref].size / count_notNaN

    plt.figure(figsize=(11, 8.5))

    plt.plot(histx, histy, '-')
    plt.title("{0} {1} (unconvolved)".format(field_stub, ref_stub))

    text_x = 0.7

    plt.plot([mode_ref, mode_ref], [0, np.max(dhist)], '--', color='blue')
    msg = "Mode: {0:.1f} MJy/sr\nAbove/Below: {1:.1%}/{2:.1%}".format(mode_ref, n_hi, n_lo)
    plt.text(text_x, 0.7, msg, transform=plt.gca().transAxes, color='blue')

    plt.plot([hi_hm_ref, hi_hm_ref], [0, partial_height], '--', color='red')
    msg = "{3} value: {0:.1f} MJy/sr\nAbove/Below: {1:.1%}/{2:.1%}".format(hi_hm_ref, n2_hi, n2_lo, hm_stub)
    plt.text(text_x, 0.4, msg, transform=plt.gca().transAxes, color='red')


    partial_height = np.max(dhist)*3/4
    hm_stub = "3/4 Max"
    try:
        spline = UnivariateSpline(bin_centers, dhist - partial_height, s=0)
        r1, r2 = spline.roots()
        fwhm = np.abs(r1 - r2)
        sigma = fwhm/2.355
    except ValueError:
        r1, r2 = np.nan, np.nan
        fwhm, sigma = np.nan, np.nan
    hi_hm_ref = r2
    n2_hi = img[img > hi_hm_ref].size / count_notNaN
    n2_lo = img[img < hi_hm_ref].size / count_notNaN
    plt.plot([hi_hm_ref, hi_hm_ref], [0, partial_height], '--', color='k')
    msg = "{3} value: {0:.1f} MJy/sr\nAbove/Below: {1:.1%}/{2:.1%}".format(hi_hm_ref, n2_hi, n2_lo, hm_stub)
    plt.text(text_x, 0.55, msg, transform=plt.gca().transAxes, color='k')

    plt.tight_layout()
    show_plot()

def quickrun_planck_obsvsmodel_ratio_stats():
    planck_bands = ['353', '545', '857']
    convolved_or_not = ["convolved", "unconvolved"]
    colors = ['b', 'r', 'k']
    styles = [':', '-']
    planck_ratio_namegen = lambda b, u: "Planck_ObsToModel_ratio_%s_%s.fits" % (b, u)
    plt.figure(figsize=(15, 10))
    vert_offset = 0.85
    for i, if_convolved in enumerate(convolved_or_not):
        for j, planck_band in enumerate(planck_bands):
            p_ratio = get_data(planck_ratio_namegen(planck_band, if_convolved))
            histxy, stats = gen_hist_and_stats(p_ratio, x_lim=(0.7, 1.5), band=planck_band)
            plt.plot(*histxy, linestyle=styles[i], alpha=0.7,
                color=colors[j], label="%sGHz %s" % (planck_band, if_convolved))
            if if_convolved == "un":
                plt.text(0.7, vert_offset, "Mode: %.2f\nSTD: %.3f" % (stats[0], stats[1]),
                    bbox=dict(fill=False, edgecolor=colors[j], linestyle=styles[i]),
                    horizontalalignment='center', verticalalignment='center',
                    transform=plt.gca().transAxes)
                vert_offset -= 0.1
    plt.legend()
    plt.xlabel("Ratio: Planck observed / model emission")
    plt.ylabel("Histogram count")
    plt.title("Observed/Model emission ratio histogram for 3 highest frequency HFI bands")
    plt.tight_layout()
    show_plot()

def quickrun_planck_obsvsmodel_ratio_mask():
    planck_bands = ['353', '545', '857']
    convolved_or_not = ["convolved10", "convolved", "unconvolved"]
    colors = ['b', 'r', 'k']
    styles = [':', '--', '-']
    planck_ratio_namegen = lambda b, u: "%sPlanck_ObsToModel_ratio_%s_%s.fits" % (CURRENT_DIR, b, u)
    plt.figure(figsize=(12, 10))
    ratios = []
    for planck_band in planck_bands:
        p_ratio = get_data(planck_ratio_namegen(planck_band, "convolved"))
        ratios.append(p_ratio)

    plt.subplot(221) # Just within 10%
    count_arrays = []
    for p_ratio in ratios:
        count_array = np.zeros(p_ratio.shape)
        count_array[np.isnan(p_ratio)] = np.nan
        count_array[(p_ratio < 1.1) & (p_ratio > 0.9)] += 1
        # count_array[(p_ratio < 1.05) & (p_ratio > 0.95)] += 1
        count_arrays.append(count_array)
    final_count = sum(count_arrays)
    plt.imshow(final_count, origin='lower')
    plt.title("Ratio between [0.9, 1.1]")
    plt.ylabel("the 3 highest frequency HFI bands")
    plt.colorbar()

    plt.subplot(222) # 1pt within 20%, additional pt within 10%
    count_arrays = []
    for p_ratio in ratios:
        count_array = np.zeros(p_ratio.shape)
        count_array[np.isnan(p_ratio)] = np.nan
        count_array[(p_ratio < 1.2) & (p_ratio > 0.8)] += 1
        count_array[(p_ratio < 1.1) & (p_ratio > 0.9)] += 1
        count_arrays.append(count_array)
    final_count = sum(count_arrays)
    plt.imshow(final_count, origin='lower')
    plt.title("[0.8, 1.2], but weighted higher if [0.9, 1.1]")
    plt.colorbar()

    plt.subplot(223) # Just within 5%
    count_arrays = []
    for p_ratio in ratios:
        count_array = np.zeros(p_ratio.shape)
        count_array[np.isnan(p_ratio)] = np.nan
        # count_array[(p_ratio < 1.1) & (p_ratio > 0.9)] += 1
        count_array[(p_ratio < 1.05) & (p_ratio > 0.95)] += 1
        count_arrays.append(count_array)
    final_count = sum(count_arrays)
    plt.imshow(final_count, origin='lower')
    plt.title("Ratio between [0.95, 1.05]")
    plt.ylabel("All plots are integrated over")
    plt.colorbar()

    plt.subplot(224) # 1pt within 10%, additional pt within 5%
    count_arrays = []
    for p_ratio in ratios:
        count_array = np.zeros(p_ratio.shape)
        count_array[np.isnan(p_ratio)] = np.nan
        count_array[(p_ratio < 1.175) & (p_ratio > .7)] += 1
        # count_array[(p_ratio < 1.05) & (p_ratio > 0.95)] += 1
        count_arrays.append(count_array)
    final_count = sum(count_arrays)
    final_count[final_count < 3] = 0
    final_count[final_count >= 3] = 1
    plt.imshow(final_count, origin='lower')
    plt.title("[0.9, 1.1], but weighted higher if [0.95, 1.05]")
    plt.colorbar()

    plt.tight_layout()
    show_plot()


def get_spire_mask(setting=3):
    # Set paths, filename generators
    per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
    gen_herschel_fn = lambda b: "%s%s-image-remapped.fits" % (per1_dir, b)
    arrays = []
    for band_stub in bands[1:]:
        data = get_data(gen_herschel_fn(band_stub))
        mask = get_referenced_mask_positive(data, setting=setting)
        arrays.append(mask)
    mask = np.all(arrays, axis=0)
    return mask


def quickrun_spire_mask():
    # Set paths, filename generators
    per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
    gen_herschel_fn = lambda b: "%s%s-image-remapped.fits" % (per1_dir, b)
    # Set up lists for masks
    masks_mode = []
    masks_12 = []
    masks_34 = []
    combo_mask = []
    all_masks = [masks_mode, masks_12, masks_34, combo_mask]
    # Loop through Herschel SPIRE bands (remapped)
    for band_stub in bands[1:]:
        data = get_data(gen_herschel_fn(band_stub))
        for i in range(1, 4):
            # 1: mode, 2: half-max, 3: 3/4  max
            mask = get_referenced_mask(data, setting=i)
            count_array = np.zeros(data.shape)
            count_array[np.isnan(data)] = np.nan
            count_array[mask] += 1
            all_masks[i-1].append(count_array)
            if band_stub in bands[2:]:
                all_masks[3].append(count_array)
    all_masks = list(map(sum, all_masks))
    titles = [
        "Mode mask counts for SPIRE bands",
        "Half-max counts for SPIRE bands",
        "3/4 max counts for SPIRE bands",
        "All three previous mask options for 350+500 bands"
    ]
    plt.figure(figsize=(12, 10))
    for i, m in enumerate(all_masks):
        plt.subplot(221 + i)
        plt.imshow(m, origin='lower')
        plt.colorbar()
        plt.title(titles[i])
        if i == 2:
            plt.ylabel("Plots integrated over")
        elif i == 0:
            plt.ylabel("multiple SPIRE bands")
    plt.tight_layout()
    show_plot()

def quickrun_spirevsplanck_offset():
    # Comparing 3/4-max SPIRE350um mask to 10% Planck ratio mask (3 highest HFI)
    # Set bands
    herschel_stub = "SPIRE350um"
    planck_stubs = ['353', '545', '857']
    # Set paths, filename generators
    per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
    gen_herschel_fn = lambda b: "%s%s-image-remapped.fits" % (per1_dir, b)
    gen_planck_fab_fn = lambda b: "%s%s-image-remapped-PLANCKfabricated.fits" % (per1_dir, b)
    planck_ratio_namegen = lambda b: "%sPlanck_ObsToModel_ratio_%s_unconvolved.fits" % (CURRENT_DIR, b)
    # Get filenames
    herschel_fn = gen_herschel_fn(herschel_stub)
    planck_fns = [planck_ratio_namegen(b) for b in planck_stubs]
    # Get the 350um 3/4-max mask
    spire_mask = get_referenced_mask(get_data(herschel_fn), setting=3)
    # Get Planck ratio masks, within 10%
    planck_mask = reduce(lambda x, y: x & y, (get_referenced_mask(get_data(x), setting=4) for x in planck_fns))
    # Get PACS data for offset calculation
    herschel_stub = "PACS160um"
    herschel_img, herschel_head = get_data(gen_herschel_fn(herschel_stub), head=True)
    planck_img = get_data(gen_planck_fab_fn(herschel_stub))
    # Convolve the Herschel map down to Planck resolution
    h_beam, p_beam = bandpass_beam_sizes[herschel_stub], 5  # arcminutes
    conv_beam = prepare_convolution(WCS(herschel_head), h_beam, p_beam, herschel_img.shape)
    herschel_img = convolve_properly(herschel_img, conv_beam)
    diff_img = planck_img - herschel_img
    # Get histograms
    hist_spiremask, stats_spiremask = gen_hist_and_stats(diff_img, spire_mask, band=herschel_stub)
    hist_planckmask, stats_planckmask = gen_hist_and_stats(diff_img, planck_mask, band=herschel_stub)
    # Make plots
    plt.figure(figsize=(18, 9))
    ax = plt.subplot(111)
    ax.plot(*hist_spiremask, '-', color='r', label="SPIRE350um mask")
    ax.plot(*hist_planckmask, '-', color='b', label="Planck ratio mask")
    plt.legend()
    show_plot()

if __name__ == "__main__":
    plotting_remotely = True
    import matplotlib
    if plotting_remotely:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    print("started...")
    quickrun_planck_obsvsmodel_ratio_mask()


    # np.warnings.filterwarnings('ignore')
    # if len(sys.argv) > 3:
    #     if sys.argv[3] == "mask":
    #         quickrun_masked_region()
    #     elif sys.argv[3] == "hist":
    #         quickrun_flux_offset()
    # else:
    #     print("nothing here")
