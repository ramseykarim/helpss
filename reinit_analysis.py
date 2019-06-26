if __name__ == "__main__":
    plotting_remotely = True
    import matplotlib
    if plotting_remotely:
        matplotlib.use('Agg')
import matplotlib.pyplot as plt
SAVE_NAME = "/home/rkarim/Downloads/Figure_X_current.png"
import numpy as np
from matplotlib.patches import Patch
from astropy.io import fits
from planck_mask import gen_hist_and_stats
import manticore_results as mtc
from planck_mask import get_spire_mask
from boolean_islands import get_mask, get_planck_mask, fill_inwards
import sys


def show_plot():
    if plotting_remotely:
        plt.savefig(SAVE_NAME)
    else:
        plt.show()


per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
soln_DL = "T4-absdiff-Per1J-3param-plus045-cDL3hpow-1000-0.1-1.80-bcreinit-Nh5E19,2E22.fits"
soln_OH = "T4-absdiff-Per1J-3param-plus045-cOH5hpow-1000-0.1-1.80-bcreinit-Nh5E19,2E22.fits"
soln_15 = "T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-1.50hpow-1000-0.1-1.80-bcreinit-Nh1E20.fits"
soln_0 = "T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-bcreinit-Nh0.0.fits"
soln_ideal = "T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-bcreinit-Nh1.4E20,5E21.fits"
soln_nominal = "T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-bcreinit-Nh1E20.fits"
soln_lowT = "T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-bcreinit-Th15.5-Nh3E20,2E21.fits"

soln_c19h16_toodense = "T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-1.90hpow-1000-0.1-1.60-bcreinit-Th16.0-Nh5e20,3e22.fits"
soln_c20h16 = "T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.00hpow-1000-0.1-1.60-bcreinit-Th16.0-Nh5e20,3e22.fits"
soln_lowerT16 = "T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.10hpow-1000-0.1-1.60-bcreinit-Th15.0-Nh5e20,3e22.fits"
soln_c21h16_tighter = "T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.10hpow-1000-0.1-1.60-bcreinit-Th16.0-Nh1e21,1e22.fits"
soln_c21h16_toodense = "T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.10hpow-1000-0.1-1.60-bcreinit-Th16.0-Nh5e20,3e22.fits"

soln_OH_16 = "T4-absdiff-Per1J-3param-plus045-cOH5hpow-1000-0.1-1.60-bcreinit-Th16.0-Nh5E19,1E22.fits"
soln_c19h16 = "T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-1.90hpow-1000-0.1-1.60-bcreinit-Th16.0-Nh5E19,1E22.fits"
soln_c21h16 = "T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.10hpow-1000-0.1-1.60-bcreinit-Th16.0-Nh5E19,1E22.fits"
soln_c21h17 = "T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.10hpow-1000-0.1-1.70-bcreinit-Th15.5-Nh5E19,1E22.fits"

soln_5pcterr = "T4-absdiff-Per1J-3param-plus045-plus05.0pct-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-bcreinit-Th15.95-Nh5E19,2E22.fits"
soln_2p_5pcterr = "T4-absdiff-Per1J-plus045-plus05.0pct-pow-1000-0.1-1.80.fits"

filament_mask_fn = per1_dir+"filament_mask_syp.fits"
# mask = fits.getdata(filament_mask_fn).astype(bool)

def get_Xs(fn):
    return mtc.load_specific_frame(fn, 9)

def plot_compare_Xs(img1, img2, label1, label2):
    img1_better = (img1 < img2).astype(float)
    img2_better = (img1 > img2).astype(float)
    for x in (img1_better, img2_better):
        x[~mask] = np.nan
    plt.figure(figsize=(14, 9))
    plt.subplot(121)
    plt.imshow(img1_better, origin='lower')
    plt.title("{} is best here".format(label1))
    plt.subplot(122)
    plt.imshow(img2_better, origin='lower')
    plt.title("{} is best here".format(label2))
    show_plot()


def try_Ngt3e21_mask():
    img_2p = mtc.load_specific_frame(soln_2p_5pcterr, 3)
    plt.figure()
    plt.imshow((img_2p > 2e21).astype(int), origin='lower')
    show_plot()
    return
    return (img_2p > 1e21)
    # the rest of this is plotting / comparison
    img_3p = mtc.load_specific_frame(soln_5pcterr, 3)
    nanmask = np.isnan(img_3p)
    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(14, 9))
    for img, ax, l, cutoff in zip((img_2p, img_3p), axes, ('2', '3'), (1e21, 3e21)):
        plt.sca(ax)
        img[nanmask] = np.nan
        plt.imshow((img<cutoff).astype(int), origin='lower')
        plt.colorbar()
        plt.title(l)
    show_plot()

try_Ngt3e21_mask()
