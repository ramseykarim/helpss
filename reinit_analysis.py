if __name__ == "__main__":
    plotting_remotely = False
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


per1_dir = "../"
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

soln_16grid = "T4-absdiff-Per1J-3param-plus045-plus05.0pct-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-1.2.4-Th15.95-Nh5E19,2E22.fits"

# ALL CROP 6
manticore_nominal_2p = "T4-absdiff-Per1J-plus045-plus05.0pct-pow-1000-0.1-1.80-crop6.fits"
manticore_nominal_3p = "T4-absdiff-Per1J-3param-plus045-plus05.0pct-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-bcreinit-Th15.95-Nh5E19,2E22-crop6.fits"
mpy_nominal_2p = "mantipyfit_save_test.fits" # T N
mpy_nominal_3p = "mantipyfit_save_test_3p.fits" # Tc Nh Nc
mpy_boundary_soln = "mantipyfit_boundary_solutions-crop6.fits" # Tc2 Nc2 Xs2 Nh1 Xs1 Tc3 Nh3 Nc3 Xs3



def get_Xs(fn):
    return fits.getdata(per1_dir+fn, 9)

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


def try_Ngt3e21_mask(frame):
    n_2p = fits.getdata(per1_dir+soln_2p_5pcterr, 3)
    n_3p = fits.getdata(per1_dir+soln_5pcterr, 3)
    T_3p = fits.getdata(per1_dir+soln_5pcterr, 1)
    # img_3p = fits.getdata(per1_dir+soln_5pcterr, frame)
    # mask = (n_3p > 3e21)
    # T_3p[~mask] = np.nan
    # n_3p[~mask] = np.nan
    plt.figure()
    # plt.subplot(121)
    plt.imshow(n_2p, origin='lower', vmin=2e21, vmax=8e21)
    plt.colorbar()
    # plt.title("Tc")
    # # plt.subplot(122)
    # plt.imshow(n_3p, origin='lower', vmin=3e21, vmax=1e22)
    # plt.colorbar()
    # plt.title("Nc")
    show_plot()
    return
    # return (img_2p > 1e21)
    # # the rest of this is plotting / comparison
    # img_3p = fits.getdata(per1_dir+soln_5pcterr, 3)
    # nanmask = np.isnan(img_3p)
    # fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(14, 9))
    # for img, ax, l, cutoff in zip((img_2p, img_3p), axes, ('2', '3'), (1e21, 3e21)):
    #     plt.sca(ax)
    #     img[nanmask] = np.nan
    #     plt.imshow((img<cutoff).astype(int), origin='lower')
    #     plt.colorbar()
    #     plt.title(l)
    # show_plot()


"""
COMPARING MANTICORE TO MANTIPYTHON (NO TRICKS, JUST NOMINAL)
"""

def compare_2p_N(): # very identical
    # manticore 2p vs mantipy nominal 2p
    mtN = np.log10(fits.getdata(per1_dir+manticore_nominal_2p, 3))
    mpN = np.log10(fits.getdata(per1_dir+mpy_nominal_2p, 2))
    plt.imshow(mtN-mpN, origin='lower')
    plt.title("MANTICORE minus PYTHON: 2param N")
    plt.colorbar()
    plt.show()

def compare_3p_Nc(): # significantly different in off-filament (close on filament)
    # manticore 3p vs mantipy nominal 3p
    mtN = (fits.getdata(per1_dir+manticore_nominal_3p, 3))
    mpN = (fits.getdata(per1_dir+mpy_nominal_3p, 3))
    plt.imshow(np.abs(mtN-mpN), origin='lower', vmin=0, vmax=5e21)
    plt.title("MANTICORE minus PYTHON: 3param Nc")
    plt.colorbar()
    plt.show()

def compare_3p_Nh(): # different in off-filament (close on filament)
    # manticore 3p vs mantipy nominal 3p
    mtN = (fits.getdata(per1_dir+manticore_nominal_3p, 7))
    mpN = (fits.getdata(per1_dir+mpy_nominal_3p, 2))
    plt.imshow(np.abs(mtN-mpN), origin='lower', vmin=0, vmax=5e20)
    plt.title("MANTICORE minus PYTHON: 3param Nh")
    plt.colorbar()
    plt.show()

"""
COMPARING MANTIPYTHON with tricks
"""

def compare_3p_Nc_7lim(): # these are identical
    # compare 3p Nc between T>0 and T>7
    mpN = fits.getdata(per1_dir+mpy_nominal_3p, 3)
    mp7N = fits.getdata(per1_dir+mpy_boundary_soln, 8)
    plt.imshow(np.abs(mpN-mp7N), origin='lower', vmin=0, vmax=5e21)
    plt.title("T>0 minus T>7: 3param Nc")
    plt.colorbar()
    plt.show()

def compare_3p_Nh_7lim(): # these are identical
    # compare 3p Nh between T>0 and T>7
    mpN = fits.getdata(per1_dir+mpy_nominal_3p, 2)
    mp7N = fits.getdata(per1_dir+mpy_boundary_soln, 7)
    plt.imshow(np.abs(mpN-mp7N), origin='lower', vmin=0, vmax=5e20)
    plt.title("T>0 minus T>7: 3param Nh")
    plt.colorbar()
    plt.show()

"""
COMPARING MANTICORE WITH BOUNDARY SOLUTIONS
"""

def compare_Nh_3p_1p():
    mtN = fits.getdata(per1_dir+manticore_nominal_3p, 7)*2 # two layers!
    mpN = fits.getdata(per1_dir+mpy_boundary_soln, 4)
    plt.imshow((mpN-mtN)/mpN, origin='lower')
    plt.title("1p minus 3p: Nh at T=15.95")
    plt.colorbar()
    plt.show()


def mask_Tc_on_Nhboundary():
    mtN = fits.getdata(per1_dir+manticore_nominal_3p, 7)*2 # two layers!
    mpN = fits.getdata(per1_dir+mpy_boundary_soln, 4)
    mask = ((mpN-mtN)/mpN > 0.1)
    mtT = fits.getdata(per1_dir+manticore_nominal_3p, 1)
    mtT[~mask] = np.nan
    plt.imshow(mtT, origin='lower')
    plt.title("manticore Tc 3p")
    plt.colorbar()
    plt.show()


"""
UNCERTAINTY MASK
based on idea that, as 3p Nh approaches 1p Nh, cold component uncertainties should blow out
specificially as compared to their values
"""
def mask_err_sidebysideplot():
    with fits.open(per1_dir+soln_5pcterr) as hdul:
        Te = hdul[2].data / hdul[1].data
        Ne = hdul[4].data / hdul[3].data
    # with fits.open(per1_dir + soln_5pcterr) as hdul:
    #     T = hdul[1].data
    #     N = hdul[3].data
    # with fits.open(per1_dir + soln_16grid) as hdul:
    #     Te = hdul[2].data / T
    #     Ne = hdul[4].data / N
    fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
    axes[0, 0].imshow(Te, origin='lower', vmin=0, vmax=1)
    axes[0, 0].set_title("Terr")
    axes[0, 1].imshow(Ne, origin='lower', vmin=0, vmax=1)
    axes[0, 1].set_title("Nerr")
    axes[1, 0].imshow(fits.getdata(per1_dir+soln_5pcterr, 1), origin='lower', vmin=5, vmax=16)
    axes[1, 0].set_title("Tc")    
    plt.show()

def mask_Nherr_plot():
    with fits.open(per1_dir+manticore_nominal_3p) as hdul:
        Te = hdul[2].data / hdul[1].data
        Ne = hdul[4].data / hdul[3].data
        Nhe = hdul[8].data / hdul[7].data
    T = fits.getdata(per1_dir+manticore_nominal_3p, 1)
    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True)
    p = axes[0].imshow(T, origin='lower', vmin=5, vmax=16)
    fig.colorbar(p, ax=axes[0])
    p = axes[1].imshow(Nhe, origin='lower', vmin=0, vmax=1)
    fig.colorbar(p, ax=axes[1])
    plt.show()

def mask_err():
    with fits.open(per1_dir+manticore_nominal_3p) as hdul:
        Te = hdul[2].data / hdul[1].data
        Ne = hdul[4].data / hdul[3].data
        Nhe = hdul[8].data / hdul[7].data
    T = fits.getdata(per1_dir+manticore_nominal_3p, 1)
    mask = (Te > .2) | (Ne > .4)
    T[mask] = np.nan
    fig, axes = plt.subplots(nrows=1, ncols=2)
    axes[0].imshow(T, origin='lower', vmin=5, vmax=16)
    axes[1].imshow(fits.getdata(per1_dir+manticore_nominal_3p, 1), origin='lower', vmin=5, vmax=16)
    plt.show()

def mask_ratio_Ns():
    with fits.open(per1_dir+manticore_nominal_3p) as hdul:
        Tc = hdul[1].data
        Nc = hdul[3].data
        Nh = hdul[7].data
    with fits.open(per1_dir+manticore_nominal_2p) as hdul:
        T1 = hdul[1].data
        N1 = hdul[3].data
    mask = (N1 < 3e21) | (Nc < 3e21)
    Tc[mask] = np.nan
    fig, axes = plt.subplots(nrows=2, ncols=3, sharex=True, sharey=True)
    axes[0, 0].imshow(Tc, origin='lower', vmin=5, vmax=16)
    axes[0, 0].set_title("Tc masked")
    axes[0, 1].imshow(fits.getdata(per1_dir+manticore_nominal_3p, 1), origin='lower', vmin=5, vmax=16)
    axes[0, 1].set_title("Tc")
    axes[1, 0].imshow(Nc, origin='lower', vmin=0, vmax=1e22)
    axes[1, 0].set_title("Nc")
    axes[1, 1].imshow(N1, origin='lower', vmin=0, vmax=1e22)
    axes[1, 1].set_title("N")
    axes[1, 2].imshow(T1, origin='lower', vmin=14, vmax=17)
    axes[1, 2].set_title("T")
    axes[0, 2].imshow(Nh, origin='lower', vmin=0, vmax=3e21)
    axes[0, 2].set_title("Nh")
    plt.show()


mask_ratio_Ns()