plotting_remotely = False
SAVE_NAME = "Figure_X_current.png"
if __name__ == "__main__":
    import matplotlib
    if plotting_remotely:
    	matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reanalyze_manticore import cleaner_core_mask, cleaner_junk_mask, masking_attempt
from boolean_islands import get_mask, fill_inwards
import manticore_results as mtc
from datetime import datetime, timezone

def show_plot():
	if plotting_remotely:
		plt.savefig(SAVE_NAME)
	else:
		plt.show()

single_comp_dir = "single_comp_beta_grid/"
per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
single_fn = "T4-absdiff-Per1J-plus045-pow-1000-0.1-1.80.fits" # 15.95
single_fn = per1_dir + single_comp_dir + single_fn
cold_mask_fn = per1_dir + "bc_cold_mask.fits"
cold_single_fn = per1_dir + "bc_single_plus045_pow1.80_soln.fits"
# updated single component run
single_fn_5pct = "./T4-absdiff-Per1J-plus045-plus05.0pct-pow-1000-0.1-1.80.fits"
single_fn_5pct = per1_dir + single_fn_5pct
# masked on 2p N > 1e21
N_mask_fn = per1_dir + "bc_2pN1e21_mask.fits"
N_single_fn = per1_dir + "bc_single_Nmask_05.0pct_plus045_pow1.80_soln.fits"


Xs_frame = 5
if __name__ == "__main__":
    Xs_img, h = fits.getdata(single_fn_5pct, Xs_frame, header=True)
    # X0, dX = .5, 0.01
    # plt.imshow(Xs_img, origin='lower', vmin=X0-dX, vmax=X0+dX)
    # show_plot()



def write_soln():
    raise RuntimeError("Should not rerun this function")
    bc_mask = fits.getdata(N_mask_fn).astype(bool)
    # Previous version: April 20, 2019
    msg1 = "Ramsey forced this to be a mask (June 28, 2019)"
    msg2 = "original: {}".format(single_fn_5pct)
    msg3 = "original found in {}".format(per1_dir)
    with fits.open(single_fn_5pct) as hdul:
        Xs_img = hdul[Xs_frame].data
        hdul[Xs_frame].data[bc_mask & ~np.isnan(Xs_img)] = 5.0
        hdul[Xs_frame].data[~bc_mask & ~np.isnan(Xs_img)] = 0.0
        hdul[Xs_frame].header['COMMENT'] = msg1
        hdul[Xs_frame].header['COMMENT'] = msg2
        hdul[Xs_frame].header['COMMENT'] = msg3
        hdul.writeto(N_single_fn)
    print("written")


def write_N_mask():
    N2p = fits.getdata(single_fn_5pct, 3)
    N2pmask = (N2p>1e21)
    plotting = False
    if plotting:
        plt.subplot(121)
        plt.imshow(N2pmask.astype(int), origin='lower')
    N2pmask = get_mask(N2pmask, n=1, min_size=60, dilation=1)
    N2pmask = fill_inwards(N2pmask, (~np.isnan(Xs_img) & mtc.get_pacs_mask()), min_size=100)
    N2pmask = get_mask(N2pmask, n=1, min_size=60, dilation=0)
    if plotting:
        plt.subplot(122)
        plt.imshow(N2pmask.astype(int), origin='lower')
        plt.show()
    phdu = fits.PrimaryHDU()
    header = fits.Header()
    header.update(WCS(h).to_header())
    header['COMMENT'] = "Cold gas mask for speeding up manticore-bc"
    header['CREATOR'] = ("Ramsey: {}".format(str(__file__)), "FITS file creator")
    header['OBJECT'] = ("per_04_nom", "Target name")
    header['DATE'] = (datetime.now(timezone.utc).astimezone().isoformat(),
        "File creation date")
    ihdu = fits.ImageHDU(N2pmask.astype(int), header=header)
    hdul = fits.HDUList([phdu, ihdu])
    hdul.writeto(N_mask_fn)
    print("written")


def write_cold_mask():
    # Dilate the "core mask" from earlier construction
    cold_mask = mtc.get_core_mask() #| mtc.get_cold_mask()
    cold_mask = get_mask(cold_mask, n=1, min_size=400, dilation=5)
    # Fill in gaps and clean it back up
    cold_mask = fill_inwards(cold_mask, (~np.isnan(Xs_img) & mtc.get_pacs_mask()))
    cold_mask = get_mask(cold_mask, n=1, min_size=400, dilation=0)
    # Add in a square to capture some thinner emission near the filaments
    ii, jj = np.meshgrid(*(np.arange(x) for x in cold_mask.shape), indexing='ij')
    rectangle_mask = (ii > 300) & (ii < 500) & (jj > 350) & (jj < 550)
    cold_mask |= rectangle_mask
    raise RuntimeError("Should not rerun this function")
    # Write the mask to a new FITS file
    phdu = fits.PrimaryHDU()
    header = fits.Header()
    header.update(WCS(h).to_header())
    header['COMMENT'] = "Cold gas mask for speeding up manticore-bc"
    header['CREATOR'] = ("Ramsey: {}".format(str(__file__)), "FITS file creator")
    header['OBJECT'] = ("per_04_nom", "Target name")
    header['DATE'] = (datetime.now(timezone.utc).astimezone().isoformat(),
        "File creation date")
    ihdu = fits.ImageHDU(cold_mask.astype(int), header=header)
    hdul = fits.HDUList([phdu, ihdu])
    hdul.writeto(cold_mask_fn)
    print("written")


if __name__ == "__main__":
    write_soln()
    # Xs_img = fits.getdata(cold_single_fn, Xs_frame)
    # plt.imshow(Xs_img, origin='lower', vmin=0.99, vmax=1.01)
    # show_plot()

