plotting_remotely = True
import numpy as np
import matplotlib
if plotting_remotely:
	matplotlib.use('Agg')
SAVE_NAME = "Figure_X_current.png"
import matplotlib.pyplot as plt
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


Xs_frame = 5

if __name__ == "__main__":
    Xs_img = fits.getdata(cold_single_fn, Xs_frame)
    plt.imshow(Xs_img, origin='lower', vmin=0.99, vmax=1.01)
    show_plot()



def write_cold_soln():
    raise RuntimeError("Should not rerun this function")
    bc_mask = fits.getdata(cold_mask_fn).astype(bool)
    msg1 = "Ramsey forced this to be a mask (April 20, 2019)"
    msg2 = "original: T4-absdiff-Per1J-plus045-pow-1000-0.1-1.80.fits"
    msg3 = "original found in {}".format(single_comp_dir)
    with fits.open(single_fn) as hdul:
        Xs_img = hdul[Xs_frame].data
        hdul[Xs_frame].data[bc_mask & ~np.isnan(Xs_img)] = 5.0
        hdul[Xs_frame].data[~bc_mask & ~np.isnan(Xs_img)] = 0.0
        hdul[Xs_frame].header['COMMENT'] = msg1
        hdul[Xs_frame].header['COMMENT'] = msg2
        hdul[Xs_frame].header['COMMENT'] = msg3
        hdul.writeto(cold_single_fn)
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
