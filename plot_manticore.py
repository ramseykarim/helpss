import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import sys
import boolean_islands as boolis
from datetime import datetime, timezone

per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
power_stub = "pow-1000-0.1-"
power_run_stub = "T4-absdiff-Per1J-plus045-"
power_run_3p_stub = "T4-absdiff-Per1J-3param-plus045-"
single_comp_dir = "single_comp_beta_grid/"
two_comp_dir_stub = "two_comp_beta_grid_betah"
fits_stub = ".fits"

GNILC_T_fn = per1_dir + "dustModel_Per1_SPIRE500umgrid_Temperature.fits"
multi_beta_fn = per1_dir + "single_component_multi_beta.fits"

def gen_power_fn(beta):
	return "{:s}{:s}{:s}{:s}{:s}{:s}".format(
		per1_dir, single_comp_dir, power_run_stub,
		power_stub, beta, fits_stub
	)

def gen_3p_power_fn(beta_c, beta_h):
	return "{:s}{:s}{:s}/{:s}c{:s}{:s}h{:s}{:s}{:s}".format(
		per1_dir, two_comp_dir_stub, beta_h,
		power_run_3p_stub,
		power_stub, beta_c,
		power_stub, beta_h,
		fits_stub
	)

local_test = "../T4-absdiff-Per1J-3param-plus046-full.fits"
b21n1e20 = "../T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-bcreinit-Nh1E20.fits"
b21n0 = "../T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-bcreinit-Nh0.0.fits"
b15n1e20 = "../T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-1.50hpow-1000-0.1-1.80-bcreinit-Nh1E20.fits"

mask = fits.getdata("../filament_mask_syp.fits").astype(bool)
data, hdr = fits.getdata(b21n1e20, 1, header=True)
data[~mask] = np.nan
w = WCS(hdr)

plt.subplot(projection=w)
plt.imshow(data, origin='lower', vmin=8, vmax=13)
plt.xlabel("Right Ascension")
plt.ylabel("Declination")
plt.xlim([50, 850])
plt.ylim([150, 850])
plt.show()


### GENERATING AND WRITING MASK TO FITS
def write_mask():
    data = fits.getdata(local_test, 1)
    nanmask = ~np.isnan(data)
    data, hdr = fits.getdata(b15n1e20, 1, header=True)
    mask = (data>11)
    mask = boolis.get_mask(mask, n=2, min_size=100, dilation=0)
    mask = boolis.fill_inwards(mask, nanmask)
    w = WCS(hdr)

    h['COMMENT'] = "BEST FILAMENT MASK TO DATE (April 23 2019)"
    h = fits.Header()
    h['CREATOR'] = ("Ramsey: {}".format(str(__file__)), "FITS file creator")
    h['OBJECT'] = ("per_04_nom", "Target name")
    h['DATE'] = (datetime.now(timezone.utc).astimezone().isoformat(),
    	"File creation date")
    h.update(w.to_header())
    fits.writeto("../filament_mask_syp.fits", mask.astype(int), header=h)
    print("WRITTEN")
