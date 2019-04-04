import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import AA_regrid_healpy as rhp

gb_dir = "/home/rkarim/Research/Filaments/Data/"
gb_fn_gen = lambda s: "{}perseus04-{}.fits".format(gb_dir, str(int(s)))
gb_T_fn = gb_dir + "HGBS_perseus_dust_temperature_map.fits"
band_shortcut = {
	160: "PACS160um",
	250: "SPIRE250um",
	350: "SPIRE350um",
	500: "SPIRE500um"
}
per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
our_fn_gen = lambda s: "{}{}-image.fits".format(per1_dir, band_shortcut[s])
our_T_fn = per1_dir + "T4-absdiff-Per1J-plus045-DL3.fits"

b = 160

# load in fits files and headers
# gb_fits, gb_head = fits.getdata(gb_fn_gen(b), header=True)
# our_fits, our_head = fits.getdata(our_fn_gen(b), header=True)
# T map
gb_fits, gb_head = fits.getdata(gb_T_fn, header=True)
our_fits, our_head = fits.getdata(our_T_fn, header=True)


# prepare box so I don't interp a massive uninteresting area
ra_center, dec_center = 52.29622, 30.923503 # degrees
ra_width, dec_width = 3.3, 3 # degrees
ra_lim = ra_center - ra_width/2, ra_center + ra_width/2
dec_lim = dec_center - dec_width/2, dec_center + dec_width/2

# use the general-use regrid function I just wrote
gb_img = rhp.regrid_to_reference(our_fits, our_head, gb_fits, gb_head)
plt.subplot(111)
plt.imshow(gb_img - our_fits, origin='lower', vmin=-1, vmax=4)
plt.title("GB minus our map, {}".format("temperature"))
plt.colorbar()
plt.show()

