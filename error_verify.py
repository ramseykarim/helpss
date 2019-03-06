import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

#dir_stub = "./Per/testregion1342190326/"
dir_stub = "./Ser/testregion1342206676/"
fn_stub = dir_stub+"T3-reldiff-Ser2-"
threeband_fn = fn_stub+"3band.fits"
#fourband_fn = fn_stub+"4bandLErr.fits"
fullfit_fn = dir_stub+"TEST3-reldiff-region2.fits"


with fits.open(threeband_fn) as hdul:
    tbT = hdul[1].data
#with fits.open(fourband_fn) as hdul:
#    fbT = hdul[1].data
"""
diff1 = tbT - fbT
diff = diff1[~np.isnan(diff1)]
print("MEAN", np.mean(diff))
print("MEDIAN", np.median(diff))
print("STD", np.std(diff))
print("min, max", np.min(diff), np.max(diff))

hdr1 = fits.Header()
hdr1['COMMENT'] = "SPIRE_Only - 4band_large_PACS_error T Map"
hdu1 = fits.PrimaryHDU(diff1, header=hdr1)
hdu1.writeto(dir_stub+"diff_T_3bandMethods.fits")
"""
with fits.open(fullfit_fn) as hdul:
    fullT = hdul[1].data

diff2 = tbT - fullT
hdr2 = fits.Header()
hdr2['COMMENT'] = "SPIRE_Only - 4band_full_fit T Map"
hdu2 = fits.PrimaryHDU(diff2, header=hdr2)
hdu2.writeto(dir_stub+"diff_T_34band.fits")

hdr3 = fits.Header()
hdr3['COMMENT'] = "4band_full_fit T Map"
hdu3 = fits.PrimaryHDU(fullT, header=hdr3)
hdu3.writeto(dir_stub+"Tmap_4bandFull.fits")

