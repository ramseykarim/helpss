import numpy as np
from astropy.io import fits
import utils_regrid_healpy as rgu
import matplotlib.pyplot as plt

directory = "/home/ramsey/Documents/Research/Filaments/"
hp_fn = directory + "HFI_SkyMap_857-field-Int_2048_R3.00_full.fits"
fits_fn = directory + "HFI_SkyMap_857_resSPIRE350.fits"

m = rgu.open_healpix(hp_fn, nest=False)
data, head = fits.getdata(fits_fn, header=True)
result = rgu.healpix2fits(m, data, head)
plt.imshow(result, origin='lower')
plt.show()
