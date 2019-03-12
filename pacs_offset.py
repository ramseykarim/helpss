plotting_remotely = True
import numpy as np
from astropy.io.fits import getdata
import matplotlib
if plotting_remotely:
	matplotlib.use('Agg')
SAVE_NAME = "Figure_X_current.png"
import matplotlib.pyplot as plt
from planck_mask import gen_hist_and_stats
import manticore_results as mtc
from compare_images import prepare_convolution, convolve_properly, plot_compare, gaussian
from AA_planck_obs  import DustModel

def show_plot():
	if plotting_remotely:
		plt.savefig(SAVE_NAME)
	else:
		plt.show()

img_stub = "-image-remapped.fits"
band_stub = "PACS160um"
region_directory = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
reference_filename = region_directory + band_stub + img_stub
h_data, h_header = getdata(reference_filename, header=True)
dust_directory = "/n/sgraraid/filaments/data/TEST4/regridding_stuff/"
sky = DustModel("Per1", band_stub, h_data, h_header, dust_directory)
planck_flux = sky.observe_planck(band_stub)

# Error
# Herschel_bands/PACS160um_fromManticore.dat not found.

plt.figure(figsize=(12, 10))
plt.imshow(planck_flux, origin='lower')
show_plot()
