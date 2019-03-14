plotting_remotely = True
import numpy as np
import matplotlib
if plotting_remotely:
	matplotlib.use('Agg')
SAVE_NAME = "Figure_X_current.png"
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from astropy.io import fits
from planck_mask import gen_hist_and_stats
import manticore_results as mtc
from planck_mask import get_spire_mask
from boolean_islands import get_mask, get_planck_mask, fill_inwards
from reanalyze_manticore import histogram, gaussian
import sys

def show_plot():
	if plotting_remotely:
		plt.savefig(SAVE_NAME)
	else:
		plt.show()

dir_per1 = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"

herschel_flux = ("PACS160um-plus046", "SPIRE250um", "SPIRE350um", "SPIRE500um")
img_stub = "-image"
rcc_stub = "-remapped-conv.fits"

gen_flux_fn = lambda i: dir_per1 + herschel_flux[i] + img_stub + rcc_stub

soln_2p_plus000 = "T4-absdiff-Per1J.fits"
soln_2p_plus043 = "T4-absdiff-Per1J-plus043.fits"
soln_2p_plus045 = "T4-absdiff-Per1J-plus045.fits"
soln_2p_plus046 = "T4-absdiff-Per1J-plus046.fits"
soln_2p_plus047 = "T4-absdiff-Per1J-plus047.fits"

soln_3p_plus043 = "T4-absdiff-Per1J-3param.fits"
soln_3p_plus045 = "T4-absdiff-Per1J-3param-plus045.fits"
soln_3p_plus045_15p5 = "T4-absdiff-Per1J-3param-plus045-15p5.fits"
soln_3p_plus046 = "T4-absdiff-Per1J-3param-plus046.fits"
soln_3p_plus046f = "T4-absdiff-Per1J-3param-plus046-full.fits"
soln_3p_plus046_15p6 = "T4-absdiff-Per1J-3param-plus046-15p6.fits"
soln_3p_plus046_15p8 = "T4-absdiff-Per1J-3param-plus046-15p8.fits"
soln_3p_plus046f_15p8 = "T4-absdiff-Per1J-3param-plus046-full-15p8.fits"
soln_3p_plus047 = "T4-absdiff-Per1J-3param-plus047.fits"
soln_3p_plus047_15p5 = "T4-absdiff-Per1J-3param-plus047-15p5.fits"

nominal_2p_soln = soln_2p_plus045
nominal_3p_soln = soln_3p_plus045


def load_manticore(filename, frames=None):
	if frames is None:
		frames = (1, 3, 5)
	with fits.open(dir_per1 + filename) as hdul:
		imgs = tuple(hdul[x].data for x in frames)
	return imgs


T, Xs = load_manticore(nominal_2p_soln, frames=(1, 5))
hist_limits = (10, 20)
img = T
mask = mtc.get_pacs_mask() & (Xs < 1)

dhist, dedges = np.histogram(img[mask].ravel(), bins=BINS, range=x_lim)
prep_arr = lambda a, b: np.array([a, b]).T.flatten()
histx, histy = prep_arr(dedges[:-1], dedges[1:]), prep_arr(dhist, dhist)
bin_centers = (dedges[:-1]+dedges[1:])/2
# try to fit the top half only

try:
	spline = UnivariateSpline(bin_centers, dhist - peak_val*factor, s=0)
	r1, r2 = spline.roots()
	fwhm = np.abs(r1 - r2)
	sigma = fwhm/2.355
except ValueError:
	fwhm, sigma = np.nan, np.nan
peak_val = np.max(dhist)
mode = bin_centers[dhist == peak_val][0]
p0 = [mode, sigma, peak_val]
popt, pcov = curve_fit(gaussian, bin_centers, dhist, p0=p0)
mode, sigma, A = popt

plt.figure(figsize=(11, 8.5))
plt.plot(histx, histy, '-')
plt.plot(bin_centers, gaussian(bin_centers, *popt), '--')
plt.tight_layout()
show_plot()
