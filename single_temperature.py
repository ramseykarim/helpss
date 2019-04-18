plotting_remotely = False
import numpy as np
import matplotlib
if plotting_remotely:
	matplotlib.use('Agg')
SAVE_NAME = "Figure_X_current.png"
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from astropy.io import fits
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
from planck_mask import gen_hist_and_stats, get_spire_mask
import manticore_results as mtc
from boolean_islands import get_mask, get_planck_mask, fill_inwards
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

nominal_2p_soln = soln_2p_plus045 # 15.7 K
nominal_3p_soln = soln_3p_plus045

dl3_2p_soln = "T4-absdiff-Per1J-plus045-DL3.fits" # this one is 14.2 K!
soln_2p_pow16 = "T4-absdiff-Per1J-plus045-pow-1000-0.1-1.60.fits" # 17.40!
soln_2p_pow17 = "T4-absdiff-Per1J-plus045-pow-1000-0.1-1.70.fits" # 16.64
soln_2p_pow18 = "T4-absdiff-Per1J-plus045-pow-1000-0.1-1.80.fits" # 15.95

def load_manticore(filename, frames=None):
	if frames is None:
		frames = (1, 3, 5)
	with fits.open(dir_per1 + filename) as hdul:
		imgs = tuple(hdul[x].data for x in frames)
	return imgs

def gaussian(x, mu, sigma, A):
	return A*np.exp(-1*(x - mu)**2 / (2*sigma*sigma))

def poly(x, a, b, c):
	return a*x*x + b*x + c

prep_arr = lambda a, b: np.array([a, b]).T.flatten()

T, Xs = load_manticore(soln_2p_pow18, frames=(1, 5))


BINS = 128
hist_limits = (10., 20.)

def gen_chisq_mask(Xs_limit):
	return (Xs < Xs_limit) & mtc.get_pacs_mask()

def center_of_mass(x, y):
	# Get center of mass of whatever is above half-max
	half_max_val = np.max(y)/2
	yh = y[y>half_max_val]
	xh = x[y>half_max_val]
	return np.sum(yh * xh) / np.sum(yh)

def parabola_peak(x, y, return_fit=False):
	# Fit a parabola to above half-max at 0.01 precision and return the peak
	half_max_val = np.max(y)/2
	yh = y[y>half_max_val]
	xh = x[y>half_max_val]
	abc = np.polyfit(xh, yh, deg=2)
	x_range = np.arange(np.around(np.min(xh), 2), np.around(np.max(xh), 2), 0.01)
	fitted_parabola = poly(x_range, *abc)
	parab_peak = x_range[fitted_parabola == np.max(fitted_parabola)][0]
	if return_fit:
		return parab_peak, (x_range, fitted_parabola)
	else:
		return parab_peak

def classic_mode(x, y):
	# Fit a spline, interpolate to 0.01 precision, and return the peak location
	half_max_val = np.max(y)/2
	spline = UnivariateSpline(x, y - half_max_val, s=0)
	r1, r2 = spline.roots()
	lo, hi = min(r1, r2), max(r1, r2)
	x_range = np.arange(np.around(lo, 2), np.around(hi, 2), 0.01)
	interp_y = spline(x_range)
	return x_range[interp_y == np.max(interp_y)]


def histogram_simple(img, mask):
	dhist, dedges = np.histogram(img[mask].ravel(), bins=BINS, range=hist_limits)
	prep_arr = lambda a, b: np.array([a, b]).T.flatten()
	histx, histy = prep_arr(dedges[:-1], dedges[1:]), prep_arr(dhist, dhist)
	bin_centers = (dedges[:-1]+dedges[1:])/2
	histx, histy = prep_arr(dedges[:-1], dedges[1:]), prep_arr(dhist, dhist)
	return (bin_centers, dhist), (histx, histy)


def plot_apply_mask(Xs_limit):
	values, histxy = histogram_simple(T, gen_chisq_mask(Xs_limit))
	bin_centers, dhist = values
	peak_val = np.max(dhist)
	peak_val_default_limits = [peak_val*0.95, peak_val*1.05]
	com = center_of_mass(*values)
	mode = classic_mode(*values)
	parab_peak, parab_values = parabola_peak(*values, return_fit=True)
	try:
		consensus = np.mean([com, mode, parab_peak])[0]
	except IndexError:
		consensus = np.mean([com, mode, parab_peak])		
	label_txt = r'$\chi^2$ < {:.2f} mask: T={:.2f}K'.format(
		Xs_limit, consensus)
	plt.plot(*histxy, '-', label=label_txt)
	plt.plot([com, com], peak_val_default_limits, '--', color='k')
	plt.plot([mode, mode], peak_val_default_limits, '--', color='k')
	plt.plot([parab_peak, parab_peak], peak_val_default_limits, '--', color='k')
	plt.plot(*parab_values, '--', color='b')

def final_temp_determination_plot():
	plt.figure(figsize=(11, 8.5))
	for Xs_limit in [5., 4., 3., 2., 1., 0.5, 0.25]:
		plot_apply_mask(Xs_limit)
	plt.legend()
	plt.xlabel("Fitted single-component temperature (K)")
	plt.ylabel("Histogram count")
	plt.title("Distribution of fitted temperatures under several single-component masks")
	plt.tight_layout()
	show_plot()

final_temp_determination_plot()
