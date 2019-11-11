import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

bandpass_beam_sizes = {
    "PACS160um": np.sqrt(11.64 * 15.65)/60,
    "SPIRE250um": 18.4/60,
    "SPIRE350um": 25.2/60,
    "SPIRE500um": 36.7/60,
    "F100": 9.66,
    "F143": 7.27,
    "F217": 5.01,
    "F353": 4.86,
    "F545": 4.84,
    "F857": 4.63,
}


def open_image(filename):
	# Opens a FITS file and returns (data, header)
	with fits.open(filename) as hdul:
		data = hdul[0].data
		header = hdul[0].header
	return data, header


def open_both_images(filename1, filename2):
	# Open both images, return the data and also a single WCS
	# Both images should be exactly on the same grid/WCS
	d1, h1 = open_image(filename1)
	d2, h2 = open_image(filename2)
	wcs_obj = WCS(h1)
	return d1, d2, wcs_obj


def prepare_convolution(w, beam_small, beam_large, data_shape):
	# Given a WCS object, two beam FWHMs in arcminutes (beam_small < beam_large), and image shape
	#  returns the Gaussian needed to bring an image at beam_small resolution down to beam_large resolution
	# Find pixel scale, in arcminutes
    dtheta_dpix_i = w.array_index_to_world(0, 0).separation(w.array_index_to_world(0, 1)).to('arcmin').to_value()
    dtheta_dpix_j = w.array_index_to_world(0, 0).separation(w.array_index_to_world(1, 0)).to('arcmin').to_value()
	dtheta_dpix_avg = (dtheta_dpix_i + dtheta_dpix_j)/2  # arcminutes
	# Using pixel scale, generate a grid for the Gaussian
	i, j = np.arange(data_shape[0]) - data_shape[0]//2, np.arange(data_shape[1]) - data_shape[1]//2
	i, j = i*dtheta_dpix_i, j*dtheta_dpix_j
	convolution_beam_width_sq = (beam_large*beam_large - beam_small*beam_small)/(2.35**2)
	i, j = np.exp(-i*i/(2*convolution_beam_width_sq)), np.exp(-j*j/(2*convolution_beam_width_sq))
	i, j = i/np.trapz(i), j/np.trapz(j)
	convolution_beam =  i[:, np.newaxis] * j[np.newaxis, :]
	return convolution_beam


def convolve_helper(image, kernel):
	# assumes no nans, just the convolution
	ft = np.fft.fft2(image)*np.fft.fft2(kernel)
	result = np.fft.ifft2(ft)
	return np.real(np.fft.fftshift(result))


def convolve_properly(image, kernel):
	# Need to preserve NaNs
	# as of July 22, 2019, corrects for edge effects, normalization from NaNs
	image = image.copy()
	nanmask = np.isnan(image)
	image[nanmask] = 0.
	result = convolve_helper(image, kernel)
	# account for edge effects, normalization
	image[~nanmask] = 1.
	norm = convolve_helper(image, kernel)
	image[:] = 1.
	norm /= convolve_helper(image, kernel)
	result /= norm
	result[nanmask] = np.nan
	return result


def quick_compare(*args):
	# Returns basic stats about the ratio
	if len(args) == 1:
		# assuming this is the ratio
		ratio = args[0]
	elif len(args) == 2:
		# assuming these are images
		ratio = args[0] / args[1]
	mean_r = np.nanmean(ratio)
	med_r = np.nanmedian(ratio)
	rms_r = np.nanstd(ratio)
	return mean_r, med_r, rms_r


def plot_compare(img1, img2, stub1, stub2):
	# Given two images already at the same resolution, plot a comparison
	#  and print out some useful statistics
	ratio = img2 / img1
	mean_r, med_r, rms_r = quick_compare(ratio)
	print("MEAN: %.3f, MEDIAN: %.3f, RMS: %.3f" % (mean_r, med_r, rms_r))
	vmin, vmax = 0.5, 3.5
	plt.figure(figsize=(19, 7))
	plt.subplot(131)
	plt.imshow(np.log10(img1), origin='lower', vmin=vmin, vmax=vmax)
	plt.title("log[%s]" % stub1)
	cbar = plt.colorbar()
	cbar.set_label("log10[MJy/sr]", rotation=270)
	plt.subplot(132)
	plt.imshow(ratio, origin='lower', vmin=0.8, vmax=1.2)
	plt.title("%s / %s" % (stub1, stub2))
	plt.text(0.01, 0.1, "Mean: %.3f\nMedian: %.3f\nRMS: %.3f" % (mean_r, med_r, rms_r),
		horizontalalignment='left', verticalalignment='center', transform=plt.gca().transAxes)
	cbar = plt.colorbar()
	cbar.set_label("MJy/sr", rotation=270)
	plt.subplot(133)
	plt.imshow(np.log10(img2), origin='lower', vmin=vmin, vmax=vmax)
	plt.title("log[%s]" % stub2)
	cbar = plt.colorbar()
	cbar.set_label("log10[MJy/sr]", rotation=270)
	plt.tight_layout()
	plt.show()


def polynomial(deg, fit, x):
	y = np.zeros(x.shape)
	for i, c in enumerate(fit):
		y += c * x**(deg - i)
	return y

def gaussian(x, mu, sigma, A):
	coeff = A/(np.sqrt(2 * np.pi)*sigma)
	exponent = -((x - mu)**2 / (2 * sigma*sigma))
	return coeff * np.exp(exponent)

def fit_gaussian(x, y, mu0=0, sigma0=1, A0=1):
	# Assumes histogram statistics
	return curve_fit(gaussian, x, y, p0=(mu0, sigma0, A0))

if __name__ == "__main__":
	from planck_mask import gen_hist_and_stats
	np.warnings.filterwarnings('ignore')
	per1_dir_stub = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
	# fn2 = "/n/sgraraid/filaments/data/HFI_SkyMap_857_resSPIRE350.fits"
	if len(sys.argv) > 1:
		band_stub = sys.argv[1]
		field_stub = sys.argv[2]
	else:
		band_stub = "SPIRE350um"
		field_stub = "Per1"
	fn1 = per1_dir_stub + "%s-image-remapped.fits" % band_stub
	fn2 = per1_dir_stub + "%s-image-PLANCKfabricated.fits" % band_stub

	img1_res = bandpass_beam_sizes[band_stub]
	img2_res = bandpass_beam_sizes['F857']
	img1, img2, w = open_both_images(fn1, fn2)

	conv_beam = prepare_convolution(w, img1_res, img2_res, img1.shape)
	img1 = convolve_properly(img1, conv_beam)
	#plot_compare(img1, img2, "Herschel-%s" % band_stub, "Planck")

	histxy, stats = gen_hist_and_stats(img2-img1, x_lim=(35, 55))
	# FIXME this is giving a bad peak value but is not well fit by a Gaussian
	# should be fit at FWHM and above by a Gaussian!
	
	popt, pcov = fit_gaussian(histxy[0], histxy[1]-np.max(histxy[1])/2)

	plt.figure(figsize=(11, 8.5))
	plt.subplot(111)
	plt.plot(*histxy, '-')
	A_approx = np.max(histxy[1])*np.sqrt(2 * np.pi)*stats[1]
	plt.plot(histxy[0], A_approx*gaussian(histxy[0], stats[0], stats[1], 1), '--')
	plt.plot(histxy[0], gaussian(histxy[0], *popt) + np.max(histxy[1])/2, '--')
	params_txt = "Mean: {0:.3f}\nSTD: {1:.5f}".format(stats[0], stats[1])
	plt.text(stats[0], np.max(histxy[1]), params_txt)
	plt.xlabel("Ratio Herschel real - fabricated")
	plt.ylabel("Histogram count")
	plt.show()

	# plt.title("Planck 857GHz obs / Herschel 350um obs")


	### FOR empirically finding the beam size of the dust model
"""
	stats_lists = [[], [], []]
	res_attempts = np.arange(2.1, 20, .5)
	for img2_res in res_attempts:
		img1, img2, w = open_both_images(fn1, fn2)
		conv_beam = prepare_convolution(w, img1_res, img2_res, img1.shape)
		img1 = convolve_properly(img1, conv_beam)
		current_stats = quick_compare(img1, img2)
		for s, l in zip(current_stats, stats_lists):
			l.append(s)
		print("Beam: %.1f" % img2_res, end="\r")
		# plot_compare(img1, img2, "857GHz", "Planck_fab")
	print()
	stats_lists = map(np.array, stats_lists)
	means, meds, rmss = stats_lists
	d = 6
	fit = np.polyfit(res_attempts, means, deg=d)
	x = np.arange(4.5, 20, 0.1)
	best_mean = x[np.argmin(polynomial(d, fit, x))]
	fit = np.polyfit(res_attempts, rmss, deg=d)
	best_rms = x[np.argmin(polynomial(d, fit, x))]
	print("[BEST] Mean: %.3f, RMS: %.3f" % (best_mean, best_rms))
	plt.subplot(311)
	plt.plot(res_attempts, means, '.')
	plt.title("mean (min~%.1f)" % best_mean)
	plt.subplot(312)
	plt.plot(res_attempts, meds, '.')
	plt.title("median")
	plt.subplot(313)
	plt.plot(res_attempts, rmss, '.')
	plt.title("rms (min~%.1f)" % best_rms)
	plt.show()
"""
