plotting_remotely = True
import numpy as np
from astropy.io.fits import getdata, writeto
from astropy.wcs import WCS
import matplotlib
if plotting_remotely:
	matplotlib.use('Agg')
SAVE_NAME = "Figure_X_current.png"
import matplotlib.pyplot as plt
from planck_mask import gen_hist_and_stats, get_referenced_mask, get_spire_mask
import manticore_results as mtc
from compare_images import prepare_convolution, convolve_properly, gaussian
from AA_planck_obs  import DustModel
from AA_regrid_healpy import project_healpix_to_fits, open_healpix

def show_plot():
	if plotting_remotely:
		plt.savefig(SAVE_NAME)
	else:
		plt.show()

XsMAX, XsMIN = 1.e0, 1.e-3
XsMOD = 2
XsMASK = lambda x: (x < XsMAX) & (x > XsMIN)
BINS = 128
diffLIM = (25, 80)
small_diffLIM = (-10, 10)
fluxLIM = (-50, 80)
XsPLOTLIM = (-5, 3.1)

img_stub = "-image-remapped.fits"
pacs_band_stub = "PACS160um"
region_directory = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
ref_stub = "SPIRE350um"
dust_directory = "/n/sgraraid/filaments/data/TEST4/regridding_stuff/"
planck_data_directory = "/n/sgraraid/filaments/data/"
planck_stubs = {
	'F353': "HFI_SkyMap_353-psb-field-IQU_2048_R3.00_full.fits",
	'F545': "HFI_SkyMap_545-field-Int_2048_R3.00_full.fits",
	'F857': "HFI_SkyMap_857-field-Int_2048_R3.00_full.fits",
}
planck_unit_conversions = {
    "F100": 244.1,
    "F143": 371.74,
    "F217": 483.690,
    "F353": 287.450,
    "F545": 1,
    "F857": 1,
}
gen_planck_fn = lambda b: planck_data_directory + planck_stubs[b]
gen_herschel_fn = lambda b: region_directory + b + img_stub

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
GNILC_resolution = 5 # arcminutes


def load_planck_band(planck_band_stub):
	# Regrid a Planck HFI band to the -remapped.fits pixel grid
	# planck_band_stub must be like 'F857'
	assert (planck_band_stub[0] == 'F') and (planck_band_stub in planck_stubs)
	reference_filename = gen_herschel_fn(ref_stub)
	ref_data, ref_header = getdata(reference_filename, header=True)
	planck_fn = gen_planck_fn(planck_band_stub)
	regridded_planck_data = project_healpix_to_fits(
		open_healpix(planck_fn, nest=True),
		ref_data, ref_header,
		nest=True
	)
	unit_conversion = planck_unit_conversions[planck_band_stub]
	return regridded_planck_data * unit_conversion

def load_herschel_band(herschel_band_stub, header=True):
	herschel_filename = gen_herschel_fn(herschel_band_stub)
	return_value = getdata(herschel_filename, header=header)
	return return_value

def generate_planck_flux(band_stub):
	# Reusable; does the regridding & GNILC dust model stuff, returns
	# the final Planck based image
	# band_stub can be anything in the bandpass_files dict
	# pulls spire350-remapped for a reference
	ref_data, ref_header = load_herschel_band(ref_stub)
	sky = DustModel("Per1", band_stub, ref_data, ref_header, dust_directory)
	return sky.observe_planck(band_stub)

def get_planck_ratio(planck_band_stub):
	assert (planck_band_stub[0] == 'F') and (planck_band_stub in planck_stubs)
	# get the ratio between the actual planck h_data
	# and the predicted dust model emission for that band
	ratio_filename = "{}Planck_ObsToModel_ratio_{}.fits".format(
		dust_directory, planck_band_stub
	)
	return getdata(ratio_filename)
	# the hard work has been done but is here for reference
	observed_planck_flux = load_planck_band(planck_band_stub)
	predicted_planck_flux = generate_planck_flux(planck_band_stub)
	return observed_planck_flux / predicted_planck_flux

def write_all_planck_masks():
	raise RuntimeError("Dont use this function anymore, you did already")
	ref_data, ref_header = load_herschel_band(ref_stub)
	ref_header['COMMENT'] = "NOT SPIRE -- PLANCK RATIO (March 12 2019)"
	for planck_band_stub in planck_stubs:
		ratio = get_planck_ratio(planck_band_stub)
		save_fn = "Planck_ObsToModel_ratio_{}.fits".format(planck_band_stub)
		writeto(save_fn, ratio, ref_header)
	return 0

def get_planck_ratio_mask(planck_band_stub, setting=4):
	ratio = get_planck_ratio(planck_band_stub)
	return get_referenced_mask(ratio, setting=setting)

def plot_cumulative_planck_ratio_mask(ax):
	masks = []
	settings = [4, 5]
	for s in settings:
		for planck_band_stub in planck_stubs:
			masks.append(get_planck_ratio_mask(planck_band_stub, setting=s))
	int_masks = [m.astype(int) for m  in masks]
	cumulative_mask = sum(int_masks)
	plt.sca(ax)
	plt.imshow(cumulative_mask, origin='lower')
	plt.colorbar()
	plt.title("Accumulated Planck Observation vs Dust Model Masks")

def accumulate_planck_ratio_masks():
	masks = []
	for planck_band_stub in planck_stubs:
		masks.append(get_planck_ratio_mask(planck_band_stub))
	masks = np.array(masks)
	final_mask = np.all(masks, axis=0)
	return final_mask

def get_pacs_difference_map():
	# Specific to PACS
	h_data, h_header = load_herschel_band(pacs_band_stub)
	planck_flux = generate_planck_flux(pacs_band_stub)

	# Convolve the Herschel map down to Planck resolution
	h_beam, p_beam = bandpass_beam_sizes[pacs_band_stub], GNILC_resolution  # arcminutes
	conv_beam = prepare_convolution(WCS(h_header), h_beam, p_beam, h_data.shape)
	h_data = convolve_properly(h_data, conv_beam)

	diff = planck_flux - h_data
	return diff


def get_pacs_offset():
	diff = get_pacs_difference_map()
	mask = accumulate_planck_ratio_masks()
	histxy, stats = gen_hist_and_stats(diff, mask, x_lim=diffLIM, setting=-1)
	return stats[0]

def plot_pacs_offset(ax):
	diff = get_pacs_difference_map()
	mask = accumulate_planck_ratio_masks()
	histxy, stats = gen_hist_and_stats(diff, mask, x_lim=diffLIM, setting=-1)
	fitxy, unused = gen_hist_and_stats(diff, mask, x_lim=diffLIM, setting=-3)
	plt.sca(ax)
	plt.plot(*histxy, '-', color='k')
	plt.plot(*fitxy, '--', color='r')
	print(stats)

def plot_multiple_offset_attempts(ax):
	settings = {4: "within 10%", 5: "within 20%"}
	colors = iter(['red', 'blue', 'orange', 'magenta', 'brown', 'green'])
	diff = get_pacs_difference_map()
	plt.sca(ax)
	for setting in settings:
		for planck_band_stub in planck_stubs:
			mask = get_planck_ratio_mask(planck_band_stub, setting=setting)
			histxy, stats = gen_hist_and_stats(diff, mask, x_lim=diffLIM, setting=-1)
			fitxy, unused = gen_hist_and_stats(diff, mask, x_lim=diffLIM, setting=-3)
			label = "{} GHz band ratio {}: offset = {:.1f}".format(planck_band_stub[1:], settings[setting], stats[0])
			color = next(colors)
			plt.plot(*histxy, '-', color=color, label=label)
			plt.plot(*fitxy, '--', color=color)
	plt.legend()
	plt.title("PACS zero-point offset distribution")
	plt.ylabel("Histogram count")
	plt.xlabel("Additive offset to PACS 160 micron (MJy/sr)")

def plot_before_and_after(ax):
	h_data = load_herschel_band(pacs_band_stub, header=False)
	planck_flux = generate_planck_flux(pacs_band_stub)
	corrected_h_data = h_data + get_pacs_offset()
	images = {
		"Uncorrected PACS 160 micron flux": (h_data, 'r-'),
		"Planck-predicted PACS 160 micron flux": (planck_flux, 'k--'),
		"Corrected PACS 160 micron flux": (corrected_h_data, 'b-'),
	}
	plt.sca(ax)
	for label in images:
		img, marker = images[label]
		histxy, stats = gen_hist_and_stats(img, ~np.isnan(img), x_lim=fluxLIM)
		plt.plot(*histxy, marker, label=label)
	plt.title("PACS 160 micron flux distribution")
	plt.xlabel("Flux (MJy/sr)")
	plt.ylabel("Histogram count")
	plt.legend()


def capstone_figure_pacs_offset():
	plt.figure(figsize=(17, 7))
	ax = plt.subplot(121)
	plot_multiple_offset_attempts(ax)
	ax = plt.subplot(122)
	plot_cumulative_planck_ratio_mask(ax)
	plt.tight_layout()
	show_plot()

def before_and_after_figure():
	plt.figure(figsize=(11.5, 8))
	ax = plt.subplot(111)
	plot_before_and_after(ax)
	plt.tight_layout()
	show_plot()

before_and_after_figure()
