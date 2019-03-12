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
import sys

def show_plot():
	if plotting_remotely:
		plt.savefig(SAVE_NAME)
	else:
		plt.show()

# THREE PARAM || TWO PARAM
##
# 01 Tc				01 T
# 02 dTc			02 dT
# 03 Nc(H2)			03 N(H2)
# 04 dNc(H2)		04 dN(H2)
# 05 Th				05 Chi2/DOF
# 06 dTh			06 nIter
# 07 Nh(H2)			07 diff160
# 08 dNh(H2)		08 diff250
# 09 Chi2/DOF		09 diff350
# 10 nIter			10 diff500
# 11 diff160		11 BAND160
# 12 diff250		12 dBAND160
# 13 diff350		13 BAND250
# 14 diff500		14 dBAND250
# 15 BAND160		15 BAND350
# 16 dBAND160		16 dBAND350
# 17 BAND250		17 BAND500
# 18 dBAND250		18 dBAND500
# 19 BAND350
# 20 dBAND350
# 21 BAND500
# 22 dBAND500

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



fns_43 = (soln_2p_plus043, soln_3p_plus043)
fns_45 = (soln_2p_plus045, soln_3p_plus045)
fns_45_15p5 = (soln_2p_plus045, soln_3p_plus045_15p5)
fns_46 = (soln_2p_plus046, soln_3p_plus046)
fns_46_dT = (soln_3p_plus046_15p6, soln_3p_plus046, soln_3p_plus046_15p8)
fns_47 = (soln_2p_plus047, soln_3p_plus047)
fns_47_15p5 = (soln_2p_plus047, soln_3p_plus047_15p5)
fns_dPACS = (soln_3p_plus045, soln_3p_plus046, soln_3p_plus047)

nominal_2p_soln = soln_2p_plus046
nominal_3p_soln = soln_3p_plus046

frames_2p = (1, 3, 5)
frames_3p = (1, 3, 7)
def load_manticore(filename, frames=None):
	if frames is None:
		frames = (1, 3, 5)
	with fits.open(dir_per1 + filename) as hdul:
		imgs = tuple(hdul[x].data for x in frames)
	return imgs

def scatter(*args, x_lim=None, y_lim=None, log_x=False, log_y=False,
	labels=None, label=None):
	if len(args) == 3:
		img1, img2, mask = args
	elif len(args) == 2:
		img1, img2 = args
		mask = np.full(img1.shape, True)
	else:
		InputError("Too many/few inputs (%d) to scatter (2 or 3 ok)" % (len(args)))
	nanmask = ~np.isnan(img1) & ~np.isnan(img2)
	img1 = img1[mask & nanmask].ravel()
	img2 = img2[mask & nanmask].ravel()
	p = plt.plot(img1, img2, 'o', markersize=2, alpha=0.01, label=label)
	plt.xlim(x_lim)
	plt.ylim(y_lim)
	if log_x:
		plt.xscale('log')
	if log_y:
		plt.yscale('log')
	if labels is not None:
		x_label, y_label = labels
	else:
		x_label, y_label = "img1", "img2"
	plt.xlabel(x_label), plt.ylabel(y_label)
	return p

def histogram(*args, x_lim=None, log=False, label=None, setting=0):
	histxy, stats = gen_hist_and_stats(*args, x_lim=x_lim,
		log=log, setting=setting)
	p = plt.plot(*histxy, '-', label=label)
	# print(stats)
	plt.text(stats[0], np.max(histxy[1]),
		"Mode: {0:.3f}\n STD: {1:.3f}".format(stats[0], stats[1]))
	return p

def pacs_vs_spire_ratio():
	fluxes = tuple(fits.getdata(gen_flux_fn(x)) for x in range(2))
	ratio = fluxes[0] / fluxes[1] # PACS160 / SPIRE250
	plt.figure(figsize=(12, 9))
	plt.subplot(121)
	plt.imshow(ratio, origin='lower', vmin=0.8, vmax=1.2)
	plt.title("PACS160 / SPIRE250 ratio")
	plt.colorbar()
	plt.subplot(122)
	plt.imshow((ratio < 1), origin='lower')
	plt.title("Ratio < 1")
	plt.tight_layout()
	show_plot()


def cleaner_core_mask():
	mask = mtc.get_better_core()
	return get_mask(mask, n=4, dilation=3)


def cleaner_envelope_mask():
	mask = mtc.get_better_envelope()
	return get_mask(mask, n=4, dilation=1, min_size=400)


def cleaner_junk_mask():
	mask = ~(mtc.get_better_envelope() | mtc.get_better_core())
	return get_mask(mask, n=5, dilation=2)


def typical_SED(mask=None):
	fluxes = tuple(fits.getdata(gen_flux_fn(x)) for x in range(4))
	fluxes = tuple(x/fluxes[1] for x in fluxes)
	nanmask = ~np.any([np.isnan(x) | (x < 0) for x in fluxes], axis=0)
	if mask is not None:
		mask &= nanmask
	else:
		mask = nanmask
	fluxes = tuple(x[mask].flatten() for x in fluxes)
	stats = tuple((np.mean(x), np.median(x), np.std(x)) for x in fluxes)
	plt.figure(figsize=(12, 9))
	plt.subplot(111)
	for i, x, s in zip(range(len(fluxes)), fluxes, stats):
		plt.plot(np.full(x.size, i+1), x, 'o', color='k', markersize=1, alpha=0.01)
		plt.errorbar([i+1], [s[1]], yerr=[s[2]], fmt='o', color='r', capsize=10)
	plt.xlim((0.5, 4.5))
	plt.ylim((0, 1.5))
	plt.tight_layout()
	show_plot()


def hisgogram_typical_SED(mask=None):
	fluxes = tuple(fits.getdata(gen_flux_fn(x)) for x in range(4))
	norm_index = 1
	fluxes = tuple(x/fluxes[norm_index] for x in fluxes)
	nanmask = ~np.any([np.isnan(x) | (x < 0) for x in fluxes], axis=0)
	if mask is not None:
		mask &= nanmask
	else:
		mask = nanmask
	fluxes = tuple(x[mask].flatten() for x in fluxes)
	stats = tuple((np.mean(x), np.median(x), np.std(x)) for x in fluxes)
	plt.figure(figsize=(16, 9))
	axes = [plt.subplot(131 + i) for i in range(3)]
	ax_iter = iter(axes)
	for i in range(len(fluxes)):
		if i == norm_index:
			continue
		elif i == 0:
			x_lim = (0.1, 1.6)
		elif i == 2:
			x_lim = (0.5, 0.8)
		elif i == 3:
			x_lim = (0.2, 0.4)
		else:
			raise ValueError("What?")
		i_mean, i_median, i_std = stats[i]
		flux = fluxes[i]
		plt.sca(next(ax_iter))
		# core_mask = mtc.get_better_core()
		# envelope_mask = mtc.get_better_envelope()
		# junk_mask = ~(mtc.get_better_core() | mtc.get_better_envelope())
		core_mask = cleaner_core_mask()
		envelope_mask = cleaner_envelope_mask()
		junk_mask = cleaner_junk_mask()
		histogram(flux, core_mask[mask].flatten(), x_lim=x_lim, label="Filament")
		histogram(flux, envelope_mask[mask].flatten(), x_lim=x_lim, label="Envelope")
		histogram(flux, junk_mask[mask].flatten(), x_lim=x_lim, label="Single-T regions")
		plt.xlabel("Ratio to {} flux".format(herschel_flux[norm_index]))
		plt.ylabel("Histogram count")
		plt.ylim((0, 15000))
		plt.title(herschel_flux[i])
		plt.legend()
	plt.tight_layout()
	show_plot()


# Investigating the little bump in the 350/250 histogram
# Doesn't look important, should write that down somewhere
# flux_spire350 = fits.getdata(gen_flux_fn(2))
# flux_spire250 = fits.getdata(gen_flux_fn(1))
# ratio = flux_spire350 / flux_spire250
# quick_mask = lambda x: (x < 0)
# ratio[quick_mask(flux_spire350) | quick_mask(flux_spire250)] = np.nan
# mask = (ratio < 0.55)# & cleaner_junk_mask()
# flux_spire350[~mask] = np.nan
# plt.imshow(flux_spire350, origin='lower', vmax=150)
# plt.colorbar()
# show_plot()

def quickrun_dT_hists():
	imgs_dT = list(load_manticore(x, frames=frames_3p) for x in fns_46_dT)
	imgs_center = imgs_dT.pop(1)
	plt.figure(figsize=(15, 10))
	axes = [plt.subplot(231 + i) for i in range(2*len(frames_3p))]
	core_mask = cleaner_core_mask()
	envelope_mask = cleaner_envelope_mask()
	lims = [(0.8, 1.2), (0.6, 1.4), (0.9, 1.1)]
	for i, frame_name in enumerate(["Tc", "Nc", "Nh"]):
		for j in range(2):
			ratio = imgs_dT[j][i] / imgs_center[i]
			txt = ("+" if j else "-") + "0.1K"
			plt.sca(axes[i])
			histogram(ratio.flatten(), core_mask.flatten(),
				x_lim=lims[i], label=txt)
			plt.title("Filament "+frame_name)
			plt.sca(axes[i + 3])
			histogram(ratio.flatten(), envelope_mask.flatten(),
				x_lim=lims[i], label=txt)
			plt.title("Envelope "+frame_name)
	for ax in axes:
		plt.sca(ax)
		plt.legend()
		plt.tight_layout()
	show_plot()

def quickrun_dPACS_hists():
	imgs_dPACS = list(load_manticore(x, frames=frames_3p) for x in fns_dPACS)
	imgs_center = imgs_dPACS.pop(1)
	plt.figure(figsize=(15, 10))
	axes = [plt.subplot(231 + i) for i in range(6)]
	core_mask = mtc.get_better_core()
	envelope_mask = mtc.get_better_envelope()
	lims = [(0.97, 1.03), (0.9, 1.1), (0.98, 1.02)]
	for i, frame_name in enumerate(["Tc", "Nc", "Nh"]):
		for j in range(2):
			ratio = imgs_dPACS[j][i] / imgs_center[i]
			txt = ("+" if j else "-") + "1MJy/sr"
			plt.sca(axes[i])
			histogram(ratio.flatten(), core_mask.flatten(),
				x_lim=lims[i], label=txt)
			if j:
				plt.title("Filament "+frame_name)
				plt.legend()
			plt.sca(axes[i + 3])
			histogram(ratio.flatten(), envelope_mask.flatten(),
				x_lim=lims[i], label=txt)
			if j:
				plt.title("Envelope "+frame_name)
				plt.legend()
	plt.tight_layout()
	show_plot()


def quickrun_examine_region_values_hist():
	imgs = load_manticore(soln_3p_plus046, frames=frames_3p)
	plt.figure(figsize=(11, 8.5))
	axes = [plt.subplot(131 + i) for i in range(3)]
	core_mask = mtc.get_better_core()
	envelope_mask = mtc.get_better_envelope()
	filament_mask = mtc.get_filament_mask()
	masks = [core_mask, envelope_mask, filament_mask]
	lims = [(3, 17), (19.5, 23), (19.5, 23)]
	mask_names = ["Old Filament", "Envelope", "New Filament"]
	for i, frame_name in enumerate(["Tc", "Nc", "Nh"]):
		for j in range(3):
			mask_name = mask_names[j]
			plt.sca(axes[i])
			value = imgs[i]
			if i > 0:
				value = np.log10(value)
			histogram(value.flatten(), masks[j].flatten(), x_lim=lims[i],
				label=mask_name)
			if j == 2:
				plt.legend()
				plt.title(frame_name)
				plt.xlabel("Value")
				plt.ylabel("Count")
	plt.tight_layout()
	show_plot()

def quickrun_examine_region_values_scatter():
	imgs = load_manticore(soln_3p_plus046, frames=frames_3p)
	plt.figure(figsize=(16, 8))
	axes = [plt.subplot(131 + i) for i in range(3)]
	core_mask = mtc.get_better_core()
	envelope_mask = mtc.get_better_envelope()
	masks = [core_mask, envelope_mask]
	lims = [(3, 17), (19.5, 23), (19.5, 23)]
	Nc = np.log10(imgs[1])
	for i, frame_name in enumerate(["Tc vs Nc", "Nh vs Nc"]):
		for j in range(2):
			mask_name = "Envelope" if j else "Core"
			plt.sca(axes[i])
			value = imgs[i*2]
			if i > 0:
				value = np.log10(value)
			scatter(Nc.flatten(), value.flatten(), masks[j].flatten(),
				y_lim=lims[i*2], x_lim=lims[1], label=mask_name,
				labels=frame_name.split(' vs ')[::-1])
			if frame_name == "Nh vs Nc":
				plt.plot(list(plt.xlim()), list(plt.xlim()), '--', alpha=0.5)
			if j:
				plt.legend()
				plt.title(frame_name)
	frame_name = "Tc vs Nh"
	for j in range(2):
		mask_name = "Envelope" if j else "Core"
		plt.sca(axes[2])
		Tc = imgs[0]
		Nh = np.log10(imgs[2])
		scatter(Nh.flatten(), Tc.flatten(), masks[j].flatten(),
			y_lim=lims[0], x_lim=lims[2], label=mask_name,
			labels=frame_name.split(' vs ')[::-1])
		if j:
			plt.legend()
			plt.title(frame_name)
	plt.tight_layout()
	show_plot()

def quickrun_3param_improvement():
	imgs2 = load_manticore(soln_2p_plus046, frames=frames_2p)
	imgs3 = load_manticore(soln_3p_plus046, frames=frames_3p)
	T, N, Xs = imgs2
	Tc, Nc, Nh = imgs3
	masks = mtc.get_better_core(), mtc.get_better_envelope()
	Xs_effmask = ~np.isnan(Tc)
	leftover_mask = Xs_effmask & ~np.any(masks, axis=0)
	# what_is_this_mask = (T > 15.7) & (Tc < 8) # interesting
	what_is_this_mask = (Nc < 10**21.4) & (Nc > 0)
	# what_is_this_mask = (Nc/N < 1) & (Nc > 10**21.25) # also interesting
	masks = (leftover_mask, masks[0], masks[1], what_is_this_mask)
	mask_labels = ["Single-T leftover", "Filament", "Envelope", "What"]
	plt.figure(figsize=(14, 9))
	axes = [plt.subplot(131 +i) for i in range(3)]
	Tclim = (3, 16)
	Tlim = (10, 20)
	Nlim = (19, 23)
	lims = [(Tlim, Tclim)] + [(Nlim, Nlim)]*2
	labels = [
		("T", "Tc"),
		("N", "Nc"),
		("N", "Nh")
	]
	for i in range(3):
		val2 = np.log10(N) if i else T
		val3 = imgs3[i]
		if i:
			val3 = np.log10(val3)
		plt.sca(axes[i])
		legend1 = []
		lim_2p, lim_3p = lims[i]
		x_range = np.linspace(*lim_2p, 10)
		plt.plot(x_range, x_range, '--', label="Unchanged",
			linewidth=2, alpha=0.3, color='k')
		for j in range(len(masks)):
			p = scatter(val2.flatten(), val3.flatten(), masks[j].flatten(),
				x_lim=lim_2p, y_lim=lim_3p, labels=labels[i])
			legend1.append(Patch(color=p[0].get_color(), label=mask_labels[j]))
		plt.legend(handles=legend1)

	plt.tight_layout()
	show_plot()


def is_original_chisq_worse_in_envelope():
	# Interesting -- the 1-T Xs distributions are ~the same!!
	Xs = load_manticore(soln_2p_plus046, frames=(5,))[0]
	Xs2 = load_manticore(soln_3p_plus046, frames=(9,))[0]
	masks = (mtc.get_better_core(), mtc.get_better_envelope())
	names = ("Filament", "Envelope")
	plt.figure(figsize=(11, 8.5))
	plt.subplot(111)
	lim = (0, 15)
	values = (Xs, Xs2)
	tags = ("1-T", "2-T")
	for v, tag in zip(values, tags):
		for mask, name in zip(masks, names):
			histogram(v.flatten(), mask.flatten(), x_lim=lim, log=0,
				label="{} ({})".format(name, tag))
	plt.legend()
	plt.title("Chi Squared from 1-Temp model")
	plt.xlabel("Chi Squared")
	plt.ylabel("Histogram count")
	plt.tight_layout()
	show_plot()


def visualize_error_fraction():
	# Cool constructs that don't take up any space!!!1!
	imgs = lambda : load_manticore(soln_3p_plus046, frames=(16, 18, 20, 22)) # the error maps
	norms = lambda : load_manticore(soln_3p_plus046, frames=(15, 17, 19, 21)) # the actual maps
	normed_images = tuple(i/n for i, n in zip(imgs(), norms()))
	i = 0
	plt.figure(figsize=(11, 8.5))
	for n in normed_images:
		plt.subplot(221 + i)
		plt.imshow(n*100, origin='lower', vmin=0 if i else 5, vmax=50 if not i else 5)
		plt.colorbar()
		plt.title(herschel_flux[i])
		i += 1
	plt.tight_layout()
	show_plot()

def check_error_maps():
	imgs = load_manticore(soln_3p_plus046, frames=(16, 18, 20, 22)) # the error maps
	norms = load_manticore(soln_3p_plus046, frames=(15, 17, 19, 21)) # the actual maps
	masks = (mtc.get_better_core(), mtc.get_better_envelope())
	names = ("Filament", "Envelope")
	plt.figure(figsize=(11, 8.5))
	plt.subplot(111)
	lims = ((0, 60), (0, 6))
	axes = (plt.subplot(221 + i) for i in range(4))
	for i in range(4):
		plt.sca(next(axes))
		value = 100*imgs[i] / norms[i]
		stub = herschel_flux[i]
		lim = lims[int(bool(i))]
		for mask, name in zip(masks, names):
			histogram(value.flatten(), mask.flatten(), x_lim=lim,
				label="{}".format(name))
		plt.legend()
		plt.title(stub)
		plt.xlabel("Stated Flux Error {%}")
		plt.ylabel("Histogram count")
	plt.tight_layout()
	show_plot()


def check_difference_maps():
	imgs1 = load_manticore(soln_3p_plus046, frames=(11, 12, 13, 14)) # the difference maps
	imgs2 = load_manticore(soln_2p_plus046, frames=(7, 8, 9, 10)) # the difference maps
	norms = load_manticore(soln_3p_plus046, frames=(15, 17, 19, 21)) # the actual maps
	masks = (mtc.get_better_core(), mtc.get_better_envelope())
	names = ("Filament", "Envelope")
	tags = ("2-T", "1-T")
	imgs = (imgs1, imgs2)
	plt.figure(figsize=(11, 8.5))
	plt.subplot(111)
	axes = (plt.subplot(221 + i) for i in range(4))
	lims = ((-20, 7), (-7, 7))
	for i in range(4):
		plt.sca(next(axes))
		lim = lims[int(bool(i))]
		for tag, img in zip(tags, imgs):
			value = 100*img[i] / norms[i]
			for mask, name in zip(masks, names):
				histogram(value.flatten(), mask.flatten(), x_lim=lim,
					label="{} ({})".format(name, tag))
		plt.legend()
		plt.title(herschel_flux[i])
		plt.xlabel("Difference")
		plt.ylabel("Histogram count")
	plt.tight_layout()
	show_plot()



def plot_SED(x, y, residuals=False):
	# INPUT IN DS9 CONVENTION (1 indexed)
	i = y - 1
	j = x - 1
	imgs1 = load_manticore(soln_3p_plus046, frames=(11, 12, 13, 14)) # the difference maps
	imgs2 = load_manticore(soln_2p_plus046, frames=(7, 8, 9, 10)) # the difference maps
	norms = load_manticore(soln_3p_plus046, frames=(15, 17, 19, 21)) # the actual maps
	herschel_wavelengths = np.array([160., 250., 350., 500.])
	values = [np.array([img[i, j] for img in img_set]) for img_set in (imgs1, imgs2, norms)]
	observed_values = values.pop()
	tags = ("2-T", "1-T")
	count = 5
	for tag, diffs in zip(tags, values):
		if residuals:
			model_values = np.zeros(4)
		else:
			model_values = diffs + observed_values
		error_bars = []
		for d in diffs:
			if d > 0:
				error_bars.append([d, 0])
			else:
				error_bars.append([0, -d])
		error_bars = np.array(list(map(list, zip(*error_bars))))
		plt.errorbar(herschel_wavelengths, model_values, yerr=error_bars,
			fmt='D', markersize=3, capsize=count, label=tag)
		count = None
	plt.legend()
	plt.title("SED of pixel X({:d}), Y({:d})".format(int(x), int(y)))
	plt.xlabel("Wavelength (micron)")
	plt.ylabel("Flux (MJy/sr)")
	plt.xscale('log')
	plt.xlim((100, 600))
	if residuals:
		plt.ylim((-20, 20))
	else:
		plt.ylim((0, 300))
	plt.tight_layout()
	show_plot()

def examine_2p_chisq_hist():
	plt.figure(figsize=(11, 8.5))
	pacs_mask = mtc.get_pacs_mask()
	T, Xs = load_manticore(soln_2p_plus046, frames=(1, 5))
	# Xsrange = 10**np.linspace(-3, 0.7, 30)
	# extra_masks = {"Xs<{:.1f}".format(x): x for x in Xsrange}
	# ax = plt.subplot(121)
	extra_masks = {"": pacs_mask,
		"not filament": ~mtc.get_filament_mask(),
		"Xs<1": Xs < 1,
	}
	lims = (3, 25)
	lim_list = []
	modes, sigmas = [], []
	mode_err, sig_err = [], []
	for m in extra_masks:
		histxy, stats = gen_hist_and_stats(T,
			((extra_masks[m]) & pacs_mask),
			x_lim=lims,
			setting=-2
		)
		print(stats)
		# lim_list.append(extra_masks[m])
		# modes.append(stats[0])
		# sigmas.append(stats[1])
		# e_m, e_s, e_A = stats[2]
		# mode_err.append(e_m)
		# sig_err.append(e_s)
		histogram(T,
			((extra_masks[m]) & pacs_mask),
			x_lim=lims,
			label="M:{}".format(m),
			setting=-1,
		)
	plt.title("T distribution from single-T")
	plt.xlabel("T")
	plt.ylabel("Histogram Count")
	plt.legend()
	# plt.subplot(222)
	# plt.plot(lim_list, mode_err, '.')
	# plt.subplot(224)
	# plt.plot(lim_list, sig_err, '.')
	plt.tight_layout()
	show_plot()


def big_change_in_T():
	fns = mtc.fns_bigchangeT
	plt.figure(figsize=(15, 12))
	axes = iter([plt.subplot(231 + i) for i in range(6)])
	assert len(fns) <= 6
	new_temperatures = sorted([*fns])
	mask = masking_attempt(filename_override=None)
	for temperature in new_temperatures:
		filename = fns[temperature]
		plt.sca(next(axes))
		# mtc.quickrun_image_masked_full(3, mask, l=(21., 22.5),
		mtc.quickrun_image_masked_full(1, mask, l=(5, 11),
			filename_override=filename, ax=1)
		plt.title(f"T = {temperature} K")
	plt.tight_layout()
	show_plot()


def definitive_hotT_method():
	"""
	Gaussian fit the Xs<1-masked single-T distribution
	and return the mean as the hot T
	This method is:
	1) very robust to the histogram limits, binning
		answer varies by less than 0.003K
	2) very robust to the mask used (i.e. Xs < what)
		answer varies by about 0.05K
	3) will give us the hot T to 0.01K, conservatively
		Though could claim down to 0.003K via the fit
	Since we know the 3-parameter fit is sensitive to
	changes in assumed T at the 0.1K level, 0.01K confidence
	should put our minds at ease

	The hot T has statistical error associated with it
	from the sample; this is impossible to correct
	That error is nearly ~1K
	This method only adds ~0.01K, so I'd approximate
	it as "perfect" within the existing statistical error
	"""
	T, Xs = load_manticore(soln_2p_plus046, frames=(1, 5))
	hist_limits = (10, 20)
	histxy, stats = gen_hist_and_stats(T,
		(mtc.get_pacs_mask() & (Xs < 1)),
		x_lim=hist_limits,
		setting=-2
	)
	hot_T = stats[0]
	err_hotT = stats[1]
	return hot_T, err_hotT

def masking_attempt(max_dT=2, n=4,
	filename_override=None):
	"""
	The mask is built on:
	1) requiring all 4 bands to be present
	2) constraining densities to be positive
	3) constraining dT to be lower than a certain value
	It is then cleaned with the boolean islands technique
	"""
	mask = mtc.get_notjunk_mask() #& mtc.get_filament_mask()
	mask &= mtc.mask_img_full(2, (0, max_dT),
		filename_override=filename_override)
	# ADDED A dNc MASK, it works very well to clean out the remaining
	# problem regions
	mask &= mtc.mask_img_full(4, (0, 10**22),
		filename_override=filename_override)
	# mask = fill_inwards(mask, mtc.get_pacs_mask(), n=n, min_size=3)
	return mask

def maxdTvsN():
	n_range = [2, 4, 6, 8]
	max_dT_range = [1.5, 2, 3, 5]
	notjunk_mask = mtc.get_notjunk_mask()
	pacs_mask = mtc.get_pacs_mask()
	img = mtc.load_specific_frame(soln_3p_plus046f, 1)
	plt.figure(figsize=(16, 16))
	axes = [plt.subplot(4, 4, 1 + x) for x in range(16)]
	for j, max_dT in enumerate(max_dT_range):
		mask = notjunk_mask & mtc.mask_img_full(2, (0, max_dT))
		for i, n in enumerate(n_range):
			mask = fill_inwards(mask, pacs_mask, n=n)
			plt.sca(axes[i*4 + j])
			img_copy = img.copy()
			img_copy[~mask] = np.nan
			plt.imshow(img_copy, origin='lower',
				vmin=5, vmax=6)
			if j == 0:
				plt.ylabel(f'n = {n}')
			if i == 3:
				plt.xlabel(f'dT limit = {max_dT}')
			print(f'finished n={n}, max_dT={max_dT}')
	plt.tight_layout()
	show_plot()

# plot_SED(531, 308, residuals=False)
np.warnings.filterwarnings('ignore')
# print(np.around(definitive_hotT_method(), 2))

# mask = masking_attempt(2, 1)
# Cold temp
#mtc.quickrun_image_masked_full(1, mask, l=(5, 11))
# Cold column
# mtc.quickrun_image_masked_full(7, mask, l=(20.5, 21.5))
big_change_in_T()
