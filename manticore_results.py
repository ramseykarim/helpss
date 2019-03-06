plotting_remotely = False
import numpy as np
import matplotlib
if plotting_remotely:
	matplotlib.use('Agg')
import matplotlib.pyplot as plt
SAVE_NAME = "Figure_X_current.png"
from astropy.io import fits
from planck_mask import gen_hist_and_stats

def show_plot():
	if plotting_remotely:
		plt.savefig(SAVE_NAME)
	else:
		plt.show()

dir_per1 = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"

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
soln_3p_plus046f_15p8 = "T4-absdiff-Per1J-3param-plus046-full-15p8.fits"
soln_3p_plus047 = "T4-absdiff-Per1J-3param-plus047.fits"
soln_3p_plus047_15p5 = "T4-absdiff-Per1J-3param-plus047-15p5.fits"

fns_43 = (soln_2p_plus043, soln_3p_plus043)
fns_45 = (soln_2p_plus045, soln_3p_plus045)
fns_45_15p5 = (soln_2p_plus045, soln_3p_plus045_15p5)
fns_46 = (soln_2p_plus046, soln_3p_plus046)
fns_47 = (soln_2p_plus047, soln_3p_plus047)
fns_47_15p5 = (soln_2p_plus047, soln_3p_plus047_15p5)

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


def load_2param_TNX(filename):
	with fits.open(filename) as hdul:
		T = hdul[1].data
		N = hdul[3].data
		Xs = hdul[5].data
	return T, N, Xs


def load_3param_TcThNcNh(filename):
	with fits.open(filename) as hdul:
		Tc = hdul[1].data
		Nc = hdul[3].data
		Th = hdul[5].data
		Nh = hdul[7].data
	return Tc, Th, Nc, Nh


def load_specific_frame(filename, frame):
	return fits.getdata(dir_per1 + filename, frame)


def load_specific_frames(filename, frames=None):
	if frames is None:
		raise RuntimeError("Need frame numbers")
	with fits.open(dir_per1 + filename) as hdul:
		imgs = tuple(hdul[x].data for x in frames)
	return imgs


def load_relevant_frames(soln2_fn, soln3_fn):
	Tsingle, Nsingle, Xsingle = load_2param_TNX(dir_per1+soln2_fn)
	Tc, Th, Nc, Nh = load_3param_TcThNcNh(dir_per1+soln3_fn)
	# Th[np.isnan(Th)] = Tsingle[np.isnan(Th)]
	# Nh[np.isnan(Nh)] = Nsingle[np.isnan(Nh)]
	return Tsingle, Nsingle, Tc, Th, Nc, Nh


field_names = ["T_single", "NH2_single", "Tc", "Th", "Nc", "Nh"]


def get_pacs_mask():
	return ~np.isnan(load_specific_frame(soln_2p_plus046, 11))


def get_cold_mask():
	base_Nc = load_specific_frame(soln_3p_plus043, 3)
	ratio_Nc = load_specific_frame(soln_3p_plus046, 3) / base_Nc
	mask = get_pacs_mask() & ~((base_Nc < 0) | (ratio_Nc < 0)) & ~np.isnan(ratio_Nc)
	return mask


def get_core_mask():
	ratio_Nc = load_specific_frame(soln_3p_plus046, 3) / load_specific_frame(soln_3p_plus043, 3)
	ratio_Nc[~get_cold_mask()] = np.nan
	return (ratio_Nc > 1)


def get_envelope_mask():
	ratio_Nc = load_specific_frame(soln_3p_plus046, 3) / load_specific_frame(soln_3p_plus043, 3)
	ratio_Nc[~get_cold_mask()] = np.nan
	mask = (ratio_Nc < 1)
	base_Tc = load_specific_frame(soln_3p_plus043, 1)
	ratio_Tc = load_specific_frame(soln_3p_plus046, 1) / base_Tc
	ratio_Tc[~mask] = np.nan
	return (ratio_Tc > 1.004) & (ratio_Tc < 1.02)


def get_mystery_envelope_mask():
	return get_cold_mask() & ~(get_core_mask() | get_envelope_mask())

def get_ratios_to_original(filenames, original=fns_43):
	original = load_relevant_frames(*original)
	newer = load_relevant_frames(*filenames)
	return tuple(map(lambda x: x[0]/x[1], zip(newer, original)))

lims_no_off = [
	(0.9, 1.2),
	(0.6, 1.1),
	(.96, 1.04),
	(.99, 1.01),
	(.85, 1.15),
	(.997, 1.005),
]


def quickrun_compare_no_offset():
	no_offset = load_2param_TNX(dir_per1 + soln_2p_plus000)
	offsets = (load_2param_TNX(dir_per1 + x) for x in (soln_2p_plus043, soln_2p_plus046, soln_2p_plus047))
	ratios = ((x[i]/no_offset[i] for i in range(2)) for x in offsets)
	pacs_mask = get_pacs_mask()
	labels = ['+43', '+46', '+47']
	plt.figure()
	axes = [plt.subplot(x) for x in (121, 122)]
	ratio_set_i = 0
	for ratio_set in ratios:
		ratio_set_j = 0
		for ratio in ratio_set:
			mask = ~np.isnan(ratio) & pacs_mask
			histxy, stats = gen_hist_and_stats(ratio, mask, x_lim=lims_no_off[ratio_set_j])
			plt.sca(axes[ratio_set_j])
			plt.plot((histxy[0] - 1)*100, histxy[1], '-',
				label="{} / 00".format(labels[ratio_set_i]))
			if ratio_set_i == 0:
				plt.title(field_names[ratio_set_j])
				plt.xlabel("Percent change")
				plt.ylabel("Histogram count")
			ratio_set_j += 1
		ratio_set_i += 1
	plt.legend()
	show_plot()


lims = [
	(1, 1.01),
	(.97, 1),
	(.96, 1.04),
	(.99, 1.01),
	(.85, 1.15),
	(.997, 1.005),
]


def quickrun_histograms():
	ratio_sets = (get_ratios_to_original(x) for x in (fns_45, fns_46, fns_47))
	pacs_mask = get_pacs_mask()
	set_labels = ['+45', '+46', '+47']
	plt.figure()
	for ratio_set_i in range(3):
		ratio_set = next(ratio_sets)
		for i in range(len(ratio_set)):
			ratio = ratio_set[i]
			mask = ~np.isnan(ratio) & pacs_mask
			histxy, stats = gen_hist_and_stats(ratio, mask, x_lim=lims[i])
			plt.subplot(231 + i)
			plt.plot((histxy[0] - 1)*100, histxy[1], '-',
				label=set_labels[ratio_set_i]+"/43")
			if ratio_set_i == 0:
				name = field_names[i]
				plt.title(name+"NORMAL")
				plt.xlabel("Percent change")
				plt.ylabel("Histogram count")
			if i == 3:
				plt.legend()
	show_plot()

def quickrun_histograms_454647(mask=None):
	original_fn = fns_45
	new_fn = (fns_46, fns_47)
	new = (get_ratios_to_original(x, original=original_fn) for x in new_fn)
	relevant_field_names = ('Tc', 'Nc', 'Nh')
	idxTc, idxNc, idxNh = tuple(field_names.index(x) for x in relevant_field_names)
	relevant_field_indices = (idxTc, idxNc, idxNh)
	set_labels = ['+46', '+47']
	if mask is None:
		mask = get_pacs_mask()
	else:
		mask &= get_pacs_mask()
	plt.figure(figsize=(14, 9))
	axes = [plt.subplot(131 + x) for x in range(len(relevant_field_indices))]
	for i, ratio_set in enumerate(new):
		for j, idx in enumerate(relevant_field_indices):
			plt.sca(axes[j])
			ratio = ratio_set[idx]
			histxy, stats = gen_hist_and_stats(ratio, mask, x_lim=lims[idx])
			plt.plot((histxy[0] - 1)*100, histxy[1], '-',
				label=set_labels[i]+"/45")
			print("i %d, j %d" % (i, j))
			if i == 1:
				name = field_names[idx]
				plt.title(name)
				plt.xlabel("Percent change")
				plt.ylabel("Histogram count")
				if j == 0:
					plt.legend()
	show_plot()


def mask_454647(i, l):
	original_fn = fns_45
	new_fn = fns_47
	ratio_set = get_ratios_to_original(new_fn, original=original_fn)
	relevant_field_names = ('Tc', 'Nc', 'Nh')
	idxTc, idxNc, idxNh = tuple(field_names.index(x) for x in relevant_field_names)
	relevant_field_indices = (idxTc, idxNc, idxNh)
	pacs_mask = get_pacs_mask()
	ratio = ratio_set[i]
	return pacs_mask & (ratio > l[0]) & (ratio < l[1])

def mask_img(i, l):
	original_fn = fns_46
	img = load_relevant_frames(*fns_46)[i]
	return get_pacs_mask() & (img > l[0]) & (img < l[1])


def get_junk_mask():
	# Cold temp is high, Cold column is negative
	return mask_img(2, (12, np.inf)) | mask_img(4, (-np.inf, 0))


def quickrun_image_454647(i, mask, l=None, label=""):
	if l is None:
		if i == 0 or i == 2:
			l = (2, 17)
		else:
			l = (19.5, 23)
	img = load_relevant_frames(*fns_46)[i]
	img[~mask] = np.nan
	if i != 0 and i != 2:
		img = np.log10(img)
	plt.figure(figsize=(12, 12))
	plt.subplot(111)
	plt.imshow(img, origin='lower', vmin=l[0], vmax=l[1])
	plt.colorbar()
	plt.title(field_names[i]+label)
	plt.tight_layout()
	show_plot()


lims_changedT = [
	(1, 1.01),
	(.97, 1),
	(.96, 1.04),
	(.99, 1.01),
	(.2, 1.8),
	(1.035, 1.1),
]

def quickrun_histograms_changedT():
	pacs_mask = get_pacs_mask() & get_core_mask()
	plt.figure()
	axes = [plt.subplot(231 + x) for x in range(6)]
	news, originals = (fns_45_15p5, fns_47_15p5), (fns_45, fns_47)
	text = ('+45', '+47')
	for j in range(2):
		ratio_set = get_ratios_to_original(news[j], original=originals[j])
		img_set = load_relevant_frames(*originals[j])
		for i in range(len(ratio_set)):
			ratio = ratio_set[i]
			img = img_set[i]
			mask = ~np.isnan(ratio) & pacs_mask
			mask &= ~((img < 0) | (ratio < 0))
			histxy, stats = gen_hist_and_stats(ratio, mask, x_lim=lims_changedT[i])
			plt.sca(axes[i])
			plt.plot((histxy[0] - 1)*100, histxy[1], '-',
				label="Th=15.5/Th=15.7 ({})".format(text[j]))
			if j == 1:
				name = field_names[i]
				plt.title(name)
				plt.xlabel("Percent change")
				plt.ylabel("Histogram count")
				if i == 3 :
					plt.legend()
	show_plot()


def quickrun_image(i, l=None):
	if l is None:
		l = lims[i]
	ratios = (get_ratios_to_original(x)[i] for x in (fns_46, fns_47))
	for j in range(2):
		plt.subplot(121 + j)
		plt.imshow(next(ratios), origin='lower', vmin=l[0], vmax=l[1])
		plt.colorbar()
		plt.title("{0}, +{1:d}/43".format(field_names[i], 46+j))
	show_plot()


def quickrun_image_changedT(i, l=None, five=False):
	if five:
		originals = fns_45
		news = fns_45_15p5
		text = "+45"
	else:
		originals = fns_47
		news = fns_47_15p5
		text = "+47"
	if l is None:
		l = lims_changedT[i]
	img = load_relevant_frames(*originals)[i]
	ratio = get_ratios_to_original(news, original=originals)[i]
	ratio[(img < 0) | (ratio < 0)] = np.nan
	plt.subplot(121)
	plt.imshow(ratio, origin='lower', vmin=l[0], vmax=l[1])
	plt.title("{0}, Th=15.5/Th=15.7".format(field_names[i]))
	plt.colorbar()
	plt.subplot(122)
	if "N" in field_names[i]:
		img = np.log10(np.abs(img))
		text = "log10({}), Th=15.7 ({})".format(field_names[i], text)
	else:
		text = "{}, Th=15.7 ({})".format(field_names[i], text)
	plt.imshow(img, origin='lower')
	plt.title(text)
	plt.colorbar()
	show_plot()


def quickrun_negativeN(n, news, originals):
	imgs = load_relevant_frames(*originals)
	ratios = tuple(x/y for x, y in zip(load_relevant_frames(*news), imgs))
	pacs_mask = get_pacs_mask()
	for i, img in enumerate(imgs):
		if n != i:
			continue
		ratio = ratios[i]
		mask = ~np.isnan(img) & pacs_mask
		pos_mask = ~((img < 0) | (ratio < 0))
		pos_ratio = np.copy(ratio)
		pos_ratio[~mask | ~pos_mask] = np.nan
		pos_img = np.copy(img)
		pos_img[~mask | ~pos_mask] = np.nan
		plt.subplot(221)
		plt.imshow(pos_ratio, origin='lower')
		plt.title("positive ratio, {}".format(field_names[i]))
		plt.subplot(222)
		plt.imshow(np.log10(pos_img), origin='lower')
		plt.title("positive image")
		neg_ratio = np.copy(ratio)
		neg_ratio[~mask | pos_mask] = np.nan
		neg_img = np.copy(img)
		neg_img[~mask | pos_mask] = np.nan
		plt.subplot(223)
		plt.imshow(neg_ratio, origin='lower')
		plt.title("negative ratio")
		plt.subplot(224)
		plt.imshow(np.log10(-neg_img), origin='lower')
		plt.title("negative image")
	show_plot()


def quickrun_N_changedT():
	imgs47_15p7 = load_3param_TcThNcNh(dir_per1 + soln_3p_plus047)
	imgs47_15p5 = load_3param_TcThNcNh(dir_per1 + soln_3p_plus047_15p5)
	relevant_field_names = ('Tc', 'Nc', 'Nh')
	idxTc, idxNc, idxNh = tuple(field_names.index(x) - 2 for x in relevant_field_names)
	relevant_field_indices = (idxTc, idxNc, idxNh)
	ratio_Nc = imgs47_15p5[idxNc] / imgs47_15p7[idxNc]
	plt.figure()
	Nlim = (20, 23)
	Tlim = (2, 17)
	ratio_lims = (1e-1, 1e2)
	mask = get_cold_mask()
	def trim_flat(img):
		img_filtered = img[mask]
		return img_filtered.flatten()
	ratio_Nc = trim_flat(ratio_Nc)
	core_mask = trim_flat(get_core_mask())
	env_mask = trim_flat(get_envelope_mask())
	mystery_mask = trim_flat(get_mystery_envelope_mask())
	for i, idx in enumerate(relevant_field_indices):
		plt.subplot(131 + i)
		current_img = trim_flat(imgs47_15p7[idx])
		if i > 0:
			current_img = np.log10(current_img)
		ms = 10
		# plt.plot(current_img, ratio_Nc, '.', markersize=ms, alpha=0.01, label='all')
		plt.plot(current_img[mystery_mask], ratio_Nc[mystery_mask], '.', markersize=ms, alpha=0.01, label='mystery envelope')
		plt.plot(current_img[env_mask], ratio_Nc[env_mask], '.', markersize=ms, alpha=0.01, label='envelope')
		plt.plot(current_img[core_mask], ratio_Nc[core_mask], '.', markersize=ms, alpha=0.01, label='core')
		plt.ylim(ratio_lims)
		plt.yscale('log')
		if i > 0:
			plt.xlim(Nlim)
		else:
			plt.xlim(Tlim)
		plt.xlabel(relevant_field_names[i])
		plt.ylabel("Nc T=15.5/T=15.7 percent change")
		plt.legend()
	show_plot()


def quickrun_image_masked(i, mask, l=None, label=""):
	if l is None:
		if i == 0 or i == 2:
			l = (2, 17)
		else:
			l = (20, 23)
	ratio = load_relevant_frames(*fns_47)[i]
	ratio[~mask] = np.nan
	if i != 0 and i != 2:
		ratio = np.log10(ratio)
	plt.figure(figsize=(12, 12))
	plt.subplot(111)
	plt.imshow(ratio, origin='lower', vmin=l[0], vmax=l[1])
	plt.colorbar()
	plt.title(field_names[i]+label)
	plt.tight_layout()
	show_plot()


def get_better_envelope():
	return ~get_junk_mask() & mask_454647(4, (0, 1))

def get_better_core():
	return ~get_junk_mask() & mask_454647(4, (1, np.inf))

def get_notjunk_mask():
	return mask_img_full(3, (0, np.inf)) & mask_img_full(7, (0, np.inf))

def get_filament_mask():
	return mask_img_full(2, (0, 2))

def get_envelope_mask_again():
	return mask_img_full(4, (10**22., np.inf))

# Now for the "full image" runs (+46 only, 15.7 and 15.8)
def quickrun_histograms_changedT_full(ratio=True):
	plt.figure(figsize=(11, 8.5))
	pacs_mask = get_pacs_mask() & get_notjunk_mask()
	extra_masks = {"": pacs_mask,
		"Mostly clean": mask_img_full(1, (6, np.inf))
			& mask_img_full(2, (0, 5)),
	}
	frame_ids = (1, 2, 3, 4, 7, 8, 9)
	frame_names = ("Tc", "dTc", "Nc", "dNc", "Nh", "dNh", "Chisq")
	axes = [plt.subplot(331 + x) for x in range(len(frame_ids))]
	originals, news = soln_3p_plus046f, soln_3p_plus046f_15p8
	originals = load_specific_frames(originals, frames=frame_ids)
	news = load_specific_frames(news, frames=frame_ids)
	if ratio:
		lims = (
			(0.96, 1.04),
			(0.5, 1.5),
			(.2, 1.8),
			(.2, 1.8),
			(0.95, 1.0),
			(0.2, 1.8),
			(.5, 1.4)
		)
	else:
		lims = (
			(3, 17),
			(0, 15),
			(20, 23),
			(19, 23),
			(20, 23),
			(19, 23),
			(0, 3)
		)
	for i in range(len(frame_ids)):
		img = originals[i]
		if ratio:
			value = news[i] /img
		else:
			value = img
			if frame_ids[i] >= 3 and frame_ids[i] <= 8:
				value = np.log10(value)
		plt.sca(axes[i])
		mask = ~np.isnan(value) & pacs_mask
		mask &= ~((img <= 0) | (value < 0))
		for m in extra_masks:
			histxy, stats = gen_hist_and_stats(value,
				(extra_masks[m] & mask),
				x_lim=lims[i])
			if ratio:
				x_value = (histxy[0] - 1)*100
			else:
				x_value = histxy[0]
			plt.plot(x_value, histxy[1], '-',
				label="M:{}".format(m))

		name = frame_names[i]
		plt.title(name)
		if ratio:
			plt.xlabel("Percent change")
		else:
			plt.xlabel("Value")
		plt.ylabel("Histogram count")
		plt.legend()
	plt.tight_layout()
	show_plot()


def mask_changedT_full(i, l):
	ratio = load_specific_frame(soln_3p_plus046f_15p8, i) / load_specific_frame(soln_3p_plus046f, i)
	return get_pacs_mask() & (ratio < l[1]) & (ratio > l[0])

def quickrun_image_masked_full(i, mask, l=None, label=""):
	if l is None:
		if i == 1:
			l = (2, 17)
		elif i <= 3 & i >= 8:
			l = (20, 23)
		else:
			l = (0, 5)
	img = load_specific_frame(soln_3p_plus046f, i)
	img[~mask] = np.nan
	if (i >= 3) & (i <= 8):
		img = np.log10(img)
	plt.figure(figsize=(12, 12))
	plt.subplot(111)
	plt.imshow(img, origin='lower', vmin=l[0], vmax=l[1])
	plt.colorbar()
	#plt.title(field_names[i]+label)
	plt.tight_layout()
	show_plot()

def quickrun_2pimage_masked_full(i, mask, l=None, label=""):
	if l is None:
		if i == 1:
			l = (2, 17)
		elif i <= 3 & i >= 8:
			l = (20, 23)
		else:
			l = (0, 5)
	img = load_specific_frame(soln_2p_plus046, i)
	img[~mask] = np.nan
	if (i >= 3) & (i <= 4):
		img = np.log10(img)
	plt.figure(figsize=(12, 12))
	plt.subplot(111)
	plt.imshow(img, origin='lower', vmin=l[0], vmax=l[1])
	plt.colorbar()
	#plt.title(field_names[i]+label)
	plt.tight_layout()
	show_plot()

def mask_img_full(i, l):
	original_fn = soln_3p_plus046f
	img = load_specific_frame(original_fn, i)
	return get_pacs_mask() & (img > l[0]) & (img < l[1])


def quickrun_2p_chi_squared_histogram():
	img = load_specific_frame(soln_2p_plus046, 5)
	mask1 = get_pacs_mask()
	mask2 = get_notjunk_mask() & mask1
	mask3 = get_filament_mask() & mask2
	masks = [mask1, mask2, mask3]
	mask_labels = ["no mask", "possibly cold", "filament"]
	plt.figure(figsize=(11, 8.5))
	plt.subplot(111)
	for m, ml in zip(masks, mask_labels):
		histxy, stats = gen_hist_and_stats(img, m, x_lim=(0, 20),
			log=True)
		plt.plot(histxy[0], histxy[1], '-',
			label="{}".format(ml))
	plt.legend()
	plt.xlabel("X^2")
	plt.ylabel("Histogram count")
	plt.title("2parameter chi squared; masks overlap")
	show_plot()

def get_2pchisq_mask():
	img = load_specific_frame(soln_2p_plus046, 5)
	return (img > 1)

def mask_2p_img(i, l):
	original_fn = soln_2p_plus046
	img = load_specific_frame(original_fn, i)
	return get_pacs_mask() & (img > l[0]) & (img < l[1])


if __name__ == "__main__":
	# quickrun_histograms_changedT_full(ratio=False)

	# mask = mask_img_full(1, (6, np.inf)) & get_notjunk_mask()

	# mask = get_filament_mask() & mask_img_full(1, (6, np.inf))
	# mask = get_notjunk_mask() & mask_img_full(2, (0, 5))
	# mask = mask_2p_img(5, (0, 0.5))
	# mask = get_notjunk_mask() & get_filament_mask() #& ~get_envelope_mask_again()
	mask = get_2pchisq_mask()
	# mask = get_pacs_mask() | ~get_pacs_mask()
	quickrun_image_masked_full(1, mask, l=(3, 15))

	# quickrun_2p_chi_squared_histogram()

	# quickrun_image_454647(2, ~get_junk_mask())
	# next task:
	# Tc mask: -3 -> -1.2 %, 0 -> 0.42 %, .42 -> 1 %
	# Nh mask: 4->4.5, 5.5->6.5
	# Nc mask: -40->-10, -4->4, 10->30 ## THE 20% DIST WAS NEGATIVE STUFF!!!
