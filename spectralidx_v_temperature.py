import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from reanalyze_manticore import cleaner_core_mask, cleaner_junk_mask, masking_attempt
from boolean_islands import get_mask, fill_inwards
from manticore_results import mask_img_full, get_pacs_mask
from astropy.stats import median_absolute_deviation as mad
import sys

per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
power_stub = "pow-1000-0.1-"
power_run_stub = "T4-absdiff-Per1J-plus045-"
power_run_3p_stub = "T4-absdiff-Per1J-3param-plus045-"
single_comp_dir = "single_comp_beta_grid/"
two_comp_dir = "two_comp_beta_grid/"
fits_stub = ".fits"

GNILC_T_fn = per1_dir + "dustModel_Per1_SPIRE500umgrid_Temperature.fits"

with open(per1_dir+"pow_todo.txt", 'r') as pow_file:
	powers = pow_file.readlines()
powers = [x.strip() for x in powers]

m_junk_fn = "{}per1_junk_mask{}".format(per1_dir, fits_stub)
m_fila_fn = "{}per1_fila_mask{}".format(per1_dir, fits_stub)

BINS = 64

def histogram(x, x_lim=None):
	if x_lim is None:
		x_lim = (np.min(x), np.max(x))
	dhist, dedges = np.histogram(x, bins=BINS, range=x_lim)
	prep_arr = lambda a, b: np.array([a, b]).T.flatten()
	histx, histy = prep_arr(dedges[:-1], dedges[1:]), prep_arr(dhist, dhist)
	bin_centers = (dedges[:-1]+dedges[1:])/2
	return histx, histy


def gen_power_fn(beta):
	return "{:s}{:s}{:s}{:s}{:s}{:s}".format(
		per1_dir, single_comp_dir, power_run_stub,
		power_stub, beta, fits_stub
	)


def gen_3p_power_fn(beta_c, beta_h):
	return "{:s}{:s}{:s}c{:s}{:s}h{:s}{:s}{:s}".format(
		per1_dir, two_comp_dir, power_run_3p_stub,
		power_stub, beta_c,
		power_stub, beta_h,
		fits_stub
	)


def Planck_T_agreement_plot():
	planck_T = fits.getdata(GNILC_T_fn)
	m_thin = fits.getdata(m_junk_fn).astype(bool)
	m_fila = mask_img_full(1, (0, 13.5), filename_override=gen_power_fn("1.50"))
	plt.figure(figsize=(14, 9))
	ax_hist = plt.subplot(121)
	ax_curve = plt.subplot(122)
	for i, beta in enumerate(powers):
		d = fits.getdata(gen_power_fn(beta))
		diff = planck_T - d
		diff_flat = diff[m_thin & ~np.isnan(diff)].flatten()
		if i % 3 == 0:
			ax_hist.plot(*histogram(diff_flat, x_lim=(-7, 7)),
				      label=r"$\beta=${}".format(beta))
		ax_curve.plot([float(beta)], [np.median(diff_flat)],
			      '+', color='k')
		ax_curve.plot([float(beta)], [np.mean(diff_flat)],
			      'o', color='k')
		diff_flat = diff[m_fila & ~np.isnan(diff)].flatten()
		ax_curve.plot([float(beta)], [np.median(diff_flat)],
			      '+', color='b')
		ax_curve.plot([float(beta)], [np.mean(diff_flat)],
			      'o', color='b')
	ax_hist.legend()
	ax_hist.set_xlabel("Planck T minus our T (K)")
	ax_hist.set_ylabel("Histogram count")
	ax_hist.set_title(r"Distributions of T deviation from Planck GNILC with varying $\beta$")
	ax_curve.set_xlabel(r"$\beta$")
	ax_curve.set_ylabel(r"Mean ($\bullet$) or median (+) T deviation (K)")
	t1 = r"Mean/Median T deviation with varying $\beta$"
	ax_curve.set_title(t1+"\n Black: Off-filament; Blue: filament")
	plt.tight_layout()
	plt.show()


def fitted_beta_plot():
	m_thin = fits.getdata(m_junk_fn).astype(bool)
	m_fila = mask_img_full(1, (0, 13.5), filename_override=gen_power_fn("1.50"))
	plt.figure(figsize=(7, 4.5))
	ax_curve = plt.subplot(111)

	def get_mean_median(beta, mask=m_thin):
		# d is X^2
		d = fits.getdata(gen_power_fn(beta), 5)
		d_flat = d[mask & ~np.isnan(d)].flatten()
		mean, median = np.mean(d_flat), np.median(d_flat)
		return float(beta), mean, median

	betas, means, medians = zip(*[get_mean_median(beta) for beta in powers])
	betas, means, medians = np.array(betas), np.array(means), np.array(medians)
	meanfit = np.polyfit(betas, means, deg=4)
	medianfit = np.polyfit(betas, medians, deg=4)
	def polynomial(x, fit):
		deg = len(fit) - 1
		arr = [f*(x**(deg-i)) for i, f in enumerate(fit)]
		return np.sum(arr, axis=0)
	beta_range = np.arange(1.3, 2.4, 0.02)
	mean_curve = polynomial(beta_range, meanfit)
	median_curve = polynomial(beta_range, medianfit)

	mean_soln = beta_range[mean_curve == np.min(mean_curve)]
	median_soln = beta_range[median_curve == np.min(median_curve)]
	print("MEAN-minimized beta: ", mean_soln)
	print("MEDIAN-minimized beta: ", median_soln)
	plt.plot(betas, means, 'o', color='k')
	plt.plot(betas, medians, '+', color='k')
	plt.plot(beta_range, mean_curve, '-', color='k', linewidth=1)
	plt.plot(beta_range, median_curve, '--', color='k', linewidth=1)
	plt.title(r"Mean/median $\chi^2$ in thin regions varying with $\beta$")
	plt.xlabel(r"$\beta$")
	plt.ylabel(r"Mean ($\bullet$) or median (+) $\chi^2$")
	plt.show()


def fitted_beta_COLD_plot():
	outer_mask, inner_mask, fila_mask = negative_Nh_masks()
	T = fits.getdata(gen_3p_power_fn("2.05", "1.80"), 1)
	nanmask = negative_Nc_mask()
	m_thin = outer_mask & nanmask
	m_cold = inner_mask & nanmask & ~fila_mask
	m_fila = fila_mask & nanmask & (T<12) & ~inner_mask
	plt.figure(1, figsize=(14, 9))
	ax_Xscurve = plt.subplot(221)
	ax_dTcurve = plt.subplot(223)
	selected_mask = m_fila

	def repeat_everything(frame, frame_label, ax, frame_function=None):
		plt.sca(ax)
		def get_mean_median(beta, mask=selected_mask):
			with fits.open(gen_3p_power_fn(beta, "1.80")) as hdul:
				d = hdul[frame].data
				Tc = hdul[1].data
				Nc = hdul[3].data
				Nh  = hdul[7].data
			if frame_function is not None:
				d = frame_function(d)
			d_flat = d[mask & ~np.isnan(d) & (Nh > 0) & (Nc > 0) & (Tc > 6) & (d > 0)].flatten()
			mean, median = np.mean(d_flat), np.median(d_flat)
			print(mean, median)
			n = d_flat.size
			return float(beta), mean, median, n

		betas, means, medians, ns = zip(*[get_mean_median(beta) for beta in powers if float(beta) <= 2.25])
		betas, means, medians, ns = np.array(betas), np.array(means), np.array(medians), np.array(ns)
		# meanfit = np.polyfit(betas, means, deg=4)
		# medianfit = np.polyfit(betas, medians, deg=4)
		# def polynomial(x, fit):
		# 	deg = len(fit) - 1
		# 	arr = [f*(x**(deg-i)) for i, f in enumerate(fit)]
		# 	return np.sum(arr, axis=0)
		# beta_range = np.arange(1.3, 2.4, 0.02)
		# mean_curve = polynomial(beta_range, meanfit)
		# median_curve = polynomial(beta_range, medianfit)

		# mean_soln = beta_range[mean_curve == np.min(mean_curve)]
		# median_soln = beta_range[median_curve == np.min(median_curve)]
		# print("MEAN-minimized beta: ", mean_soln)
		# print("MEDIAN-minimized beta: ", median_soln)
		plt.plot(betas, means, 'o', color='k')
		plt.plot(betas, medians, '+', color='k')
		# plt.plot(beta_range, mean_curve, '-', color='k', linewidth=1)
		# plt.plot(beta_range, median_curve, '--', color='k', linewidth=1)
		plt.title(r"Mean/median {} in filamentary regions varying with $\beta$".format(frame_label))
		plt.xlabel(r"$\beta$")
		plt.ylabel(r"Mean ($\bullet$) or median (+) {}".format(frame_label))
		# plt.yscale('log')
		return betas, ns

	betas, ns_Xs = repeat_everything(9, r"$\chi^2$", ax_Xscurve)
	betas, ns_dT = repeat_everything(2, "dT", ax_dTcurve)

	ax_count = plt.subplot(222)
	plt.plot(betas, ns_Xs, '-', label="Xs")
	plt.plot(betas, ns_dT, '--', label="dT")
	plt.title("number of included pixels for stats \n in left panel (Nh>0 requirement)")
	plt.xlabel(r"$\beta$")
	plt.ylabel("Number of pixels")
	plt.legend()

	ax_mask = plt.subplot(224)
	T = fits.getdata(gen_3p_power_fn("2.05", "1.80"), 1)
	T[~selected_mask] = np.nan
	plt.imshow(T, origin='lower', vmin=9, vmax=12)
	plt.title("Mask \n Shown: T_cold; beta_c = 2.05, beta_h = 1.8")

	plt.tight_layout()
	plt.show()


def main_beta_plot():
	plt.figure(figsize=(14, 9))

	d, h = fits.getdata(gen_power_fn("2.50"), header=1)
	m_thin = fits.getdata(m_junk_fn).astype(bool) # has nan areas
	m_fila = mask_img_full(1, (0, 13.5), filename_override=gen_power_fn("1.50"))
	# m_fila = masking_attempt(max_dT=1.25, max_dN=5e21)
	# m_fila = fill_inwards(m_fila, ~np.isnan(d))
	# for i in range(3):
	# 	m_fila = get_mask(m_fila, n=6, min_size=15, dilation=0)
	m_fila = get_mask(m_fila, n=1, min_size=2, dilation=1)

	# d[~m_fila] = np.nan
	# plt.imshow(d, origin='lower', vmin=10, vmax=14)
	# plt.show()
	# sys.exit()
	# m_fila = fits.getdata(m_fila_fn).astype(bool)

	def mask_T_and_plot(data, mask, ax, x_lim=(8, 23), **kwargs):
		d_flat = data[mask & ~np.isnan(data)].flatten()
		plt.sca(ax)
		plt.plot(*histogram(d_flat, x_lim=x_lim), **kwargs)


	def typical_Xs_plot(data, mask, ax, beta, **kwargs):
		avg, med = np.nanmean(data[mask]), np.nanmedian(data[mask])
		plt.sca(ax)
		plt.plot([beta], [avg], 'o', **kwargs)
		plt.plot([beta], [med], '+', **kwargs)


	plt.figure(1)
	ax_fila = plt.subplot2grid((2, 5), (0, 0), colspan=3)
	ax_thin = plt.subplot2grid((2, 5), (1, 0), colspan=3)
	ax_m_fila = plt.subplot2grid((2, 5), (0, 3), colspan=2)
	ax_m_thin = plt.subplot2grid((2, 5), (1, 3), colspan=2)
	plt.figure(2)
	ax_Xs_hist = plt.subplot(111)
	for beta in powers:
		d, h = fits.getdata(gen_power_fn(beta), header=1)
		mask_T_and_plot(d, mask=m_fila, ax=ax_fila,
				label=r"$\beta=${}".format(beta))
		mask_T_and_plot(d, mask=m_thin, ax=ax_thin,
				label=r"$\beta=${}".format(beta))
		d, h = fits.getdata(gen_power_fn(beta), 5, header=1)
		typical_Xs_plot(d, m_thin, ax_Xs_hist, float(beta), color='k')
		# typical_Xs_plot(d, m_fila, ax_Xs_hist, float(beta), color='b')
		#	mask_T_and_plot(d, mask=m_fila, ax=ax_Xs_hist, x_lim=(0, 5),
#			label=r"$\beta=${}".format(beta))
	ax_fila.legend(), ax_thin.legend()
	for ax, txt in zip((ax_fila, ax_thin), ("Filament histograms", "Off-filament histograms")):
		ax.set_xlabel("T (K)")
		ax.set_ylabel("Histogram count")
		ax.set_title(txt)

	d, h = fits.getdata(gen_power_fn("1.80"), header=1)

	dc = d.copy()
	dc[~m_fila] = np.nan
	plt.sca(ax_m_fila)
	plt.imshow(dc, origin='lower', vmin=10, vmax=15)
	plt.colorbar()
	plt.title("Filament mask")

	dc = d.copy()
	dc[~m_thin] = np.nan
	plt.sca(ax_m_thin)
	plt.imshow(dc, origin='lower', vmin=14, vmax=17)
	plt.colorbar()
	plt.title(r"Off-filament mask (Pictured: $\beta=1.8$)")
	# plt.sca(ax_Xs_hist)
	# plt.yscale('log')
	# plt.legend()

	plt.tight_layout()
	plt.show()


def negative_Nh_masks():
	m_fila = masking_attempt()
	pacs_mask = get_pacs_mask()
	# Xs_stack = []
	Nh_lowbeta_mask_stack = []
	Nh_highbeta_mask_stack = []
	for beta in powers:
		if float(beta) > 2.1:
			continue
		with fits.open(gen_3p_power_fn(beta, "1.80")) as hdul:
			# Tc = hdul[1].data
			Nh = hdul[7].data
			# Xs = hdul[9].data
			if float(beta) >= 1.8:
				Nh_lowbeta_mask_stack.append((Nh < 0).astype(int))
			elif float(beta) < 1.8:
				Nh_highbeta_mask_stack.append((Nh < 0).astype(int))
	Nh_lowbeta_cumulative = np.sum(Nh_lowbeta_mask_stack, axis=0)
	Nh_highbeta_cumulative = np.sum(Nh_highbeta_mask_stack, axis=0)
	Th_fixed = np.nanmedian(fits.getdata(gen_3p_power_fn("1.80", "1.80"), 5))
	Th = fits.getdata(gen_power_fn("1.80"), 1)
	outer_mask = (Nh_lowbeta_cumulative < 1) & pacs_mask
	Th_outer = Th[outer_mask].flatten()
	inner_mask = (Nh_lowbeta_cumulative >= 1) & pacs_mask
	Th_inner = Th[inner_mask].flatten()
	fila_mask = (Nh_highbeta_cumulative >= 1) & pacs_mask
	Th_fila = Th[fila_mask].flatten()

	# plt.figure(figsize=(14, 9))


	#### PLOT MASKS, SEE WHERE GAS IS
	"""
	plt.subplot(221)
	Thc = Th.copy()
	Thc[~outer_mask] = np.nan
	plt.imshow(Thc, origin='lower', vmin=13, vmax=18)
	plt.title("Outer gas")
	plt.subplot(222)
	Thc = Th.copy()
	Thc[~inner_mask] = np.nan
	plt.imshow(Thc, origin='lower', vmin=13, vmax=18)
	plt.title("Inner gas")
	plt.subplot(223)
	Thc = Th.copy()
	Thc[~fila_mask] = np.nan
	plt.imshow(Thc, origin='lower', vmin=13, vmax=18)
	plt.title("Filament gas")
	plt.tight_layout()
	plt.show()
	return
	"""

	#### PLOT HISTOGRAMS, SEE HOW GAS FIT IN SINGLE-T
	"""
	plt.plot(*histogram(Th_outer, x_lim=(13, 18)), label="Outer Gas")
	plt.plot(*histogram(Th_inner, x_lim=(13, 18)), label="Inner Gas")
	plt.plot(*histogram(Th_fila, x_lim=(13, 18)), label="Filament Gas")
	plt.plot([Th_fixed, Th_fixed], plt.ylim(), '--')
	plt.legend()
	plt.show()
	return
	"""

	#### PLOT SEDS, SEE WHAT PIXEL DATA LOOK LIKE
	"""
	with fits.open(gen_power_fn("1.80")) as hdul:
		band_data = [hdul[x].data for x in [11, 13, 15, 17]]
	inner_band_data = []
	outer_band_data = []
	fila_band_data = []
	for i in range(4):
		bd = band_data[i]
		ibd = bd[inner_mask].flatten()
		notnan = ~np.isnan(ibd)
		inner_band_data.append((np.median(ibd[notnan]), mad(ibd[notnan])))
		obd = bd[outer_mask].flatten()
		notnan = ~np.isnan(obd)
		outer_band_data.append((np.median(obd[notnan]), mad(obd[notnan])))
		fbd = bd[fila_mask].flatten()
		notnan = ~np.isnan(fbd)
		fila_band_data.append((np.median(fbd[notnan]), mad(fbd[notnan])))
	bands = [160, 250, 350, 500]
	ibd_v, ibd_e = zip(*inner_band_data)
	obd_v, obd_e = zip(*outer_band_data)
	fbd_v, fbd_e = zip(*fila_band_data)
	print(ibd_v)
	print(obd_v)
	print(fbd_v)
	dx = 2
	plt.errorbar([x+dx for x in bands], obd_v, yerr=obd_e, fmt='.', capsize=3,
		label="Outer Gas")
	plt.errorbar([x-dx for x in bands], ibd_v, yerr=ibd_e, fmt='.', capsize=3,
		label="Inner Gas")
	plt.errorbar([x-(2*dx) for x in bands], fbd_v, yerr=fbd_e, fmt='.', capsize=3,
		label="Filament Gas")
	plt.legend()
	plt.show()
	return
	"""

	return outer_mask, inner_mask, fila_mask


def negative_Nc_mask():
	pacs_mask = get_pacs_mask()
	spire_mask = ~np.isnan(fits.getdata(gen_power_fn("1.80"), 1))
	nanmask = pacs_mask & spire_mask
	Nc_mask_stack = []
	for beta in powers:
		if float(beta) > 2.1:
			continue
		with fits.open(gen_3p_power_fn(beta, "1.80")) as hdul:
			Nc = hdul[3].data
			Nc_mask_stack.append((Nc < 0).astype(int))
	Nc_cumulative = np.sum(Nc_mask_stack, axis=0)
	"""
	plt.imshow(Nc_cumulative, origin='lower')
	plt.show()
	return
	"""
	Nc_mask = (Nc_cumulative < 5) & (Nc_cumulative >= 0) & nanmask
	"""
	plt.imshow(Nc_mask, origin='lower')
	plt.show()
	return
	"""
	return Nc_mask


def best_beta_per_pixel():
	nanmask = negative_Nc_mask()
	Xs_stack = []
	for beta in powers:
		if float(beta) > 2.3:
			continue
		with fits.open(gen_3p_power_fn(beta, "1.80")) as hdul:
			Xs = hdul[9].data
			dT = hdul[2].data
			Nc = hdul[3].data
			Nh = hdul[7].data
			bad_mask = (dT <= 0)
			bad_mask |= (Nc < 0) | (Nh < 0)
			bad_mask |= ~nanmask
			Xs[bad_mask] = np.inf
			Xs_stack.append(Xs)
	Xsmin = np.argmin(Xs, axis=0)
	plt.imshow(Xsmin, origin='lower')
	plt.show()


if __name__ == "__main__":
	# Planck_T_agreement_plot()
	fitted_beta_COLD_plot()
