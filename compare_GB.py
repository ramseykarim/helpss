import numpy as np
import sys
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import AA_regrid_healpy as rhp

from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord, FK5
from astropy import units as u

import planck_mask
planck_mask.BINS = int(6/0.3)
from reanalyze_manticore import histogram

gb_dir = "/home/rkarim/Research/Filaments/Data/"
per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
gb_fn_gen = lambda s: "{}perseus04-{}.fits".format(gb_dir, str(int(s)))
gb_T_fn = gb_dir + "HGBS_perseus_dust_temperature_map.fits"
gb_N_fn = gb_dir + "HGBS_perseus_column_density_map.fits"
gb_TN_fn = per1_dir+"GB_TN_map_regrid.fits"

zari_fn = "/n/pleiades/lgm/Dropbox/Perseus/zari-planck_herschel_fit.fits"

band_shortcut = {
	160: "PACS160um",
	250: "SPIRE250um",
	350: "SPIRE350um",
	500: "SPIRE500um"
}
our_fn_gen = lambda s: "{}{}-image.fits".format(per1_dir, band_shortcut[s])
DL3_T_fn = per1_dir + "T4-absdiff-Per1J-plus045-DL3.fits"
OH5_T_fn = per1_dir + "T4-absdiff-Per1J-plus045.fits"
pow2_T_fn = per1_dir + "T4-absdiff-Per1J-plus045-pow-1000-0.1-2.0.fits"
pow18_T_fn = per1_dir + "T4-absdiff-Per1J-plus045-pow-1000-0.1-1.8.fits"
pow17_T_fn = per1_dir + "T4-absdiff-Per1J-plus045-pow-1000-0.1-1.7.fits"
pow16_T_fn = per1_dir + "T4-absdiff-Per1J-plus045-pow-1000-0.1-1.6-3band.fits"
pow15_T_fn = per1_dir + "T4-absdiff-Per1J-plus045-pow-1000-0.1-1.5.fits"
power_fns = {2.0: pow2_T_fn, 1.8: pow18_T_fn,
             1.7: pow17_T_fn, 1.6: pow16_T_fn,
	     1.5: pow15_T_fn}

DL3_crop_fn = per1_dir + "T4-absdiff-Per1J-plus045-DL3-crop1.fits"
DL3_CMB_crop_fn = per1_dir + "T4-absdiff-Per1J-plus045-DL3-CMB-crop1.fits"
OH5_CMB_crop_fn = per1_dir + "T4-absdiff-Per1J-plus045-OH5-CMB-crop1.fits"
OH5_crop_fn = per1_dir + "T4-absdiff-Per1J-plus045-OH5-crop1.fits"

def plot_frame(filename, frame, vmin, vmax, ax, title="",
	       f=None):
	try:
		f1, f2 = filename
		image1 = fits.getdata(f1, frame)
		image2 = fits.getdata(f2, frame)
		if f is not None:
			image1, image2 = f(image1), f(image2)
		image = image1 - image2
	except TypeError:
		image = fits.getdata(filename, frame)
		if f is not None:
			image = f(image)
	except ValueError:
		image = fits.getdata(filename, frame)
		if f is not None:
			image = f(image)
	plt.sca(ax)
	plt.imshow(image, origin='lower', vmin=vmin, vmax=vmax)
	plt.title(title)
	plt.colorbar()


def get_zari():
	with fits.open(zari_fn) as hdul:
		data = hdul[0].data
		head = hdul[0].header
		T = data[2, :, :]
		beta = data[4, :, :]
		Xs = data[6, :, :]
	return beta, head

# Examine Sadavoy paper B1-E region (and Zari 2016)
sadavoy_coord = SkyCoord("3:36:20.00 +31:12:00.00", frame=FK5, unit=(u.hourangle, u.deg))
sadavoy_dimensions = (14*u.arcminute, (105./3600)*u.hourangle) # Dec, RA (reversed)

gb_T_map, head = fits.getdata(gb_TN_fn, 1, header=True)
gb_w = WCS(head)
our20_T_map, head = fits.getdata(pow2_T_fn, 1, header=True)
our20_w = WCS(head)
our16_T_map, head = fits.getdata(pow16_T_fn, 1, header=True)
our16_w = WCS(head)
zari_T_map, head = get_zari()
zari_w = WCS(head, naxis=2)

gb_cut = Cutout2D(gb_T_map, sadavoy_coord, sadavoy_dimensions, wcs=gb_w)
our20_cut = Cutout2D(our20_T_map, sadavoy_coord, sadavoy_dimensions, wcs=our20_w)
our16_cut = Cutout2D(our16_T_map, sadavoy_coord, sadavoy_dimensions, wcs=our16_w)
zari_cut = Cutout2D(zari_T_map, sadavoy_coord, sadavoy_dimensions, wcs=zari_w)
del gb_T_map, our20_T_map, our16_T_map, zari_T_map

plt.figure(figsize=(14, 9))
plt.subplot(111)
plt.imshow(zari_cut.data, origin='lower', vmin=1.59, vmax=1.62)
plt.colorbar()
plt.title(r"Zari+2016 $\beta$")
plt.show()

"""
plt.figure(figsize=(14, 9))
plt.subplot(221)
plt.imshow(gb_cut.data, origin='lower', vmin=12, vmax=20)
plt.title("Gould Belt T map")
plt.colorbar()
plt.subplot(222)
plt.imshow(our20_cut.data, origin='lower', vmin=12, vmax=20)
plt.title(r"HELPSS T map, $\beta=2$")
plt.colorbar()
plt.subplot(223)
plt.imshow(zari_cut.data, origin='lower', vmin=12, vmax=20)
plt.title(r"Zari+2016 T map")
plt.colorbar()
plt.subplot(224)
plt.imshow(our16_cut.data, origin='lower', vmin=12, vmax=20)
plt.title(r"HELPSS T map, $\beta=1.6$")
plt.colorbar()
"""
"""
plt.subplot(325)
plt.imshow(gb_cut.data - our20_cut.data, origin='lower')
plt.title("GB map minus HELPSS 2.0 map")
plt.colorbar()
plt.subplot(326)
plt.imshow(zari_cut.data - our16_cut.data, origin='lower')
plt.title("Zari map minus HELPSS 1.6 map")
plt.colorbar()
"""
"""
plt.figure(figsize=(14, 9))
plt.subplot(111)
histogram(gb_cut.data.flatten(), x_lim=(12, 20), label="GB", text=0)
histogram(our20_cut.data.flatten(), x_lim=(12, 20), label="HELPSS 2.0", text=0)
histogram(our16_cut.data.flatten(), x_lim=(12, 20), label="HELPSS 1.6", text=0)
histogram(zari_cut.data.flatten(), x_lim=(12, 20), label="Zari+2016", text=0)
plt.xlabel("T (k)")
plt.ylabel("Histogram Count")
plt.title("Temperature histograms in cutout region")
plt.legend()
"""
plt.show()









# plt.figure(figsize=(12, 10))
#axes = iter([plt.subplot(221 + i) for i in range(4)])
# comparison_fn = gb_TN_fn
#comparison_fn = OH5_T_fn
#comparison_fn = DL3_T_fn


# plot_frame((gb_TN_fn, pow17_T_fn), 1, vmin=-2, vmax=2,
# 	   ax=plt.subplot(111),
# 	   title=r'GB-1.8 T map')

"""
for beta in power_fns:
    plot_frame(power_fns[beta], 1,
	       vmin=8, vmax=20,
               ax=next(axes),
	       title=r'$\beta=${:.1f} (T map)'.format(beta))
"""
"""
for beta in power_fns:
    plot_frame((comparison_fn, power_fns[beta]), 1,
	       vmin=-3, vmax=3, #GB
#	       vmin=-5, vmax=0, #DL3
#	       vmin=-3, vmax=1.1, #OH5
               ax=next(axes),
	       title=r'GB result minus $\beta=${:.1f} (T map)'.format(beta))
"""
"""
for beta in power_fns:
    plot_frame(power_fns[beta], 3,
	       vmin=20.5, vmax=22,
               ax=next(axes),
	       title=r'$\beta=${:.1f} (N map)'.format(beta),
	       f=np.log10)
"""
"""
for beta in power_fns:
    plot_frame((comparison_fn, power_fns[beta]), 3,
	       vmin=-0.1, vmax=0.3, #GB
#	       vmin=0.6, vmax=1.0, #DL3
#	       vmin=-0.25, vmax=.1, #OH5
               ax=next(axes),
	       title=r'GB result minus $\beta=${:.1f} (T map)'.format(beta),
	       f=np.log10)
"""


"""
## CMB
plt.figure(figsize=(14, 10))
plt.subplot(121)
#plt.imshow(CMB_crop_fits - DL3_crop_fits, origin='lower', vmin=-1, vmax=3)
plt.imshow(np.log10(DL3_CMB_crop_fits), origin='lower', vmin=19, vmax=23)#, vmin=10, vmax=35)
#plt.title("beta=1.7 map minus beta=2 map, {}".format("temperature"))
plt.title("Temperature map, single component DL3.1 fit with CMB")
plt.colorbar()

plt.subplot(122)
plt.imshow(np.log10(DL3_crop_fits), origin='lower', vmin=19, vmax=23)#, vmin=10, vmax=17)
plt.title("same as left panel, but without CMB")
plt.colorbar()

plt.tight_layout()
plt.show()
# plt.savefig("/n/sgraraid/filaments/data/TEST4/regridding_stuff/DL3_CMB_hiT.png")
"""
