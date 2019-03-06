import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import sys


def convolve(a1, a2):
	ft = np.fft.fft2(a1)*np.fft.fft2(a2)
	return np.fft.ifft2(ft)


per1_dir_stub = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
spire_fn = "SPIRE350um-image.fits"
spire_fn = per1_dir_stub + spire_fn

new_planck_fn = "/n/sgraraid/filaments/data/HFI_SkyMap_857_resSPIRE350.fits"
new_planck_fn = per1_dir_stub + "SPIRE350um-image-PLANCKfabricated_TESTint_cc.fits"
spire_fn = "/n/sgraraid/filaments/data/HFI_SkyMap_857_resSPIRE350.fits"


# HERSCHEL SPIRE 350
with fits.open(spire_fn) as hdul:
    h_data = hdul[0].data
    h_header = hdul[0].header
    w = WCS(hdul[0].header)
h_data[np.isnan(h_data)] = 0

# PLANCK (on Herschel350 pixels)
with fits.open(new_planck_fn) as hdul:
    p_data = hdul[0].data
p_data[np.isnan(p_data)] = 0


# LOOKING AT HOW MUCH BRIGHTER HERSCHEL IS THAN PLANCK
# diff_data = h_data - p_data
# tolerance = 1
# zero_fraction = 100 * np.sum(np.ones(h_data.shape)[(np.abs(diff_data) <= tolerance) & (h_data != 0) & (p_data != 0)]) / diff_data.size
# print("%.3f pct of pixels consistent with zero within %.5f MJy/sr" % (zero_fraction, tolerance))
# negative_fraction = np.sum(np.ones(h_data.shape)[(diff_data < 0) & (h_data != 0) & (p_data != 0)]) / diff_data.size
# print("%.3f pct of pixels below zero" % (negative_fraction))
# plt.hist(diff_data[(h_data != 0) & (p_data != 0)].flat, log=True, range=(-10, 10), bins=128)
# plt.show()
# sys.exit()

#""" # DELETE
#CONVOLVING HERSCHEL UP TO PLANCK (calculate kernel)
dtheta_dpix1 = np.sqrt(np.sum([(x2 - x1)**2 for (x1, x2) in zip(w.wcs_pix2world(0, 0, 0), w.wcs_pix2world(0, 1, 0))]))*60
dtheta_dpix2 = np.sqrt(np.sum([(x2 - x1)**2 for (x1, x2) in zip(w.wcs_pix2world(0, 0, 0), w.wcs_pix2world(1, 0, 0))]))*60
dtheta_dpix = (dtheta_dpix1 + dtheta_dpix2)/2  # arcminutes
print("%f arcminute / %f arcsecond pixel scale" % (dtheta_dpix, dtheta_dpix*60))
x, y = np.arange(h_data.shape[1]) - h_data.shape[1]//2, np.arange(h_data.shape[0]) - h_data.shape[0]//2
x, y = x*dtheta_dpix2, y*dtheta_dpix1
xx, yy = np.meshgrid(x, y)
# h_beam, p_beam = 0.415, 4.63 # arcminutes
# h_beam, p_beam = 0.415, 5 # arcminutes
h_beam, p_beam = 5, 10 # arcminutes # this is for planck vs planck. "h_beam" refers to 857GHz data.
conv_beam_width = np.sqrt(p_beam*p_beam - h_beam*h_beam)/2.35
conv_beam = np.exp(-(xx*xx + yy*yy)/(2 * conv_beam_width*conv_beam_width))
del x, y, xx, yy
conv_beam /= np.sum(conv_beam)  # NORM TO 1

# 3D PLOT FOR BEAM
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.plot_surface(xx, yy, conv_beam)
# plt.show()

# CONVOLVE HERSCHEL
h_data_convolved = np.real(np.fft.fftshift(convolve(h_data, conv_beam)))
del conv_beam

# PLOT CONVOLVED HERSCHEL
# plt.subplot(131)
# plt.imshow(np.log10(h_data), origin='lower', vmin=0.5, vmax=3.5)
# plt.title("Hersch")

# ZERO PAD CONVOLVED HERSCHEL
h_data_convolved[h_data == 0] = 0 # ensure 0 padded regions remain 0 padded
# del h_data

# plt.subplot(132)
# plt.imshow(np.log10(h_data_convolved), origin='lower', vmin=0.5, vmax=3.5)
# plt.title("H conv")
# plt.subplot(133)
# plt.imshow(np.log10(p_data), origin='lower', vmin=0.5, vmax=3.5)
# plt.title("Planck")
# plt.show()

#""" # DELETE
# h_data_convolved = h_data

# CONVOLVED HERSCHEL MINUS PLANCK
diff_data = h_data_convolved / p_data

# LOOKING AT HOW MUCH BRIGHTER HERSCHEL (convolved) IS THAN PLANCK
# tolerance = 1
# zero_fraction = 100 * np.sum(np.ones(h_data.shape)[(np.abs(diff_data) <= tolerance) & (h_data != 0) & (p_data != 0)]) / diff_data.size
# print("%.3f pct of pixels consistent with zero within %.5f MJy/sr" % (zero_fraction, tolerance))
# negative_fraction = np.sum(np.ones(h_data.shape)[(diff_data < 0) & (h_data != 0) & (p_data != 0)]) / diff_data.size
# print("%.3f pct of pixels below zero" % (negative_fraction))
# plt.hist(diff_data[(h_data != 0) & (p_data != 0)].flat, log=True, range=(-10, 10), bins=128)
# plt.show()
# sys.exit()

mean_diff = np.nanmean(diff_data)
med_diff = np.nanmedian(diff_data)
rms_diff = np.nanstd(diff_data)
print("MEAN: %.3f, MEDIAN: %.3f, RMS: %.3f" % (mean_diff, med_diff, rms_diff))

plt.subplot(131)
plt.imshow(np.log10(h_data_convolved), origin='lower', vmin=0.5, vmax=3.5)
plt.title("log[convolved HFI857]")
cbar = plt.colorbar()
cbar.set_label("log10[MJy/sr]", rotation=270)
plt.subplot(132)
plt.imshow(diff_data, origin='lower', vmin=1., vmax=1.5)
plt.title("HFI857 / Planck_fab")
cbar = plt.colorbar()
cbar.set_label("MJy/sr", rotation=270)
plt.subplot(133)
plt.imshow(np.log10(p_data), origin='lower', vmin=0.5, vmax=3.5)
plt.title("log[Planck_fab]")
cbar = plt.colorbar()
cbar.set_label("log10[MJy/sr]", rotation=270)
plt.show()


# THIS METHOD COULD BE USED FOR ROTATION, BUT NEED TO RETRIEVE WCS INFO!! probably should do it a different way
# h_data_rotated = rotate(h_data, 30.)
# plt.subplot(121)
# plt.imshow(np.log10(h_data), origin='lower', vmin=0.5, vmax=3.5)
# plt.colorbar()
# plt.subplot(122)
# plt.imshow(np.log10(h_data_rotated), origin='lower', vmin=0.5, vmax=3.5)
# plt.show()
