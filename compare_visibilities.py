#!/astromake/opt/python/anaconda3/bin/python
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.ndimage.interpolation import rotate
import sys

np.warnings.filterwarnings('ignore')

def convolve(a1, a2):
    ft = np.fft.fft2(a1)*np.fft.fft2(a2)
    return np.fft.ifft2(ft)

band_crop_params = {
    "PACS160um": (1704, 1729, 2700),
    "SPIRE250um": (858, 858, 1625),
    "SPIRE350um": (501, 512, 1000),
    "SPIRE500um": (379, 360, 690),
}

bandpass_beam_sizes = {
    "PACS160um": np.sqrt(11.64 * 15.65)/60,
    "SPIRE250um": 18.4/60,
    "SPIRE350um": 25.2/60,
    "SPIRE500um": 36.7/60,
}


home_stub = "/n/sgraraid/filaments/data/"

# Per 1 specific filenames and stuff
dir_stub = "%sTEST4/Per/testregion1342190326JS/" % home_stub

fn_Herschel = lambda band: "%s%s-image.fits" % (dir_stub, band)
fn_Planck = lambda band: "%s%s-image-PLANCKfabricated.fits" % (dir_stub, band)
# fn_Planck = lambda x: "/n/sgraraid/filaments/data/HFI_SkyMap_857_resSPIRE350.fits"

if len(sys.argv) > 1:
	band_stub = sys.argv[1]
	field_stub = sys.argv[2]
else:
	band_stub = "SPIRE350um"
	field_stub = "Per1"


H_fn = fn_Herschel(band_stub)
P_fn = fn_Planck(band_stub)
print(H_fn)
print(P_fn)


with fits.open(H_fn) as hdul:
	H_d = hdul[0].data
	H_h = hdul[0].header
	H_w = WCS(H_h)
with fits.open(P_fn) as hdul:
	P_d = hdul[0].data
	P_h = hdul[0].header

P_d[np.isnan(P_d)] = 0
H_d[np.isnan(H_d)] = 0
H_d = H_d[1:, 1:]
P_d = P_d[1:, 1:]
print(H_d.shape, P_d.shape)

MJysr_to_Jyas2 = 4*np.pi/4.25e4
H_d *= MJysr_to_Jyas2
P_d *= MJysr_to_Jyas2

dtheta_dpix1 = np.sqrt(np.sum([(x2 - x1)**2 for (x1, x2) in zip(H_w.wcs_pix2world(0, 0, 0), H_w.wcs_pix2world(0, 1, 0))]))*60
dtheta_dpix2 = np.sqrt(np.sum([(x2 - x1)**2 for (x1, x2) in zip(H_w.wcs_pix2world(0, 0, 0), H_w.wcs_pix2world(1, 0, 0))]))*60
dtheta_dpix = (dtheta_dpix1 + dtheta_dpix2)/2  # arcminutes
print("%f arcminute / %f arcsecond pixel scale" % (dtheta_dpix, dtheta_dpix*60))

# ROTATION
# From original use:
# 	[THIS METHOD COULD BE USED FOR ROTATION, BUT NEED TO RETRIEVE WCS INFO!! probably should do it a different way]
rotation_value = 31.
P_d = rotate(P_d, rotation_value)
H_d = rotate(H_d, rotation_value)
x1lo, x2lo, box_length = band_crop_params[band_stub]
crop_quick = lambda arr: arr[x1lo:x1lo+box_length, x2lo:x2lo+box_length]
P_d = crop_quick(P_d)
H_d = crop_quick(H_d)

H_ft = np.fft.fft2(H_d)
P_ft = np.fft.fft2(P_d)

h_beam, p_beam = bandpass_beam_sizes[band_stub], 10 # arcminutes
# CHOOSE ONE OF THE FOLLOWING
deconvolve_planck = True
convolve_herschel = False
#----------------------------------------------------------------------------------------
if deconvolve_planck and not convolve_herschel:
	# DECONVOLVE PLANCK (setup)
	x, y = np.arange(H_d.shape[1]) - H_d.shape[1]//2, np.arange(H_d.shape[0]) - H_d.shape[0]//2
	x, y = x*dtheta_dpix2, y*dtheta_dpix1
	xx, yy = np.meshgrid(x, y)
	deconv_beam_width = p_beam/2.35
	deconv_beam = np.exp(-(xx*xx + yy*yy)/(2 * deconv_beam_width*deconv_beam_width))
	deconv_beam /= np.sum(deconv_beam)  # NORM TO 1
	conv_beam_width = h_beam/2.35
	conv_beam = np.exp(-(xx*xx + yy*yy)/(2 * conv_beam_width*conv_beam_width))
	conv_beam /= np.sum(conv_beam)  # NORM TO 1

	# DECONVOLVE PLANCK (actually do it)
	# NOTE: ok this sucks, so we probably need to rotate/crop the image before this happens. TODO!!!!
	ft_conv_beam = np.fft.fft2(conv_beam)
	ft_deconv_beam = np.fft.fft2(deconv_beam)
	P_ft *= ft_conv_beam / ft_deconv_beam

#----------------------------------------------------------------------------------------
if convolve_herschel and not deconvolve_planck:
	# CONVOLVE HERSCHEL UP TO PLANCK (setup)
	x, y = np.arange(H_d.shape[1]) - H_d.shape[1]//2, np.arange(H_d.shape[0]) - H_d.shape[0]//2
	x, y = x*dtheta_dpix2, y*dtheta_dpix1
	xx, yy = np.meshgrid(x, y)
	conv_beam_width = np.sqrt(p_beam*p_beam - h_beam*h_beam)/2.35
	conv_beam = np.exp(-(xx*xx + yy*yy)/(2 * conv_beam_width*conv_beam_width))
	conv_beam /= np.sum(conv_beam)  # NORM TO 1

	# CONVOLVE HERSCHEL (actually do it)
	ft_conv_beam = np.fft.fft2(conv_beam)
	H_ft *= ft_conv_beam

#----------------------------------------------------------------------------------------

H_d = np.real(np.fft.fftshift(np.fft.ifft2(H_ft)))
P_d = np.real(np.fft.ifft2(P_ft))

central_f_ratio = np.absolute(H_ft[0, 0]) / np.absolute(P_ft[0, 0])

P_ft = np.fft.fftshift(P_ft)
H_ft = np.fft.fftshift(H_ft)


# # # PLOT IMAGES
# plt.subplot(121)
# plt.imshow(H_d, origin='lower')
# plt.subplot(122)
# plt.imshow(P_d, origin='lower')
# plt.show()
# sys.exit()

# PLOT VISIBILITIES
# plt.subplot(221)
# plt.imshow(np.real(H_ft), origin='lower')
# plt.subplot(222)
# plt.imshow(np.imag(H_ft), origin='lower')
# plt.subplot(223)
# plt.imshow(np.real(P_ft), origin='lower')
# plt.subplot(224)
# plt.imshow(np.imag(P_ft), origin='lower')
# plt.show()

H_ft, P_ft = np.absolute(H_ft), np.absolute(P_ft)

print("Central frequencies ratio Herschel/Planck = %f" % central_f_ratio)

# ANGULAR SCALE GRID
x, y = np.arange(H_ft.shape[0]) - H_ft.shape[0]//2, np.arange(H_ft.shape[1]) - H_ft.shape[1]//2
x = x / (2 * (H_ft.shape[0]//2) * dtheta_dpix)
y = y / (2 * (H_ft.shape[1]//2) * dtheta_dpix)
xx, yy = np.meshgrid(x, y, indexing='ij')
radial_distance_grid = np.sqrt(xx*xx + yy*yy)  # units: inverse arcminutes
ang_scale = 1/radial_distance_grid  # Yeah this is better now, matches up with the paper


# plt.plot(ang_scale[ang_scale <= 0.5].flat, P_ft[ang_scale <= 0.5].flat, 'b,', alpha=0.01)
# plt.plot(ang_scale[ang_scale <= 0.5].flat, H_ft[ang_scale <= 0.5].flat, 'r,', alpha=0.01)

# plt.plot(ang_scale[ang_scale > 0.5].flat, P_ft[ang_scale > 0.5].flat, 'b^', alpha=0.1)
# plt.plot(ang_scale[ang_scale > 0.5].flat, H_ft[ang_scale > 0.5].flat, 'r^', alpha=0.1)

cross_cal_mask = (ang_scale < 100) & (ang_scale > 10)
vis_ratio = H_ft[cross_cal_mask] / P_ft[cross_cal_mask]
cross_cal_factor, cross_cal_error = np.mean(vis_ratio), np.std(vis_ratio)
del vis_ratio
print("Cross-calibration factor Herschel/Planck = %.3f +/- %.3f" % (cross_cal_factor, cross_cal_error))

plt.figure(figsize=(15, 15))
plt.subplot(111)
plt.plot(ang_scale.flat, H_ft.flat, 'r,', alpha=1.)
ang_scale_min = 3  # arcminutes (the scale the paper uses)
plt.plot(ang_scale[ang_scale > ang_scale_min].flat, P_ft[ang_scale > ang_scale_min].flat, 'b,', alpha=1.)


BINS = 128
P_hist, Pxedges, Pyedges = np.histogram2d(np.log10(ang_scale.flat), np.log10(P_ft.flat),
	bins=BINS, range=[[-2, 1], [-3, 4]])
H_hist, Hxedges, Hyedges = np.histogram2d(np.log10(ang_scale.flat), np.log10(H_ft.flat),
	bins=BINS, range=[[-2, 1], [-3, 4]])

assert np.all(Pxedges == Hxedges)
assert np.all(Pyedges == Hyedges)

xbin_centers = (Pxedges[1:] + Pxedges[:-1])/2
ybin_centers = (Pyedges[1:] + Pyedges[:-1])/2

# rgb = np.stack([H_hist, np.zeros(H_hist.shape), P_hist])
# rgb = np.swapaxes(rgb, 0, 2)

# plt.subplot(121)
# plt.imshow((P_hist.transpose()), origin='lower', cmap='inferno')
# plt.gca().invert_xaxis()
# plt.subplot(122)
# plt.imshow((H_hist.transpose()), origin='lower', cmap='inferno')
# plt.gca().invert_xaxis()
# plt.show()
plt.title("Herschel %s / Planck C-C factor: %.2f +/- %.2f" % (sys.argv[1], cross_cal_factor, cross_cal_error))

plt.yscale('log')
plt.xscale('log')
plt.xlim([1e2, 1e-1])
plt.ylim([1e-4, 1e4])

# plt.gca().invert_xaxis()
plt.ylabel("flux (Jy/as$^2$?)")
plt.xlabel("angular scale (arcminutes)")
plt.legend(handles=[mpatches.Patch(color='red', label='Herschel'), mpatches.Patch(color='blue', label='Planck')])
plt.tight_layout()
plt.show()

