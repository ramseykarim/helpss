import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.io import fits
import healpy as hp
import AA_regrid_healpy as rhp
import AA_planck_obs as pobs
import planck_mask as pm
from astropy.wcs import WCS


per1_dir_stub = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
planck_dir = "/n/sgraraid/filaments/data/"
# The filename below is correct, even after the directory move
dust_directory = "/n/sgraraid/filaments/data/TEST4/regridding_stuff/"
spire_fn = "SPIRE350um-image-remapped.fits"
spire_fn = per1_dir_stub + spire_fn
planck_fns = {
    "F100": 3,
    "F143": 4,
    "F217": 5,
    "F353": "HFI_SkyMap_353-psb-field-IQU_2048_R3.00_full.fits", # IN Kcmb UNITS????
    "F545": "HFI_SkyMap_545-field-Int_2048_R3.00_full.fits",
    "F857": "HFI_SkyMap_857-field-Int_2048_R3.00_full.fits", #"HFI_SkyMap_857_2048_R2.02_full.fits",
}
field_stub = "Per1"

planck_unit_conversions = {
    "F100": 244.1,
    "F143": 371.74,
    "F217": 483.690,
    "F353": 287.450,
    "F545": 1,
    "F857": 1,
}


with fits.open(spire_fn) as hdul:
	h_data = hdul[0].data
	h_header = hdul[0].header

p_band_stub = sys.argv[1]


p_map = hp.read_map(planck_dir+planck_fns["F"+p_band_stub], nest=True)
p_img_real = rhp.project_healpix_to_fits(p_map, h_data, h_header, nest=True)
p_img_real *= planck_unit_conversions["F"+p_band_stub]

convolve = False

if convolve:
	p_real_beam = pm.bandpass_beam_sizes["F"+p_band_stub]
	p_fabr_beam = 5
	conv_beam = pm.prepare_convolution(WCS(h_header), p_real_beam, p_fabr_beam, p_img_real.shape)
	p_img_real = pm.convolve_properly(p_img_real, conv_beam)


sky = pobs.DustModel(field_stub, "F"+p_band_stub, h_data, h_header, dust_directory)
p_img_fabr = sky.observe_planck("F"+p_band_stub, islarge=False, row_step=128)

if p_band_stub in ['545', '857']:
	vmin, vmax = 0.5, 3.5
else:
	vmin, vmax = -1, 2

plt.figure(figsize=(18, 9))
plt.subplot(131)
plt.imshow(np.log10(p_img_real), origin='lower', vmin=vmin, vmax=vmax)
if convolve:
	not_convolved = ""
else:
	not_convolved = "not "
plt.title("Planck %sGHz band data (%sconvolved up to 5')" % (p_band_stub, not_convolved))

plt.subplot(132)
p_ratio = p_img_real / p_img_fabr
plt.imshow(p_ratio, origin='lower', vmin=0.95, vmax=1.25)
plt.colorbar()
plt.title("Observed / Model Emission")

plt.subplot(133)
plt.imshow(np.log10(p_img_fabr), origin='lower', vmin=vmin, vmax=vmax)
plt.title("Planck %sGHz band model emission" % p_band_stub)

if convolve:
	un = ""
else:
	un = "un"

print("save_this() to save, save_ratio() to save as fits")
plt.tight_layout()
plt.show()

def save_this():
	plt.figure(figsize=(18, 9))
	plt.subplot(131)
	plt.imshow(np.log10(p_img_real), origin='lower', vmin=vmin, vmax=vmax)
	plt.title("Planck %sGHz band data (%sconvolved up to 5')" % (p_band_stub, not_convolved))

	plt.subplot(132)
	p_ratio = p_img_real / p_img_fabr
	plt.imshow(p_ratio, origin='lower', vmin=0.95, vmax=1.25)
	plt.colorbar()
	plt.title("Observed / Model Emission")

	plt.subplot(133)
	plt.imshow(np.log10(p_img_fabr), origin='lower', vmin=vmin, vmax=vmax)
	plt.title("Planck %sGHz band model emission" % p_band_stub)

	plt.tight_layout()
	plt.savefig("Figure_%scomp_%sconvolved.png" % (p_band_stub, un))

def save_ratio():
	save_fn = "Planck_ObsToModel_ratio_%s_%sconvolved.fits" % (p_band_stub, un)
	h_header['COMMENT'] = "THIS IS NOT HERSCHEL"
	h_header['COMMENT'] = "this is the ratio between Planck Obs / dust model at %s GHz" % p_band_stub
	fits.writeto(save_fn, p_ratio, h_header)


save_this()
save_ratio()

