import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

wavelengths = [160, 250, 350, 500]

band_stubs = {
	160: "PACS160um",
	250: "SPIRE250um",
	350: "SPIRE350um",
	500: "SPIRE500um",
}

err_mod = 5.

per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
img = "-image"
err = "-error"
suffix = "-remapped-conv.fits"
def get_img_stub(wavelength):
	band_stub = band_stubs[wavelength]
	if wavelength == 160:
		band_stub += "-plus045"
	return f"{per1_dir}{band_stub}{img}{suffix}"

def get_new_err_stub(wavelength):
	band_stub = band_stubs[wavelength] + "-plus{:04.1f}pct".format(err_mod)
	return f"{per1_dir}{band_stub}{err}{suffix}"

def get_err_stub(wavelength):
	band_stub = band_stubs[wavelength]
	return f"{per1_dir}{band_stub}{err}{suffix}"

for wl in wavelengths:
	i_data = fits.getdata(get_img_stub(wl))
	e_data, hdr = fits.getdata(get_err_stub(wl), header=True)
	hdr['comment'] = "Added {:04.1f} percent of flux to this map (RLK)".format(err_mod)
	e_data = np.sqrt(e_data**2 + (i_data*(err_mod/100.))**2)
	fits.writeto(get_new_err_stub(wl), e_data, hdr)
	print("Wrote ", get_new_err_stub(wl))
