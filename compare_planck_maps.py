#!/jupiter/rkarim/anaconda3/envs/py36/bin/python
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import healpy as hp
import sys

planck_dir = "/n/sgraraid/filaments/data/"


filename_formulas = {
	"raw": ("HFI_SkyMap_", "-field-Int_2048_R3.00_full.fits"),
	"CMB-subtracted": ("HFI_CompMap_Foregrounds-nilc-", "_R3.00.fits"),
	"dust-only": ("COM_CompMap_Dust-GNILC-F", "_2048_R2.00.fits"),
}
def gen_planck_fn(frequency, map_type):
	first, second = filename_formulas[map_type]
	return f"{planck_dir}{first}{frequency}{second}"


"""
for s in filename_formulas:
	# confirmed; all are in MJy/sr
	head = fits.getheader(gen_planck_fn(857, s), 1)
	unit = head['TUNIT1']
	print(f"{s} data comes in {unit}")
sys.exit()
"""

#raw8 = hp.read_map(gen_planck_fn(857, "raw"), nest=False)
#cmbsub8 = hp.read_map(gen_planck_fn(857, "CMB-subtracted"), nest=False)
# dust8 = hp.read_map(gen_planck_fn(857, "dust-only"), nest=False)

raw5 = hp.read_map(gen_planck_fn(545, "raw"), nest=False)
#cmbsub5 = hp.read_map(gen_planck_fn(545, "CMB-subtracted"), nest=False)
# dust5 = hp.read_map(gen_planck_fn(545, "dust-only"), nest=False)

print(fits.getheader(gen_planck_fn(545, "raw"), 1))
hp.mollview(raw5, xsize=800) #, norm='log')#, min=-1, max=1)
plt.show()
sys.exit()

lon, lat = (156.5, 160.5), (-23., -19.)
projection = hp.projector.CartesianProj(xsize=800, ysize=800,
                                        lonra=np.array(lon),
                                        latra=np.array(lat))
projection.set_flip('astro')
print("center", projection.get_center(lonlat=1))
print("extent", projection.get_extent())
n_side = hp.get_nside(cmbsub8)
vec2pix_func = lambda x, y, z: hp.pixelfunc.vec2pix(n_side, x, y, z, nest=0)

img8 = projection.projmap((raw8 - cmbsub8)/raw8, vec2pix_func)
img5 = projection.projmap((raw5 - cmbsub5)/raw5, vec2pix_func)

"""
plt.figure(figsize=(18, 12))
plt.subplot(121)
# plt.imshow(np.log10(np.abs(img8)), origin='lower', vmin=-7, vmax=-4)
# plt.title("log10 of absolute value [(CMB-subtracted $-$ dust model)/CMB-subtracted] fraction (857 GHz)")
plt.imshow(img8, origin='lower')
plt.title("raw sky minus CMB-subtracted, MJy/sr (857 GHz)")
plt.colorbar()
plt.subplot(122)
# plt.imshow(np.log10(np.abs(img5)), origin='lower', vmin=-7, vmax=-4)
# plt.title("log10 of absolute value [(CMB-subtracted $-$ dust model)/CMB-subtracted] fraction (545 GHz)")
plt.imshow(img5, origin='lower')
plt.title("raw sky minus CMB-subtracted, MJy/sr (545 GHz)")
plt.colorbar()
plt.show()
"""
