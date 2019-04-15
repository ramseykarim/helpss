import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import sys

planck_857 = "/n/sgraraid/filaments/data/"
planck_857 += "HFI_SkyMap_857-field-Int_2048_R3.00_full.fits"

planck_857 = hp.read_map(planck_857, nest=False)
#hp.mollview(planck_857)
#plt.show()

lon = float(sys.argv[1]), float(sys.argv[2])
lat = float(sys.argv[3]), float(sys.argv[4])


projection = hp.projector.CartesianProj(xsize=400, ysize=400,
                                        lonra=np.array(lon),
                                        latra=np.array(lat))
projection.set_flip('astro')
print("center", projection.get_center(lonlat=1))
print("extent", projection.get_extent())
n_side = hp.get_nside(planck_857)
vec2pix_func = lambda x, y, z: hp.pixelfunc.vec2pix(n_side, x, y, z, nest=0)
img = projection.projmap(planck_857, vec2pix_func)
plt.imshow(np.log10(img), origin='lower')
plt.show()
