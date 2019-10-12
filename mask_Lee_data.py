import numpy as np
import matplotlib.pyplot as plt
import glob, sys
from astropy.io import fits
from pipeline.calc_offset import flquantiles

"""
Testing masks on Lee's data (small frames)
Creation date: Oct 12, 2019
"""
__author__ = "Ramsey Karim"

lee_dir = "/n/sgraraid/filaments/lgm/manticore-test/L723/"
comp_stub = "TEST1-"
comp1_dir = comp_stub + "1comp/"
comp2_dir = comp_stub + "2comp*/"
fits_stub = "*.fits"
comp1_fn = glob.glob(lee_dir+comp1_dir+fits_stub)[0]
comp2_fns = glob.glob(lee_dir+comp2_dir+fits_stub)

nominal_2p_fn = comp1_fn
nominal_3p_fn = comp2_fns[0]


with fits.open(nominal_2p_fn) as hdul:
    T2 = hdul[1].data
    N2 = hdul[3].data
    # pacs_mask True where VALID
    pacs_mask = ~np.isnan(hdul[12].data)

with fits.open(nominal_3p_fn) as hdul:
    Nh = hdul[3].data
    Tc = hdul[5].data
    Nc = hdul[7].data

def imglim(img):
    lims = flquantiles(img[~np.isnan(img)].ravel(), 10)
    return {v: l for v, l in zip(('vmin', 'vmax'), lims)}

pltimg = Tc
maskN = (N2 > 2.5e21) & (Nc > 0)
lims = imglim(pltimg[pacs_mask&maskN])
pltimg_c = pltimg.copy()
pltimg_c[~(maskN&pacs_mask)] = np.nan
"""
It seems like in single_temperature_varying.py,
I used a cutoff of N2 = 1.5e21.
This does not work here, since most of the frame (in L723)
is above 1.5e21 to begin with.

How about N2 > 2e21? This sort of works.
It still lets in a lot of the stripey noise.

I found that Nc == 0 for a lot of the stripey noise, so & Nc>0 helps.
This still lets in some Tc~5K.

How about examining the N2 map under these masks? Reveals N2>2.5e21
could help.
However, there is STILL Tc~5K. Maybe we need the Nh smoothing iteration.
"""

plt.figure(figsize=(16, 9))
plt.subplot(121)
plt.imshow(pltimg, origin='lower', **lims)
plt.colorbar()
plt.subplot(122)
plt.imshow(pltimg_c, origin='lower', **lims)
plt.colorbar()
plt.show()
