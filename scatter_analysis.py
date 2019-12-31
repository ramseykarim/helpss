"""
Need to document all my dumb script files better!!!!
This file was created on December 30 2019 on Amtrak westbound Capitol Corridor train 543
Purpose: sequel to reinit_analysis, whose purpose was to make a lot of scatter plots
of pixels' solutions under various conditions.
Tracy (Dec 30 2019 emails) was looking at similar plots and it inspired me to
take a second look at these, with updated data and a fresh perspective.

Big reminder to just take freakin' notes in comments here so I don't run around
figuring out what I've already done 50x every time I open my laptop i stg
"""
__author__ = "Ramsey Karim"

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

singlecomp_fn = "./pipeline/full-1.5.1-Per1-pow-750-0.05625-2.10.fits"
twocomp_fn = "./pipeline/full-1.5.1-Per1-pow-750-0.05625-2.10,pow-750-0.05625-1.80-Th15.95.fits"

"""
Lee has been making "total N" maps: N_tot = N_cold + 2*N_hot
I'll try taking a look at one of these first.
Then I'll recreate Tracy's 1-comp N vs total N

Christ this looks like a jupyter notebook. still not giving into the peer pressure.
"""

with fits.open(singlecomp_fn) as hdul:
    N1 = hdul[3].data
    header = hdul[3].header
    w = WCS(header)
with fits.open(twocomp_fn) as hdul:
    Nh = hdul[3].data
    Nc = hdul[7].data

Ntot = Nc + 2*Nh

"""
I'll make one mask to rule them all
Uses some N-single cutoff in the 1-3e21 range. Lee uses 3e21, I use 1.5e21, no one
has any justification either way. Lee and I both use ours because
they look good.
"""
MASK = N1 > 3e21

def plot_total_N():
    plt.figure(figsize=(16, 9))
    kwargs = dict(origin='lower', vmin=1e21, vmax=10e21)
    plt.subplot(121)
    plt.imshow(N1, **kwargs)
    plt.title("Single component N")
    plt.subplot(122)
    plt.imshow(Ntot, **kwargs)
    plt.title("Total N")
    plt.show()

def plot_save_combined_N():
    Ncombined = np.full(Ntot.shape, np.nan)
    Ncombined[MASK] = Ntot[MASK]
    Ncombined[~MASK] = N1[~MASK]
    plt.figure(figsize=(16, 9))
    kwargs = dict(origin='lower', vmin=1e21, vmax=10e21)
    plt.subplot(111)
    plt.imshow(Ncombined, **kwargs)
    plt.title("Combined N")
    h = fits.Header()
    h.update(w.to_header())
    h['COMMENT'] = "Combined single/total N at cutoff singleN>3e21"
    fits.writeto("./pipeline/combinedN_ramsey.fits", Ncombined, header=h, overwrite=True)
    plt.show()

print('nothing doing')
