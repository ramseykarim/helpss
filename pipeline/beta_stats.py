"""
Quickly write out some stats given the Planck spectral index map saved in
a given directory.

This script is used similarly to calibrate.py: just use it in the directory
with the files planck-beta.fits and calib-mask.fits

Created: May 26, 2021
"""
__author__ = "Ramsey Karim"

import os
import numpy as np
from astropy.io import fits

beta = fits.getdata("planck-beta.fits")
mask = fits.getdata("calib-mask.fits")

mean_full = np.nanmean(beta)
median_full = np.nanmedian(beta)
std_full = np.nanstd(beta)

mean_mask = np.nanmean(beta[mask])
median_mask = np.nanmedian(beta[mask])
std_mask = np.nanstd(beta[mask])

with open('beta-stats.txt', 'w') as f:
    f.write(f"MEAN (FULL): {mean_full:.2f}")
    f.write(f"MEAN (MASK): {mean_mask:.2f}")

    f.write(f"MEDIAN (FULL): {median_full:.2f}")
    f.write(f"MEDIAN (MASK): {median_mask:.2f}")

    f.write(f"STDDEV (FULL): {std_full:.2f}")
    f.write(f"STDDEV (MASK): {std_mask:.2f}")
