"""
Created: August 3, 2021
Checking the PACS 160 filter curves.

The ones I've been using were copied and pasted out of the .cc files in Kevin's
Manticore code.
There are some on this website: http://archives.esac.esa.int/hsa/legacy/ADP/PACS/PACS-P_filter_curves/
And there are some on the SVO filter service: http://svo2.cab.inta-csic.es/svo/theory/fps/index.php?id=Herschel/Pacs.red&&mode=browse&gname=Herschel&gname2=Pacs#filter

I will be comparing all of them right now.
"""
__author__ = "Ramsey Karim"

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.io import votable
from astropy import units as u


def load_local():
    fn_local = "/home/ramsey/Documents/Research/Filaments/filterInfo_PlanckHerschel/PACS160um_fromManticore.dat"
    data_local = np.genfromtxt(fn_local)
    freq_local, R_local = data_local[:, 0], data_local[:, 1]
    freq_local = freq_local * u.Hz
    wl_local = freq_local.to(u.micron, equivalencies=u.spectral()).to_value()
    return wl_local, R_local, "Manticore"


def load_archive():
    fn_archive = "/home/ramsey/Downloads/PACS_phot_filterTrans_red.fits.gz"
    with fits.open(fn_archive) as hdul:
        wl_archive = hdul[1].data['wavelength']
        R_archive = hdul[1].data['transmission']
    return wl_archive, R_archive, "Archive"


def load_SVO():
    fn_SVO = "/home/ramsey/Downloads/Herschel.Pacs.red.xml"
    table = votable.parse_single_table(fn_SVO).array
    wl_SVO = table['Wavelength'].data # Angstroms
    R_SVO = table['Transmission'].data
    wl_SVO = wl_SVO * u.Angstrom.to(u.micron)
    return wl_SVO, R_SVO, "SVO"


curves = [load_local, load_archive, load_SVO]
ls = ['-.', '-', '--']
ax = plt.subplot(111)
for i, load_func in enumerate(curves):
    wl, R, label = load_func()
    ax.plot(wl, R, linestyle=ls[i], label=label)
ax.set_xlabel("wavelength (micron)")
ax.set_ylabel("PACS 160 transmission")
ax.legend()
ax.set_xlim([100, 275])
plt.show()
# plt.savefig("/home/ramsey/Pictures/2021-08-03-work/PACS160_transmission_SVO.png")
