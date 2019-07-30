import numpy as np
from astropy.io import fits
import utils_regrid as rgu
import path_config as cfg
import matplotlib.pyplot as plt


"""
These tests are currently designed for my laptop but will
be re-written for the department computers. (July 29, 2019)
"""


def test_healpy_regrid():
    """
    Test HEALPix regrid via healpy/scipy
    TEST PASSED by visual inspection
    """
    # set some paths (laptop local)
    directory = "/home/ramsey/Documents/Research/Filaments/"
    hp_fn = directory + "HFI_SkyMap_857-field-Int_2048_R3.00_full.fits"
    fits_fn = directory + "HFI_SkyMap_857_resSPIRE350.fits"
    # open the healpix map with the convenience function
    m = rgu.open_healpix(hp_fn, nest=False)
    # get fits
    data, head = fits.getdata(fits_fn, header=True)
    result = rgu.healpix2fits(m, data, head, method='scipy')
    plt.imshow(result, origin='lower')
    plt.show()


def test_bandpass_config():
    # get center of PACS 160 micron band
    f0 = cfg.HerschelConfig.bandpass_center('PACS160um')
    # check that the result is within 10 microns of 160 (should be)
    assert (160 - (1e6 * 3e8/f0)) < 10
    try:
        cfg.HerschelConfig.bandpass_center('F857')
    except RuntimeError:
        # Should raise this error
        pass
    else:
        assert False
    finally:
        # now make sure it gets the Planck bandpass center alright
        f0 = cfg.PlanckConfig.bandpass_center('F857')
        assert ((f0 * 1e-9) - 857) < 1
    frq, wgt = cfg.PlanckConfig.bandpass_profile('F545')
    plt.plot(frq, wgt)
    plt.xscale('log')
    plt.show()


if __name__ == "__main__":
    test_bandpass_config()
