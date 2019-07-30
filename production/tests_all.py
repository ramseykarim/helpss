import numpy as np
from astropy.io import fits
import utils_regrid as rgu
import path_config as cfg
import calc_pacs_offset as cpo
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
    hp_fn = cfg.PlanckConfig.light_map_filename('F353')
    fits_fn = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
    fits_fn += "PACS160um-image-remapped.fits"
    # open the healpix map with the convenience function
    m = rgu.open_healpix(hp_fn, nest=False)
    # get fits
    data, head = fits.getdata(fits_fn, header=True)
    result = rgu.healpix2fits(m, data, head, method='scipy')
    result = cfg.PlanckConfig.unit_conversion('F353', result)
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


def test_predict():
    pacs_fn = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
    pacs_fn += "PACS160um-image-remapped.fits"
    # pdata, phead = fits.getdata(pacs_fn, header=True)
    model = cpo.GNILCModel(pacs_fn)
    model.accumulate_masks()
    model.difference_to_target()
    return model

if __name__ == "__main__":
    model = test_predict()
