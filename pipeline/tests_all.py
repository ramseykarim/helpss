import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import utils_regrid as rgu
import path_config as cfg
import calc_offset as cpo
import matplotlib.pyplot as plt
from modify_fits import add_offset, add_systematic_error


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
    pdata, phead = fits.getdata(pacs_fn, header=True)
    # set pixel_scale_arcsec to 300 so this runs quickly for testing purposes
    # this is just the pixel scale for an intermediate step; the final grid
    #  is the same as the input PACS grid
    model = cpo.GNILCModel(pdata, phead, pixel_scale_arcsec=300)
    model.get_offset(full_diagnostic=True)
    return model


def test_add_flux_mod_error():
    pacs_path = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
    pacs_fn = pacs_path+"PACS160um-image-remapped-conv.fits"
    pacs_err = pacs_path+"PACS160um-error-remapped-conv.fits"
    add_offset(45, pacs_fn)
    add_systematic_error(0.05, pacs_err, pacs_fn)
    spire_fn = lambda x: pacs_path+"SPIRE{:d}um-image-remapped-conv.fits".format(x)
    spire_err = lambda x: pacs_path+"SPIRE{:d}um-error-remapped-conv.fits".format(x)
    for b in [250, 350, 500]:
        add_systematic_error(0.015, spire_err(b), spire_fn(b))

if __name__ == "__main__":
    test_add_flux_mod_error()
    # test_predict()
