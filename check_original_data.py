import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys


psrcPSW = "/n/sgraraid/filaments/data/Per/1342190326/level2_5/psrcPSW/"
psrcPSW += "hspirepsw271_25pmp_0328_p3051_1342190326_1342190327_1462439413361.fits.gz"

extdPSW = "/n/sgraraid/filaments/data/Per/1342190326/level2_5/extdPSW/"
extdPSW += "hspirepsw271_25pxmp_0328_p3051_1342190326_1342190327_1462439341039.fits.gz"

extdPMW = "/n/sgraraid/filaments/data/Per/1342190326/level2_5/extdPMW/"
extdPMW += "hspirepmw271_25pxmp_0328_p3051_1342190326_1342190327_1462439337908.fits.gz"

extdPLW = "/n/sgraraid/filaments/data/Per/1342190326/level2_5/extdPLW/"
extdPLW += "hspireplw271_25pxmp_0328_p3051_1342190326_1342190327_1462439335742.fits.gz"

PER_250 = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
PER_250 += "SPIRE250um-image.fits"
PER_350 = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
PER_350 += "SPIRE350um-image.fits"
PER_500 = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
PER_500 += "SPIRE500um-image.fits"

GB_250 = "/home/rkarim/Research/Filaments/Data/perseus04-250.fits"

JSR = "/n/sgraraid/filaments/data/Per/1342190326/level2_5/HPPJSMAPR/"
JSR += "hpacs_25HPPJSMAPR_0330_p3056_00_v1.0_1470647587537.fits.gz"

PER_160 = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
PER_160 += "PACS160um-image.fits"

def sqArcsecToSterad():
    degToRad = np.pi/180.
    sqDegToSterad = degToRad*degToRad
    sqArcsec_per_sqDeg = 3600.*3600.
    return sqDegToSterad / sqArcsec_per_sqDeg

def beamToSterad(beam_arcsec):
    return beam_arcsec * sqArcsecToSterad()

SPIRE_BEAM_SIZES_ARCSEC = {
    250: 469.35,
    350: 831.27,
    500: 1804.31,
}

def jyPerBeam_to_MJyPerSr(wavelength, jy_per_beam):
    sterad_per_beam = beamToSterad(SPIRE_BEAM_SIZES_ARCSEC[wavelength])
    return jy_per_beam / (sterad_per_beam * 1e6)

def jyPerPix_to_MJyPerSr(delta_pix_deg, jy_per_pix):
    pixel_unit_conv = (delta_pix_deg**2) * ((np.pi/180.)**2.) * 1e6
    return jy_per_pix / pixel_unit_conv


#pdata = fits.getdata(psrcPSW)
#pdata = jyPerBeam_to_MJyPerSr(250, pdata)
edata = fits.getdata(extdPLW)
sdata = fits.getdata(PER_500)
#gdata = fits.getdata(GB_250)
"""
jdata, jhead = fits.getdata(JSR, header=True)
pix_width = (abs(jhead['CDELT1']) + abs(jhead['CDELT2']))/2.
jdata = jyPerPix_to_MJyPerSr(pix_width, jdata)
pdata = fits.getdata(PER_160)
diff_map = pdata - jdata
print("MAX {:.3E}, MIN {:.3E}".format(
        np.nanmin(diff_map), np.nanmax(diff_map)))
"""
plt.imshow(edata - sdata, origin='lower') #, vmin=33, vmax=36)
plt.show()
