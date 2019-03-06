#!/astromake/opt/python/anaconda3/bin/python
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord, FK5
from astropy import units as u
import matplotlib.pyplot as plt
import sys
from collections import deque

band_stubs = ["PACS160um", "SPIRE250um", "SPIRE350um", "SPIRE500um"]
band_stub = band_stubs[0]
fnPer1 = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"+band_stub+"-plus043-image-remapped-conv.fits"
fnPer1_2p = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/T4-absdiff-Per1J-4bandLErr.fits"
fnPer1_c1 = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"+band_stub+"-plus043-image-remapped-conv-cropped1.fits"
fnPer1_c2 = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"+band_stub+"-plus043-image-remapped-conv-cropped2.fits"
fnSer1 = "/n/sgraraid/filaments/data/TEST4/Ser/testregion1342186278JS/"+band_stub+"-image-remapped.fits"
fnHFI = "/n/sgraraid/filaments/data/HFI_SkyMap_857_2048_R2.02_full.fits"
fnC = "crop.fits"
fcCHFI = "cropHFI.fits"

width = 0.5
Per1_c1 = SkyCoord("3:29:07.738 +31:20:28.05", frame=FK5, unit=(u.hourangle, u.deg))
Per1_c2 = SkyCoord("3:36:23.014 +31:13:04.54", frame=FK5, unit=(u.hourangle, u.deg))
xbox, ybox = width*u.deg, width*u.deg
with fits.open(fnPer1) as hdul:
    d = hdul[0].data
    h = hdul[0].header
    w = WCS(h)
    w.cunit1, w.cunit2 = "Degrees", "Degrees"
    d1 = Cutout2D(d, Per1_c1, (ybox, xbox), wcs=w)
    h.update(d1.wcs.to_header())
    fits.writeto(fnPer1_c1, d1.data, h)
    d2 = Cutout2D(d, Per1_c2, (ybox, xbox), wcs=w)
    h.update(d2.wcs.to_header())
    fits.writeto(fnPer1_c2, d2.data, h)

"""
hdu_l = fits.HDUList()
count = 0
with fits.open(fnPer1_2p) as hdul:
    for hdu in hdul:
        if count == 0:
            h = fits.PrimaryHDU(header=hdu.header)
        else:
            w = WCS(hdu.header)
            w.cunit1, w.cunit2 = "Degrees", "Degrees"
            d1 = Cutout2D(hdu.data, Per1_c1, (ybox, xbox), wcs=w)
            hdu.header.update(d1.wcs.to_header())
            h = fits.ImageHDU(data=d1.data, header=hdu.header)
            hdu_l.append(h)
        count += 1
hdu_l.writeto(fnPer1_c1)
print("wrote c1")

hdu_l = fits.HDUList()
count = 0
with fits.open(fnPer1_2p) as hdul:
    for hdu in hdul:
        if count == 0:
            h = fits.PrimaryHDU(header=hdu.header)
        else:
            w = WCS(hdu.header)
            w.cunit1, w.cunit2 = "Degrees", "Degrees"
            d2 = Cutout2D(hdu.data, Per1_c2, (ybox, xbox), wcs=w)
            hdu.header.update(d2.wcs.to_header())
            h = fits.ImageHDU(data=d2.data, header=hdu.header)
            hdu_l.append(h)
        count += 1
hdu_l.writeto(fnPer1_c2)
print("wrote c1")
"""
"""
hdu, hdr = fits.getdata(fnHFI, header=True)
"""

#print((ybox, xbox))
#print(Per1_c1)
#print(Per1_c2)


#plt.imshow(hdu, origin='lower', vmin=0, vmax=100)
#plt.show()
"""
for d, h in zip(all_data, all_headers):
    w = WCS(h)
    wcss.append(w)
    w.cunit1, w.cunit2 = "Degrees", "Degrees"

    d1 = Cutout2D(d, Per1_c1, (ybox, xbox), wcs=w)
    h.update(d1.wcs.to_header())
    fits.writeto(fnPer1_c1, d1.data, h, overwrite=True)

    d2 = Cutout2D(d, Per1_c2, (ybox, xbox), wcs=w)
    h.update(d2.wcs.to_header())
    fits.writeto(fnPer1_c2, d2.data, h, overwrite=True)
"""
#hdu_crop1 = Cutout2D(hdu, Per1_c1, (ybox, xbox), wcs=WCS(hdr))
#hdu_crop2 = Cutout2D(hdu, Per1_c2, (ybox, xbox), wcs=WCS(hdr))
#wcs_cropped1 = hdu_crop1.wcs
#wcs_cropped1.cunit1, wcs_cropped1.cunit2 = "Degrees", "Degrees"
#wcs_cropped2 = hdu_crop2.wcs
#wcs_cropped2.cunit1, wcs_cropped2.cunit2 = "Degrees", "Degrees"

#hdr.update(wcs_cropped1.to_header())
#hdr['COMMENT'] = "I cropped this to what should be a (0.5 deg)^2 box"
#fits.writeto(fnPer1_c1, hdu_crop1.data, hdr, overwrite=True)
#hdr.update(wcs_cropped2.to_header())
#hdr['COMMENT'] = "I cropped this to what should be a (0.5 deg)^2 box"
#fits.writeto(fnPer1_c2, hdu_crop2.data, hdr, overwrite=True)

# remap:
# "/astromake/opt/wcstools/bin/"+remap
# spawn,dirREMAP+"remap -n -9999 -f "+reference+" -o "+fileRESAMP+" "+fileNORESAMP
# Tracy fixes header at this point
