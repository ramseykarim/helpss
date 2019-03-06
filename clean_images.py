import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

rcc = "-remapped-conv.fits"
rccl = "-remapped-conv-clean.fits"
bands =[
    "PACS160um-",
    "SPIRE250um-",
    "SPIRE350um-",
    "SPIRE500um-",
]
img, err = "image", "error"
ie = [img, err]

def constr_fn(b, m_type):
    return b + m_type + rcc

def constr_sn(b, m_type):
    return b + m_type + rccl

def load_map(fn):
    data, header = fits.getdata(fn, header=True)
    return data, header

images = []
headers = []
masks = []
masks2 = []
for b in bands:
    ds, hs = [], []
    for m_type in ie:
        d, h = load_map(constr_fn(b, m_type))
        ds.append(d), hs.append(h)
    images.append(ds), headers.append(hs)
    mask = np.isnan(ds[0]) | np.isnan(ds[1]) | (ds[1] == 0)
    masks.append(mask)
    mask2 = 1*np.isnan(ds[0]) + 1*np.isnan(ds[1]) + 1*(ds[0]==0) + 1*(ds[1]==0)
    masks2.append(mask2)

masks = np.array(masks)
masks2 = np.array(masks2)
mask = np.any(masks, axis=0)
mask2 = np.sum(masks2, axis=0)

for i, b in enumerate(bands):
    ds, hs = images[i], headers[i]
    for j, m_type in enumerate(ie):
        d, h = ds[j], hs[j]
        d[mask] = np.nan
        sn = constr_sn(b, m_type)
        print("Saving %s" % sn)
        fits.writeto(sn, d, h, overwrite=True)


