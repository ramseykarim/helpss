"""
Investigate the Ser-Main temperature gradient by checking the Planck GNILC model T maps

Created: August 24, 2022
"""
__author__ = "Ramsey Karim"

import os.path

import numpy as np

from astropy.io import fits
from astropy.wcs import WCS

from reproject import reproject_from_healpix

import matplotlib.pyplot as plt


gnilc_path = "/n/sgraraid/filaments/data/filterInfo_PlanckHerschel"

gnilc_t_filename = "COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.00.fits"
gnilc_beta_filename = "COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.00.fits"
gnilc_tau_filename = "COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.00.fits"
component_filenames = [gnilc_t_filename, gnilc_beta_filename, gnilc_tau_filename]
components = ['t', 'beta', 'tau']
full_component_name = ['temperature', 'spectral index', 'opacity']
units = ['T (K)', 'spectral index', 'opacity (353um)']


choose_component = 'tau'
gnilc_component_filename = component_filenames[components.index(choose_component)]
gnilc_component_full_filename = os.path.join(gnilc_path, gnilc_component_filename)


ser_aqu_path = "/n/sgraraid/filaments/Serpens-Aquila/Herschel/processed"
ser_main_obsID = "1342206676"
aqu_east_obsID = "1342206694"

# For footprint / coords
pacs_image_filename = "PACS160um-image-remapped-conv.fits"
pacs_image_full_filename = os.path.join(ser_aqu_path, ser_main_obsID, pacs_image_filename)
# Get WCS
pacs_data, pacs_hdr = fits.getdata(pacs_image_full_filename, header=True)
pacs_wcs = WCS(pacs_hdr)

ser_main_component_map, _ = reproject_from_healpix(gnilc_component_full_filename, pacs_wcs, shape_out=pacs_data.shape)

if True: # already did this on 2022-08-24
    new_hdr = pacs_wcs.to_header()
    new_hdr['BUNIT'] = 'K'
    new_hdr['COMMENT'] = f'Planck GNILC dust model {full_component_name[components.index(choose_component)]}'
    new_hdr['COMMENT'] = f'regridded to SPIRE500um obsID {ser_main_obsID}'
    new_hdr['HISTORY'] = 'Created by Ramsey Karim, 2022-08-24'

    savename = os.path.join(ser_aqu_path, ser_main_obsID, f"{full_component_name[components.index(choose_component)].replace(' ', '-')}_map_planck_gnilc.fits")
    hdu = fits.PrimaryHDU(data=ser_main_component_map, header=new_hdr)
    hdu.writeto(savename)


pacs_footprint = np.isfinite(pacs_data).astype(float)

fig = plt.figure(figsize=(15, 5))
ax1 = plt.subplot(121, projection=pacs_wcs)
ax2 = plt.subplot(122, projection=pacs_wcs)

im1 = ax1.imshow(pacs_data, origin='lower', vmin=-100, vmax=120)
fig.colorbar(im1, ax=ax1, label='(MJy/sr)')
ax1.set_title("PACS 160um map")
im2 = ax2.imshow(ser_main_component_map, origin='lower')
fig.colorbar(im2, ax=ax2, label=units[components.index(choose_component)])
ax2.contour(pacs_footprint, levels=[0.5])
ax2.set_title(f"Planck GNILC dust {choose_component.upper()}")

plt.savefig(f"/home/rkarim/Pictures/2022-08-24/ser_main_{choose_component.upper()}_map.png")
