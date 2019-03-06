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
import pickle
from scipy.interpolate import griddata
import matplotlib.patches as mpatches


"""
HEADS UP
This file was written when we were still using a txt file database for the Planck values
Now we use healpy, which does not run on the department machines (yet..)
This file has some useful functions but should not be used as a dependency
"""


herschel_band_stubs = ["PACS160um", "SPIRE250um", "SPIRE350um", "SPIRE500um"]
herschel_band_stub = herschel_band_stubs[2]
home_stub = "/n/sgraraid/filaments/data/"

# Per 1 specific filenames and stuff
dir_stub = "%sTEST4/Per/testregion1342190326JS/" % home_stub
field_stub = "Per1J"
c1_coord = SkyCoord("3:29:07.738 +31:20:28.05", frame=FK5, unit=(u.hourangle, u.deg))
c2_coord = SkyCoord("3:36:23.014 +31:13:04.54", frame=FK5, unit=(u.hourangle, u.deg))
c_coords = {0: None}

# Generator expressions for FITS image names
fn_2paramfit = "%sT4-absdiff-%s-4bandLErr.fits" % (dir_stub, field_stub)
fn_Herschel = lambda band: "%s%s-image.fits" % (dir_stub, band)
fn_cropped_Herschel = lambda band, n: "%s%s-image-cropped%d.fits" % (dir_stub, band, n) if n > 0 else fn_Herschel(band)
# Note we are using the FULL image, not regridded/convolved/anything
# We will find the corrective factor using this unperturbed image

# =====================================================
# Crop the original Herschel images
# =====================================================


def crop_original(herschel_band_index):
	"""
	Crop the images!
	Saves the fits file
	"""
	cropped_width = 3
	xbox, ybox = cropped_width*u.deg, cropped_width*u.deg

	print("load up", fn_Herschel(herschel_band_stubs[herschel_band_index]))

	with fits.open(fn_Herschel(herschel_band_stubs[herschel_band_index])) as hdul:
		d = hdul[0].data
		h = hdul[0].header
		w = WCS(h)
	w.cunit1, w.cunit2 = "Degrees", "Degrees"

	for j in c_coords:
		print("wrote to", fn_cropped_Herschel(herschel_band_stubs[herschel_band_index], j))
		d_c = Cutout2D(d, c_coords[j], (ybox, xbox), wcs=w)
		h.update(d_c.wcs.to_header())
		fits.writeto(fn_cropped_Herschel(herschel_band_stubs[herschel_band_index], j),
			d_c.data, h, overwrite=True)

	print("")
	# we have run the above function only on Per1 SPIRE 350 (Dec 12, 2018 2:49 pm)

# =====================================================
# Now regrid Planck to those images
# =====================================================

planck_band_stubs = ["857GHz-350um"]
planck_band_stub = planck_band_stubs[0]

# pickle file of Planck HFI 857 GHz data
fn_HFI_pkl = "%sPlanck-%s-FLUX.pkl" % (home_stub, planck_band_stub)
# Destinations for interpolated values
fn_cropped_HFI = lambda Pband, Hband, n: "%sHFI_SkyMap_857_2048_R2_%sres_crop%d.fits" % (home_stub, Hband, n)


def get_HFI():
	"""
	Get the HFI data; shape is (n, 3) where n is like 50 million
	3 columns are: RA, DEC, FLUX (ra/dec in degrees) (flux presumed to be MJy/sr)
	"""
	with open(fn_HFI_pkl, 'rb') as handle:
		HFI_data = pickle.load(handle)
	return HFI_data


def regrid_Planck(herschel_band_index, planck_band_index):
	"""
	Regrid Planck, of any band, to a Herschel field of any band
	Does both cutouts
	Saves the FITS file
	"""
	data = get_HFI()
	print("Loaded HFI map")
	for i in c_coords:
		with fits.open(fn_cropped_Herschel(herschel_band_stubs[herschel_band_index], i)) as hdul:
			d = hdul[0].data
			h = hdul[0].header
			w = WCS(h)
		print("Loaded %s" % fn_cropped_Herschel(herschel_band_stubs[herschel_band_index], i))
		# Create pixel coordinate list by looping through all the indices
		x_pix_array, y_pix_array = np.arange(d.shape[0]), np.arange(d.shape[1])
		destination_pixel_list = deque()
		for xp in x_pix_array:
			for yp in y_pix_array:
				destination_pixel_list.append(np.array([xp, yp]))
		destination_pixel_list = np.array(destination_pixel_list)
		# NOTE: numpy 0 indexes the pixels, hence the 0 argument
		destination_coord_list = w.wcs_pix2world(destination_pixel_list, 0)
		extra_width = 0.25
		# Figure out sky coordinate bounds of the desired region
		minmax_list = []
		for j in range(2):
			# Loop through 2 coords: RA, DEC
			# This yields an RA, DEC min max pair: [RA[min, max], DEC[min, max]]
			coords = destination_coord_list[:, j]
			coord_min, coord_max = np.min(coords), np.max(coords)
			diff = coord_max - coord_min
			minmax_list.append([coord_min - diff*extra_width, coord_max + diff*extra_width])
		interp_box = np.where(
			(data[:, 0] > min(minmax_list[0])) & (data[:, 0] < max(minmax_list[0])) &
			(data[:, 1] > min(minmax_list[1])) & (data[:, 1] < max(minmax_list[1])) )
		data_subset = data[interp_box]
		# All 3 available interpolation methods
		methods = ['nearest', 'linear', 'cubic']
		# NEAREST seems most honest
		method = methods[0]
		print("Set up interpolation, interpolating...", end="")
		interpolated_values = griddata(data_subset[:, :2],
			data_subset[:, 2], destination_coord_list, method=method)
		interpolated_data = np.empty(d.shape)
		print("done. Assigning values.")
		for j in range(destination_pixel_list.shape[0]):
			xy_pix = destination_pixel_list[j, :]
			x_pix, y_pix = xy_pix[0], xy_pix[1]
			interpolated_data[x_pix, y_pix] = interpolated_values[j]
		h['HISTORY'] = "THIS IS NOT HERSCHEL --- THIS IS PLANCK 857 GHz"
		h['COMMENT'] = "Ramsey edited this file"
		# Write the nearest neighbor version out to file
		fits.writeto(fn_cropped_HFI(planck_band_stubs[planck_band_index],
			herschel_band_stubs[herschel_band_index], i),
			interpolated_data, h, overwrite=True)
		print("Saved interpolated map as %s" % fn_cropped_HFI(planck_band_stubs[planck_band_index],
			herschel_band_stubs[herschel_band_index], i))

if __name__ == '__main__':
	# put indices from 0-3 in the to-do-list
	to_do_list = [2]
	for i in to_do_list:
		# crop_original(i)
		regrid_Planck(i, 0)
