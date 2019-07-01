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
# these two are older; not sure what parts of the sky they are
c1_coord = SkyCoord("3:29:07.738 +31:20:28.05", frame=FK5, unit=(u.hourangle, u.deg))
c2_coord = SkyCoord("3:36:23.014 +31:13:04.54", frame=FK5, unit=(u.hourangle, u.deg))
# this is newer (April 5, 2019) and is most of the Per1 region
c3_coord = SkyCoord("3:30:18.406 +30:53:09.30", frame=FK5, unit=(u.hourangle, u.deg))
# this is also newer, smaller box to try to run the CMB manticore
c4_coord = SkyCoord("3:30:26.584 +29:58:56.26", frame=FK5, unit=(u.hourangle, u.deg)) # 0
c5_coord = SkyCoord("3:33:06.408 +31:08:03.71", frame=FK5, unit=(u.hourangle, u.deg)) # 1
# even newer (July 1, 2019), these two highlight some low-Tc/hi-Nc noise holdouts
# this one is just above B1
c6_coord = SkyCoord("3:33:53.142 +31:19:30.83", frame=FK5, unit=(u.hourangle, u.deg)) # 6
# this is right at NGC 1333
c7_coord = SkyCoord("3:29:09.435 +31:13:13.51", frame=FK5, unit=(u.hourangle, u.deg)) # 7

c_coords = {6: (c6_coord, 0.14, 0.16), 7: (c7_coord, 0.12, 0.17)}

# Generator expressions for FITS image names
fn_2paramfit = "%sT4-absdiff-%s-4bandLErr.fits" % (dir_stub, field_stub)
fn_Herschel = lambda band: "%s%s-image.fits" % (dir_stub, band)
fn_cropped_Herschel = lambda band, n: "%s%s-image-cropped%d.fits" % (dir_stub, band, n) if n > 0 else fn_Herschel(band)
# Note we are using the FULL image, not regridded/convolved/anything
# We will find the corrective factor using this unperturbed image

# =====================================================
# Crop the original Herschel images
# =====================================================


def crop_original(coords, source_filename, target_filename, cropped_width=3):
	"""
	Crop the images!
	Saves the fits file
	"""
	try:
		# can use tuple (dx, dy) for width
		# dy should be d_dec, dx should be d_ra
		xbox, ybox = cropped_width[0]*u.deg, cropped_width[1]*u.deg
	except TypeError:
		xbox, ybox = cropped_width*u.deg, cropped_width*u.deg


	print("load up", source_filename)

	with fits.open(source_filename) as hdul:
		d = hdul[0].data
		h = hdul[0].header
		w = WCS(h)
	w.cunit1, w.cunit2 = "Degrees", "Degrees"

	print("wrote to", target_filename)
	d_c = Cutout2D(d, coords, (ybox, xbox), wcs=w)
	h.update(d_c.wcs.to_header())
	fits.writeto(target_filename,
		     d_c.data, h, overwrite=True)
	print("")
	# we have run the above function only on Per1 SPIRE 350 (Dec 12, 2018 2:49 pm)
	# this is definitely no longer true and hasn't been for a while (Jul 1, 2019 2:16 pm)

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
	img1 = "/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/PACS160um-plus045-image-remapped-conv"
	img2 = "/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/SPIRE250um-image-remapped-conv"
	img3 = "/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/SPIRE350um-image-remapped-conv"
	img4 = "/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/SPIRE500um-image-remapped-conv"
	err1 = "/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/PACS160um-plus05.0pct-error-remapped-conv"
	err2 = "/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/SPIRE250um-plus05.0pct-error-remapped-conv"
	err3 = "/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/SPIRE350um-plus05.0pct-error-remapped-conv"
	err4 = "/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/SPIRE500um-plus05.0pct-error-remapped-conv"
	all_files = [img1, img2, img3, img4, err1, err2, err3, err4]
	for f in all_files:
		for i in c_coords:
			coord, dx, dy = c_coords[i]
			source_filename = "/n" + f + ".fits"
			target_filename = "/n" + f + "-crop{}.fits".format(str(int(i)))
			crop_original(coord, source_filename, target_filename, cropped_width=(dx, dy))


