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



band_stubs = ["PACS160um", "SPIRE250um", "SPIRE350um", "SPIRE500um"]
band_stub = band_stubs[0]

# # wcs sources (x/y positions on which to interpolate)
fnPer1 = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"+band_stub+"-plus043-image-remapped-conv.fits"
fn_Per1_c1 = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"+band_stub+"-image-remapped-conv-cropped1.fits"
fn_Per1_c2 = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"+band_stub+"-image-remapped-conv-cropped2.fits"

# # table of positions and values to interpolate
fn_HFI_tbl = '/n/sgraraid/filaments/data/Planck-857GHz-350um-FLUX.tbl'
# # same table but a pickle file (half the size)
fn_HFI_pkl = '/n/sgraraid/filaments/data/Planck-857GHz-350um-FLUX.pkl'

# # destination for interpolated values
fnC_HFI_c1 = "/n/sgraraid/filaments/data/HFI_SkyMap_857_2048_R2_crop1.fits"
fnC_HFI_c2 = "/n/sgraraid/filaments/data/HFI_SkyMap_857_2048_R2_crop2.fits"

# # get the HFI data
# # shape is (n, 3) where n is like 50 million
# # 3 columns are: RA, DEC, FLUX (ra/dec in degrees) (flux presumed to be MJy/sr)
def get_HFI():
	with open(fn_HFI_pkl, 'rb') as handle:
		HFI_data = pickle.load(handle)
	return HFI_data

# # convenient lists
source_filenames = [fn_Per1_c1, fn_Per1_c2]
destination_filenames = [fnC_HFI_c1, fnC_HFI_c2]



for i in range(2):
	# # open the file with the desired positions
	with fits.open(source_filenames[i]) as hdul:
		d = hdul[0].data
		h = hdul[0].header
		w = WCS(h)

	# # pixel coordinate list
	x_pix_array, y_pix_array = np.arange(d.shape[0]), np.arange(d.shape[1])
	destination_pixel_list = deque()
	for xp in x_pix_array:
		for yp in y_pix_array:
			destination_pixel_list.append(np.array([xp, yp]))
	destination_pixel_list = np.array(destination_pixel_list)
	# # NOTE: numpy 0 indexes the pixels, hence the 0 argument
	destination_coordinate_list = w.wcs_pix2world(destination_pixel_list, 0)

	extra_width = 0.25

	minmax_list = []
	for j in range(2):
		coords = destination_coordinate_list[:, j]
		coord_min, coord_max = np.min(coords), np.max(coords)
		diff = coord_max - coord_min
		minmax_list.append([coord_min - diff*extra_width, coord_max + diff*extra_width])
	 # # this yields an RA, DEC min max pair: [RA[min, max], DEC[min, max]]

	data = get_HFI()
	interp_box = np.where(
		(data[:, 0] > min(minmax_list[0])) & (data[:, 0] < max(minmax_list[0])) &
		(data[:, 1] > min(minmax_list[1])) & (data[:, 1] < max(minmax_list[1])) )
	data = data[interp_box]

	# # compare all 3 available interpolation methods
	methods = ['nearest', 'linear', 'cubic']
	destination_data_list = []
	for method in methods:
		destination_values = griddata(data[:, :2], data[:, 2], destination_coordinate_list, method=method)
		destination_data = np.empty(d.shape)
		for j in range(destination_pixel_list.shape[0]):
			xy_pix = destination_pixel_list[j, :]
			x_pix, y_pix = xy_pix[0], xy_pix[1]
			destination_data[x_pix, y_pix] = destination_values[j]
		destination_data_list.append(destination_data)

	h['HISTORY'] = "THIS IS NOT HERSCHEL --- THIS IS PLANCK 857 GHz"
	# # write the nearest neighbor version out to file
	fits.writeto(destination_filenames[i], destination_data_list[0], h, overwrite=True)


	plt.figure()
	plt.subplot(221)
	plt.imshow(d, origin='lower')
	plt.title('true sky')
	plt.subplot(222)
	plt.imshow(destination_data_list[0], origin='lower')
	plt.title(methods[0])
	plt.subplot(223)
	plt.imshow(destination_data_list[1], origin='lower')
	plt.title(methods[1])
	plt.subplot(224)
	plt.imshow(destination_data_list[2], origin='lower')
	plt.title(methods[2])
	plt.show()
