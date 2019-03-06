import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

directory_stub = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"

fns = [
	("PACS160um-", "-remapped-conv.fits"),
	("SPIRE250um-", "-remapped-conv.fits"),
	("SPIRE350um-", "-remapped-conv.fits"),
	("SPIRE500um-", "-remapped-conv.fits"),
]

pacs_offset = "plus046-"
data_types = ("image", "error")
cutout_stub = "-CUTOUT"

def gen_filenames(band_index):
	partA, partB = fns[band_index]
	img, err = data_types
	if "PACS160" in partA:
		img = pacs_offset + img
	img = directory_stub + partA + img + partB
	err = directory_stub + partA + err + partB
	return img, err

def gen_writenames(band_index):
	partA, partB = fns[band_index]
	fits_stub = partB[-5:]
	partB = partB[:-5]
	return tuple(directory_stub+partA+data_type+partB+cutout_stub+fits_stub for data_type in data_types)


xlim = (250, 280)
ylim = (560, 610)
def cutout(array):
	# X and Y are FITS standard
	# Numpy would use "[Y, X]" indexing
	return array[ylim[0]:ylim[1], xlim[0]:xlim[1]]


def cutout_all():
	for i in range(len(fns)):
		both_filenames = gen_filenames(i)
		both_writenames = gen_writenames(i)
		for filename, writename in zip(both_filenames, both_writenames):
			print("READING FROM\n{}".format(filename))
			data, header = fits.getdata(filename, header=True)
			data = cutout(data)
			print("WRITING TO\n{}".format(writename))
			fits.writeto(writename, data, header, overwrite=True)


cutout_all()
