from astropy.io import fits

"""
Example for writing out multi-extension FITS files.
I offer 2 options for writing:
    1) reuse the old headers completely (faster)
    2) make new headers (more control)
These methods ASSUME that the image grids are the EXACT SAME.
They will reassign the EXACT SAME WCS info for the writeouts.

Created: Oct 15, 2019
Author: Ramsey Karim
"""

# This can be multi-extension
read_filename = ""
write_filename = ""

# For example
extensions_i_care_about = [3, 5, 7]
ext_data = []
ext_hdrs = []
with fits.open(read_filename) as hdul:
    primary_header = hdul[0].header
    for i in extensions_i_care_about:
        ext_data.append(hdul[i].data)
        ext_hdrs.append(hdul[i].header)

def operation(array):
    return array * 2
# Operate on data as numpy arrays
# Could be masking and combining multiple images
new_data = [operation(data) for data in ext_data]


#################################################
# Option 1: pretty quick multi-extension writout
# Reuses headers, so assumes the same frames (not necessarily same order)
#################################################
# Make list of relevant image HDUs using their original headers.
# These headers have the WCS info in them already
hdu_list = [fits.ImageHDU(d, header=h) for d, h in zip(new_data, ext_hdrs)]
# Prepend the primary HDU from the original file
# (contains no data, just global header)
primary_hdu = fits.PrimaryHDU(header=primary_header)
hdu_list = [primary_hdu] + hdu_list
# Turn the HDU list into a fits.HDUList object
hdul = fits.HDUList(hdu_list)
# Write it out using HDUList writeto method
# (if overwrite=False, it will throw an error if write_filename exists)
hdul.writeto(write_filename, overwrite=True)
# This method doesn't offer much flexibility with the headers,
# but is fairly quick.


#################################################
# Option 2: more code but more flexible multi-extension writout
# Create new headers
#################################################
# Create a fresh primary HDU and a fresh image header template
primary_hdu = fits.PrimaryHDU()
img_hdr = fits.Header()
# Need to import WCS
from astropy.wcs import WCS
# Get WCS from old headers (they all should have the same WCS)
wcs_info = WCS(ext_hdrs[0])
# Assign this WCS to the image header template
img_hdr.update(wcs_info.to_header())
# Add some other information if you'd like
from datetime import datetime, timezone
img_hdr['COMMENT'] = "CREATED AS EXAMPLE"
img_hdr['CREATOR'] = ("Ramsey: {}".format(str(__file__)), "FITS file creator")
img_hdr['OBJECT'] = ("my_object_name", "Target name")
img_hdr['DATE'] = (datetime.now(timezone.utc).astimezone().isoformat(),
    "File creation date")
# Create list of HDUs using this template header
hdu_list = [fits.ImageHDU(d, header=img_hdr) for d in new_data]
# Create extension labels describing the content of each frame
ext_labels = ['x', 'x^2', 'sin(x)']
# Create unit labels for each field
ext_units = ['cm', 'cm2', '']
# Assign labels and units to each field
for label, unit, hdu in zip(ext_labels, ext_units, hdu_list):
    hdu.header['EXTNAME'] = label
    hdu.header['BUNIT'] = (unit, "Data unit")
# Prepend primary HDU
hdu_list = [primary_hdu] + hdu_list
# Create HDUList object
hdul = fits.HDUList(hdu_list)
# Write out
hdul.writeto(write_filename)
