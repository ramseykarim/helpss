import numpy as np
import matplotlib.pyplot as plt

# import pipeline functions (written by Ramsey)
import calc_offset
import modify_fits


# Let's say you've got all your data in this directory
working_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
# This is where your PACS/SPIRE data is, and it's where your Manticore solutions will go

# PACS CALIBRATION

# Tracy's procedures have been run, so there are PACS/SPIRE maps named
#    "-remapped.fits" and "-remapped-conv.fits"
# We will use the "-remapped.fits" PACS map to calculate the necessary calibration offset
# Were we to use "-remapped-conv.fits", we'd simply need to recalculate
#    the beam needed to convolve PACS160 up to GNILC
pacs_flux_filename = working_dir + "PACS160um-image-remapped.fits"
# This file needs correct WCS information; Tracy's procedures handle that just fine

# Create the "GNILCModel" object that Ramsey wrote
model = calc_offset.GNILCModel(pacs_flux_filename)
# The above call does most of the heavy lifting. Some useful things are stored in memory now.
# Use the model to calculate the offset and show some diagnostic figures
model.get_offset()
# get_offset() can also be run with the keyword arg "full_diagnostic=True",
#    showing three plots rather than just one. See documentation for details.
# The offset is about 45.5 in this test example, and I've been rounding down to 45.

# Now that we have the offset, let's add it to the "-remapped-conv.fits" version
pacs_flux_filename = working_dir + "PACS160um-image-remapped-conv.fits"

calibrated_pacs_flux_filename = modify_fits.add_offset(45, pacs_flux_filename)
# These modify_fits functions all return the new filenames (complete with directory path)
# NOTE: these will overwrite calibrated maps of the same filename (should be identical anyway)


# PACS+SPIRE SYSTEMATIC UNCERTAINTY

# The PACS error map needs 5%, and the SPIRE errors need 1.5%
# PACS error map from Tracy's procedures should be called
pacs_error_filename = working_dir + "PACS160um-error-remapped-conv.fits"
# Add 5% error to the PACS error map, with the calibrated map as a reference
modify_fits.add_systematic_error(0.05, pacs_error_filename, calibrated_pacs_flux_filename)

# The SPIRE maps follow the same name convention and all need the same fraction of their flux
#    added as systematic uncertainty
# A simple loop will do
def generate_spire_flux_filename(wavelength_um):
    return working_dir + "SPIRE{:d}um-image-remapped-conv.fits".format(wavelength_um)
def generate_spire_error_filename(wavelength_um):
    return working_dir + "SPIRE{:d}um-error-remapped-conv.fits".format(wavelength_um)
for band in [250, 350, 500]:
    modify_fits.add_systematic_error(0.015, generate_spire_error_filename(band),
        generate_spire_flux_filename(band))

# That should be it!
