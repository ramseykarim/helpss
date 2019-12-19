import numpy as np
import matplotlib.pyplot as plt

# import pipeline functions (written by Ramsey)
from .. import calc_offset
from .. import modify_fits


# Let's say you've got all your data in this directory
# This is where your PACS/SPIRE data is, and it's where your Manticore solutions will go
working_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
# You may want to save the edited files to a new directory
# You can specify a directory rather than a full filename with a string ending in "/"
new_dir = "./"


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
derived_offset = model.get_offset(full_diagnostic=True)
# get_offset() can also be run with the keyword arg "full_diagnostic=True",
#    showing three plots rather than just one. See documentation for details.
# The offset is about 45.5 in the Per1 test example, and I've been rounding down to 45.


# Now that we have the offset, let's add it to the "-remapped-conv.fits" version
pacs_flux_filename = working_dir + "PACS160um-image-remapped-conv.fits"

# Apply the offset, rounded to nearest int, to the convolved PACS image.
# If there was an issue with the automatic derivation, you can write this number in yourself.
calibrated_pacs_flux_filename = modify_fits.add_offset(int(derived_offset), pacs_flux_filename, savename=new_dir)
# These modify_fits functions all return the new filenames (complete with directory path)
# NOTE: these will overwrite calibrated maps of the same filename (should be identical anyway)

# There used to be a step in which we add systematic uncertainty to the error maps.
# We no longer do this.

# That should be it!
