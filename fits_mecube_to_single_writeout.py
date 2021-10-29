#!/n/astromake/opt/python/anaconda3/bin/python3
"""
Function and script to open a multi-extension FITS file, extract a single plane, and
re-save it as its own single-plane FITS image.

Created: October 27, 2021
"""
__author__ = "Ramsey Karim"


import os
import argparse

from astropy.io import fits


def super_simple_example():
    """
    Much simpler script which can be copied into another block of code.
    """
    source_filename = "/n/sgraraid/filaments/Perseus/Herschel/results/full-1.5.3-Per1-pow-750-0.05625-1.6.fits"
    target_filename = "/home/rkarim/Desktop/test.fits"
    extension = 'T' # Can be key or integer extension index
    # Open the file and grab the extension
    with fits.open(source_filename) as hdul:
        hdu = hdul[extension]
        # Convert from ImageHDU to Primary HDU so the saved file is really 2D
        hdu = fits.PrimaryHDU(data=hdu.data, header=hdu.header)
        # Save the extension
        hdu.writeto(target_filename, overwrite=True)


"""
I also turned this file into a script to do this, which you can run from command line!
"""


def extract_and_save_extension(source_filename, target_filename, extension=None):
    """
    Function to open a (presumed) multi-extension FITS cube,
    extract a single plane, and re-save it under a new filename.
    The union of the Primary Header and the extension Header will be saved
    to the new file (exension Header will be given priority).
    :param source_filename: path to the existing FITS cube to be read.
        The path needs to be valid, but can be absolute or relative to the
        working directory
    :param target_filename: path to the FITS file to be written.
        Same as above, the path needs to be valid but can be relative or
        absolute.
        If a file already exists at that path, THIS FUNCTION WILL OVERWRITE IT!!
        It will warn you with a message, but it will still overwrite it.
    :param extension: the extension number or key to be extracted.
        Default is the first extension to have valid data.
    """
    # Load HDUList
    # Validate source path implicitly, an error will be thrown if not found
    hdul = fits.open(source_filename)
    hdu = None
    if extension is not None:
        # Grab the requested extension
        hdu = hdul[extension]
    else:
        # Grab the first extension with data
        for hdu in hdul:
            if hdu.data:
                break
    # Get the Primary HDU Header
    hdr = hdul[0].header
    # Get the extension Header
    hdr.update(hdu.header)
    # Make a new HDU with the data and new header
    new_hdu = fits.PrimaryHDU(data=hdu.data, header=hdr)
    # Close the HDUList
    hdul.close()
    # Check if target filename is occupied
    if os.path.exists(target_filename):
        print(f"File already exists at {target_filename}; overwriting...")
    # Write the new HDU to the target path
    # Iimplicitly validates the path, will throw an error if it doesn't exist
    new_hdu.writeto(target_filename, overwrite=True)


def collect_arguments():
    """
    Collect arguments to the extraction function from the command line.
    :returns: dict() with keys source_filename, target_filename, extension
    """
    parser = argparse.ArgumentParser(description="Command line script to extract a single extension from a multi-extension cube and save it.")
    parser.add_argument('source_filename', type=str, help='Multi-extension FITS file to read from.')
    parser.add_argument('target_filename', type=str, help='Filename to write to. If file exists at that location already, it will be overwritten.')
    parser.add_argument('--extension', type=str, help='Extension key or index. If numeric (positive integer), interpreted as index. If the key has special characters (like the parentheses in N(H2)), use quotes. Default is the first extension with valid data.')
    args = parser.parse_args()
    if args.extension and args.extension.isnumeric():
        args.extension = int(args.extension)
    return vars(args)


# Wrapper script so this can be run from command line.
if __name__ == "__main__":
    extract_and_save_extension(**collect_arguments())

