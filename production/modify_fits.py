from astropy.io import fits


def add_offset(offset, filename, savename=None, extension=0):
    """
    Add a fixed offset to a FITS image. Will only save a single extension.
    Does not overwrite the original (unless savename==filename)
    :param offset: int or float additive offset
    :param filename: file path to FITS to modify
    :param savename: file path to which to write the modified data
        Default is to add "-plus#####" right before ".fits"
        Will overwrite anything already saved as savename
    :param extension: FITS extension of data to read; default is 0
    """
    # Make tidy offset string based on float vs int
    offset_str = "{:05d}".format(offset) if type(offset) == int else "{:05.1f}".format(offset)
    if savename is None:
        # Add "-plus#####.fits" for a default
        filename_first, fits_stub = filename[:-5], filename[-5:]
        assert fits_stub == ".fits"
        savename = f"{filename_first}-plus{offset_str}{fits_stub}"
    with fits.open(filename) as hdul:
        data = hdul[extension].data
        head = hdul[extension].header
    data += offset
    head['COMMENT'] = "Added {:s} MJy/sr offset".format(offset_str)
    # Overwrite if previous version exists, but inform the user
    try:
        fits.writeto(savename, data, header=head)
        print(f"Wrote to {savename} with {offset_str} offset")
    except OSError:
        fits.writeto(savename, data, header=head, overwrite=True)
        print(f"Wrote to {savename} with {offset_str} offset (overwriting existing)")



