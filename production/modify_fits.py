from astropy.io import fits


def add_offset(offset, filename, savename=None, extension=0):
    """
    Add a fixed offset to a FITS image. Will only save a single extension.
    Does not overwrite the original (unless savename==filename)
    :param offset: int or float additive offset
    :param filename: file path to FITS to modify
    :param savename: file path to which to write the modified data
        Default is to add "-plus######" right before ".fits"
        Will overwrite anything already saved as savename
    :param extension: FITS extension of data to read; default is 0
    """
    # Make tidy offset string based on float vs int
    offset_str = "{:06d}".format(offset) if type(offset) == int else "{:06.1f}".format(offset)
    if savename is None:
        # Add "-plus######.fits" for a default
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


def add_systematic_error(flux_fraction, error_filename, flux_filename,
                         error_extension=0, flux_extension=0, savename=None):
    """
    Adds the specified fraction of the flux map to the error map as
        uncorrelated systematic uncertainty.
    Specify the flux fraction as a decimal fraction (NOT percent)
    Does not overwrite the existing error map unless you specify
        the original error map's name as the savename here.
    :param flux_fraction: decimal fraction of flux to be added to error
        NOT PERCENTAGE. 1.5% should be input here as 0.015
    :param error_filename: file path to error map to modify.
        Does not overwrite this error map (unless savename==error_filename)
    :param flux_filename: file path to flux map referenced by flux_fraction
        Flux and error need to be on grids of the same shape
        They should also be aligned in WCS, but this function will not check
    :param error_extension: FITS extension where error map is found. default=0
    :param flux_extension: FITS extension where flux map is found. default=0
    :param savename: file path to which to write the modified error map
        Default is to add "-plus#.#pct" right before ".fits"
    """
    # Make flux percentage string
    pct_string = "{:03.1f}pct".format(flux_fraction*100)
    if savename is None:
        # Add "-plus###pct.fits" for a default
        filename_first, fits_stub = error_filename[:-5], error_filename[-5:]
        assert fits_stub == ".fits"
        savename = f"{filename_first}-plus{pct_string}{fits_stub}"
    with fits.open(error_filename) as hdul:
        error = hdul[error_extension].data
        head = hdul[error_extension].header
    with fits.open(flux_filename) as hdul:
        flux = hdul[flux_extension].data
    error += flux * flux_fraction
    head['COMMENT'] = "Added {:s} % flux as uncorrelated systematic uncertainty".format(pct_string)
    try:
        fits.writeto(savename, error, header=head)
        print(f"Wrote to {savename} with {pct_string} offset")
    except OSError:
        fits.writeto(savename, error, header=head, overwrite=True)
        print(f"Wrote to {savename} with {pct_string} offset (overwriting existing)")
