import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, FK5, Angle
from astropy import units as u
import sys
from collections import deque
try:
    import healpy as hp
except ModuleNotFoundError:
    hp = None
from scipy.interpolate import interpn, griddata#, RectBivariateSpline # (could also work!)


FN_T = "COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.00.fits"
FN_tau = "COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.00.fits"
FN_beta = "COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.00.fits"


def get_galactic_from_WCS(w, ij_list):
    # ij_list should be an array of shape (n x m, 2), for an original image of dimensions (n, m)
    # ij_list should be calculated from the numpy array
    # Convert pixel X/Y list to RA/DEC coordinates
    coord_list = w.wcs_pix2world(ij_list, 0)
    ra_list, dec_list = coord_list[:, 0], coord_list[:, 1]
    # Convert RA/DEC to galactic lon/lat
    coord_list = SkyCoord(Angle(ra_list, unit=u.deg), dec_list*u.deg, frame=FK5)
    l, b = coord_list.galactic.l.degree, coord_list.galactic.b.degree
    return l, b


def get_ij_list_slower(array_shape, mask=None):
    # Given the array.shape, create an array of pixel coordinates
    # Resulting array will be shape (n x m, 2) for a 2D array of shape (n, m)
    # NOTE: could also do this with meshgrid and flatten
    # Mask should be a boolean array of same dimensions as the image
    # Make separate x/y pixel index arrays
    pix_arrays = tuple(np.arange(x) for x in array_shape)
    # Configure mask (probably a nan mask)
    # Mask == True means INCLUDE this. False means DO NOT INCLUDE
    if mask is None:
        mask_fn = lambda yp, xp: True
    else:
        mask_fn = lambda yp, xp: mask[yp, xp]
    # Make list for every possible x/y pixel index pair
    pix_list = deque()
    for yp in pix_arrays[0]:
        for xp in pix_arrays[1]:
            # Skipping nans
            if mask_fn(yp, xp):
                # Note the reversed order of X/Y
                pix_list.append(np.array([xp, yp]))
    pix_list = np.array(pix_list)
    return pix_list


def get_ij_list(array_shape, mask=None):
    # Given the array.shape, create a list of pixel coordinates
    # Resulting array will be shape (n x m, 2) for a 2D array of shape (n, m)
    #   AND coordinates will be REVERSED
    # More pythonic version of the above method.
    # Mask should be a boolean array of same dimensions as the image
    if mask is None:
        mask = np.full(array_shape, True)
    # Line by line:
    # pix_arrays = tuple(np.arange(x) for x in array_shape)
    # xxyy = np.meshgrid(*pix_arrays)
    # xym = (x.ravel() for x in (*xxyy, mask))
    # return np.array(tuple((y, x) for x, y, m in zip(*xym) if m))
    # Just to say I can do it:
    return np.array(tuple((y, x) for x, y, m in zip(*(x.ravel() for x in (*np.meshgrid(*(np.arange(x) for x in array_shape), indexing='ij'), mask))) if m))


def assign_interpolated_values(interpolated_vals_list, pixel_list, target_shape):
    # Given a list of interpolated values, a list of coordinates, and the shape of the final image,
    #  generate that final image
    # interpolated_vals_list should be len(n x m)
    # pixel_list should be shape (n x m, 2)
    # target_shape should thus be (n, m)
    # The lists of interpolated values and pixel coordinates should match up
    # The target_shape should match the pixel coordinates in the pixel_list
    # Assumes the pixel coordinates are swapped (j, i) in the list (i->n, j->m)
    interpolated_data = np.full(target_shape, np.nan)
    for i in range(pixel_list.shape[0]):
        xy_pix = pixel_list[i, :]
        # Again, note the X/Y order is flipped to revert the first flip
        interpolated_data[xy_pix[1], xy_pix[0]] = interpolated_vals_list[i]
    return interpolated_data


def regrid_healpix_to_fits(source_hp, dest_fits_data, dest_fits_header, nest=False):
    """
    Designed for Planck Archive HEALPix to Herschel PACS/SPIRE grids
    Inputs are:
        HEALPix map, already opened with healpy
        Herschel data array and fits header

    Set "nest" parameter to True if it's the frequency map
    """
    pix_list = get_ij_list(dest_fits_data.shape, mask=(~np.isnan(dest_fits_data)))
    ##### BUG pix_list may be empty!

    # Convert this pixel X/Y list to galactic lon/lat coordinates
    l, b = get_galactic_from_WCS(WCS(dest_fits_header), pix_list)
    # Interpolate from HEALPix map to these galactic coordinates
    interpolated_vals_list = hp.get_interp_val(source_hp, l, b, lonlat=True, nest=nest)
    # Generate the final array (same shape as PACS/SPIRE) and populate
    interpolated_data = assign_interpolated_values(interpolated_vals_list, pix_list, dest_fits_data.shape)
    return interpolated_data


def get_galactic_limits(target_fits_data, target_fits_header):
    # Returns tuple((min_l, max_l), (min_b, max_b))
    l, b = get_galactic_from_WCS(WCS(target_fits_header), get_ij_list(target_fits_data.shape))
    return (np.min(l), np.max(l)), (np.min(b), np.max(b))


def get_pixel_count(coord_limits_galactic, pixel_scale_as):
    # Given box limits in gal (l, b) degrees and desired pixel scale in arcseconds,
    #  calculate the necessary pixel count.
    # First number is for l direction, second is for b
    # Calculate distances in degrees, convert to arcseconds, divide by pixel scale
    return tuple(int(np.ceil(3600*(tup[1] - tup[0])/pixel_scale_as)) for tup in coord_limits_galactic)


def get_lb_grid_arrays(projection):
    # Given a healpy CartesianProj object, get the defining galactic l and b arrays (like meshgrid inputs)
    # The projection appears to work like meshgrid, so this is fine.
    # I compared using meshgrid. This is correct.
    pix_l, pix_b = projection.arrayinfo['xsize'], projection.arrayinfo['ysize']
    # THERE IS SOMETHING WRONG WITH HEALPY
    # probably because there is weirdness with l>180 ranges
    # healpy uses negative versions of l, which is reflected in get_extent()
    # BUT if you use get_center(), that only worked on my ~160 l. NOT on my ~270 l,
    #  where I got the complement (~90)
    # I think this fools healpy into thinking it's generating invalid coordinates, when queried with i,j
    # So, my hacky fix is to just unmask the arrays even though healpy thinks they're junk
    # This may cause problems somewhere, but it works for now.
    l_range = projection.ij2xy(i=np.full(pix_l, pix_b//2), j=np.arange(pix_l))[0].data
    b_range = projection.ij2xy(i=np.arange(pix_b), j=np.full(pix_b, pix_l//2))[1].data
    l_range = -l_range  # FLIP SIGN ON LONGITUDE (healpy records l as negative?? needs fixing on their end)
    return l_range, b_range


class ProjectionWrapper:

    def __init__(self, target_fits_data, target_fits_header, pixel_scale_as=75):
        self.target_fits_data = target_fits_data
        self.target_fits_header = target_fits_header
        self._projection = self.prepare_projection(target_fits_data, target_fits_header, pixel_scale_as)

    def switch_target(self, new_target_fits_data, new_target_fits_header):
        # Make sure this is completely contained within the CartesianProj limits
        # If it isn't, just make a new object
        self.target_fits_data = new_target_fits_data
        self.target_fits_header = new_target_fits_header

    def prepare_projection(self, target_fits_data, target_fits_header, pixel_scale_as):
        box_l_lim, box_b_lim = get_galactic_limits(target_fits_data, target_fits_header)
        pix_l, pix_b = get_pixel_count((box_l_lim, box_b_lim), pixel_scale_as)
        projection = hp.projector.CartesianProj(xsize=pix_l, ysize=pix_b, lonra=np.array(box_l_lim), latra=np.array(box_b_lim))
        projection.set_flip('astro')
        return projection

    def generate_image(self, healpix_map, nest=False):
        n_side = hp.get_nside(healpix_map)
        vec2pix_func = lambda x, y, z: hp.pixelfunc.vec2pix(n_side, x, y, z, nest=nest)
        return self._projection.projmap(healpix_map, vec2pix_func)

    def interpolate_to_target(self, image, method='nearest'):
        # Interpolates the input image to the target FITS image grid (the one used in __init__)
        # Assumes the input image is on the CartesianProj grid
        # Get information about the input image grid, from which we will interpolate (treat as regular grid)
        l_range, b_range = get_lb_grid_arrays(self._projection)
        # Get information about the target image pixel locations (treat as unstructured points)
        pix_list = get_ij_list(self.target_fits_data.shape, mask=(~np.isnan(self.target_fits_data)))
        l_target, b_target = get_galactic_from_WCS(WCS(self.target_fits_header), pix_list)
        # Note inversion of l, b: l is the "x" coordinate, which takes the "j" index at axis=1
        bl_target_pairs = np.stack((b_target, l_target), axis=1)
        del l_target, b_target
        # Interpolate. Method options are ("linear", "nearest", "splinef2d")
        # Note that we invert the order of the longitudes to make them strictly increasing
        interpolated_vals_list = interpn((b_range, l_range[::-1]), image[:, ::-1], bl_target_pairs, method=method)
        interpolated_data = assign_interpolated_values(interpolated_vals_list, pix_list, self.target_fits_data.shape)
        return interpolated_data


def constrain_coordinates(coordinate_list, pixel_list, limits):
    # filter coordinates to be between RA and DEC limits
    # coordinate array is (N, 2) in RA, DEC order
    # limits are ((ra_lo, ra_hi), (dec_lo, dec_hi))
    # ------
    # get limits
    ra_lo, ra_hi = limits[0]
    dec_lo, dec_hi = limits[1]
    # split RA and DEC
    ra_list, dec_list = coordinate_list[:, 0], coordinate_list[:, 1]
    # make target mask
    ra_mask = (ra_list >= ra_lo) & (ra_list <= ra_hi)
    dec_mask = (dec_list >= dec_lo) & (dec_list <= dec_hi)
    coord_mask = ra_mask & dec_mask
    # apply mask
    ra_list = ra_list[coord_mask]
    dec_list = dec_list[coord_mask]
    # apply mask to pixel list too
    pix_i, pix_j = pixel_list[:, 0], pixel_list[:, 1]
    pix_i, pix_j = pix_i[coord_mask], pix_j[coord_mask]
    return np.stack((ra_list, dec_list), axis=1), np.stack((pix_i, pix_j), axis=1)


def regrid_to_reference(target_array, target_pix2world, source_array, source_pix2world,
    limits=None):
    # general-use regrid function for standard FITS files
    # pix2world should take in 0-indexed X, Y (col j, row i) and return RA, DEC
    # regrids the source values to the target grid
    # assumes that you've selected overlapping regions
    # if limits is not None, should be ((ra_lo, ra_hi), (dec_lo, dec_hi))
    # ------
    # get pixel lists for both target and source
    t_pix_list = get_ij_list(target_array.shape, mask=~np.isnan(target_array))
    s_pix_list = get_ij_list(source_array.shape, mask=~np.isnan(source_array))
    # we did some dark magic with the coordinate ordering so it has to be like this
    # coordinate list for target
    t_coord_list = target_pix2world(t_pix_list)
    # coordinate list for source
    s_coord_list = source_pix2world(s_pix_list)
    # clean target coordinates if limits exist
    if limits is not None:
        t_coord_list, t_pix_list = constrain_coordinates(t_coord_list, t_pix_list, limits)
        s_coord_list, s_pix_list = constrain_coordinates(s_coord_list, s_pix_list, limits)
        # find cropped target grid
        min_i, min_j = np.min(t_pix_list[:, 1]), np.min(t_pix_list[:, 0])
        max_i, max_j = np.max(t_pix_list[:, 1]), np.max(t_pix_list[:, 0])
    s_value_list = source_array[s_pix_list[:, 1], s_pix_list[:, 0]]
    # (though this is also cool) vv (does the same as above but with zip/map)
    # s_value_list = source_array[tuple(zip(*map(tuple, s_pix_list)))]
    # interpolate!
    interpolated_vals_list = griddata(s_coord_list, s_value_list, t_coord_list, method='linear')
    if limits is not None:
        t_value_list = target_array[t_pix_list[:, 1], t_pix_list[:, 0]]
        # grid the data
        t_pix_list[:, 1] -= min_i
        t_pix_list[:, 0] -= min_j
        interpolated_data = assign_interpolated_values(interpolated_vals_list,
            t_pix_list, (max_i - min_i + 1, max_j - min_j + 1))
        # return the smaller target section too
        cropped_target_data = assign_interpolated_values(t_value_list,
            t_pix_list, (max_i - min_i + 1, max_j - min_j + 1))
        return interpolated_data, cropped_target_data
    else:
        # grid the data
        interpolated_data = assign_interpolated_values(interpolated_vals_list,
            t_pix_list, target_array.shape)
        return interpolated_data

def open_healpix(filename, nest=False):
    # Shortcut to opening a HEALPix map, so healpy need not be imported for a single call
    # Nest is TRUE for the frequency maps and FALSE for the dust component maps
    return hp.read_map(filename, nest=nest)


def project_healpix_to_fits(source_hp, target_fits_data, target_fits_header, nest=False, extraction_pixel_scale_as=75):
    # Regrids a HEALPix map to a standard target FITS image
    # Steps through an intermediate projection stage, extracting data from HEALPix at a particular pixel scale
    # Uses linear interpolation
    projection_wizard = ProjectionWrapper(target_fits_data, target_fits_header, pixel_scale_as=extraction_pixel_scale_as)
    intermediate_image = projection_wizard.generate_image(source_hp, nest=nest)
    final_image = projection_wizard.interpolate_to_target(intermediate_image)
    return final_image


# |  projmap(self, map, vec2pix_func, rot=None, coord=None)
# |      Create an array containing the projection of the map.
# |
# |      Input:
# |        - vec2pix_func: a function taking theta,phi and returning pixel number
# |        - map: an array containing the spherical map to project,
# |               the pixelisation is described by vec2pix_func
# |      Return:
# |        - a 2D array with the projection of the map.
# |
# |      Note: the Projector must contain information on the array.
# get_center(self, lonlat=False)

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    planck_857_fn = "/n/sgraraid/filaments/data/HFI_SkyMap_857_2048_R2.02_full.fits"
    spire_fn = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/SPIRE350um-image.fits"
    save_fn = "/n/sgraraid/filaments/data/HFI_SkyMap_857_resSPIRE350_intStep.fits"
    p_map = hp.read_map(planck_857_fn, nest=True)
    with fits.open(spire_fn) as hdul:
        fits_data = hdul[0].data
        fits_header = hdul[0].header
    arr = project_healpix_to_fits(p_map, fits_data, fits_header, nest=True)
    plt.imshow(arr, origin='lower')
    plt.show()
    fits_header['HISTORY'] = "this is Planck data, intermediate-step interpolation method -- Ramsey (1/22/2019)"
    print("use fits.writeto(save_fn, arr, fits_header) to save")
