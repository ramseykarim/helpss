import numpy as np


def discretize(data, n_points=None, scale_f=None,
    lims=None, beam_size=1, SNR0=3, SNR1=3, within_pixel=True):
    """
    Convert a pixelated (regularly sampled), N-dimensional density function
    into a representative discrete set of points whose point-to-point density
    indicates the original pixel values.

    DATA: numpy.ndarray of any dimension.
        DATA must already be all valid, real floats above zero. Rescaling can
        be done with SCALE_F. NaNs should already be set to zero (they won't
        be sampled)
    N_POINTS (optional): int number of sample points in returned point set.
        If not set, calculated based on LIMS, BEAM_SIZE, and desired SNR0
    SCALE_F (optional): ufunc/callable(float) rescaling function for DATA
    LIMS (optional): tuple(float low, float high) lowest and highest limits to
        resolve well
    BEAM_SIZE (optional): int number of pixels across a beam
        only matters if you 1) didn't specify N_POINTS or
        2) set WITHIN_PIXEL to False, indicating you want to sample the beam
    SNR0, SNR1 (optional): desired signal-to-noise at beams with LIMS values
    WITHIN_PIXEL (optional): sample within a pixel. If false, sample a beam
        around a pixel
    """
    data_filter = (data > 0)
    minmax = (np.min(data[data_filter]), np.max(data[data_filter]))
    if lims is None:
        lims = minmax
    elif None in lims:
        lims = (x if x is not None else m for x, m in zip(lims, minmax))
    if scale_f is not None:
        lims = tuple(scale_f(l) for l in lims)
        data[data_filter] = scale_f(data[data_filter])
    if n_points is None:
        n_points = optimal_N(np.sum(data), lims[0], SNR0, beam_size)
    cumulative_data = np.cumsum(data.flatten(order='C'))
    cdf = cumulative_data / cumulative_data[-1]
    rand_cdf_values = np.random.uniform(size=(n_points, 1))
    pixel_indices = np.searchsorted(cdf, rand_cdf_values)
    s = None if within_pixel else beam_size
    pixel_centers = get_position(pixel_indices, data.shape)
    pixel_scatters = sample_kernel(size=n_points, sigma=s)
    return pixel_centers + pixel_scatters

def optimal_N(total_weight, smallest_weight, SNR0, beam_size):
    # Calculates a number of points that will provide at least a signal-to-noise
    #  of SNR0 in beams with values of smallest_weight
    return int((total_weight * SNR0**2) / (smallest_weight * beam_size**2))

def get_position(i, shape):
    # Get pixel indices from flattened indices
    # i must be a 1-D array
    return np.concatenate(np.unravel_index(i, shape), axis=1).astype(float)

def sample_kernel(size=1, sigma=None):
    if sigma is None:
        # Within pixel
        return np.random.uniform(size=(size, 2))
    else:
        # Within beam
        return np.random.normal(scale=sigma, size=(size, 2))


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import pickle
    with open("MantiPython/grid_file.pkl", 'rb') as pfl:
        data = pickle.load(pfl)
    data[np.isnan(data)] = 0
    if np.min(data) < 0:
        data += -1*np.min(data)
    print(data.min(), data.max())


    # print("testing discretize")
    # working_dir = "/home/ramsey/Documents/Research/Filaments/"
    # soln = "T4-absdiff-Per1J-plus045-pow-1000-0.1-1.80.fits"
    # from astropy.io import fits
    # img = fits.getdata(working_dir + soln, 3)
    # img[np.isnan(img)] = 0
    # print("N_POINTS: ", end="")
    # print(optimal_N(np.sum(np.sqrt(img)), np.sqrt(1e20), 1., 3.))
    # points = discretize(img, scale_f=np.sqrt, lims=(1e20, None), SNR0=3.,
    #     beam_size=3)
    # print(points.shape)
    # plt.plot(points[:, 1], points[:, 0], ',', alpha=0.01)
    # plt.show()
