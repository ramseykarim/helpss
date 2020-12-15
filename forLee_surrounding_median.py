"""

Created: October 28, 2020
"""

def replace_medians(array, pct_valid=70):
    """
    Replace isolated NaNs with the median of the 8 surrouding pixels.
    The 8 pixels are all but the center in a 3x3 square centered on
    the NaN pixel.

    :param array: a 2D numpy array representing an image. The values of
        the array should be floats or something, and the array can have NaNs.
    :param pct_valid: the percentage of the 8 surrounding pixels that must be
        non-NaN in order for the pixel to be filled in.
    """
    # Trim off 1 pixel on all four sides so that indexing is easier.
    # We will never want to deal with the edges
    trimmed_array = array[1:-1, 1:-1]

    # Gather the shifted images
    # Shift right
    shifted_arrays.append(array[:, :-1])
    # Shift left
    shifted_arrays.append(array[:, 1:])
    # Shift up
    shifted_arrays.append(array[:-1, :])
    # Shift down
    shifted_arrays.append(array[1:, :])

    # Shift up-right
    shifted_arrays.append(array[:-1, :-1])
    # Shift down-right
    shifted_arrays.append(array[1:, :-1])
    # Shift down-left
    shifted_arrays.append(array[1:, 1:])
    # Shift up-left
    shifted_arrays.append(array[:-1, 1:])

    # Stack arrays
    shifted_arrays = np.array(shifted_arrays)
    # Find percentage of valid pixels
    surrounding_nans = None # TODO LEFT OFF HERE (nanmask?)
