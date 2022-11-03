"""

Created: October 28, 2020

Updated: September 21, 2022
For Rachel, who needs to know how  many True/valid neighbors each pixel has
"""
import numpy as np
import matplotlib.pyplot as plt

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
    shifted_arrays = []
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

    count = shifted_arrays.astype(int).sum(axis=0)


    # Find percentage of valid pixels
    surrounding_nans = None # TODO LEFT OFF HERE (nanmask?)



def better_pixel_shift(original_array):
    """
    September 21, 2022
    not a complete function, but a guide on how to count neighbors
    Assume array is Boolean (or 1/0 int or float)
    """
    # Pad with 0s
    array = np.pad(original_array, 1)

    shifted_arrays = []
    # Shift right (i.e. check left neighbors)
    shifted_arrays.append(array[1:-1, :-2])
    # Shift left
    shifted_arrays.append(array[1:-1, 2:])
    # Shift up
    shifted_arrays.append(array[:-2, 1:-1])
    # Shift down
    shifted_arrays.append(array[2:, 1:-1])

    # Shift up-right
    shifted_arrays.append(array[:-2, :-2])
    # Shift down-right
    shifted_arrays.append(array[2:, :-2])
    # Shift down-left
    shifted_arrays.append(array[2:, 2:])
    # Shift up-left
    shifted_arrays.append(array[:-2, 2:])

    shifted_arrays_cube = np.array(shifted_arrays)
    neighbor_count = shifted_arrays_cube[:, 1:-1, 1:-1].sum(axis=0)

    plt.subplot(121)
    plt.imshow(original_array, origin='lower')
    plt.subplot(122)
    plt.imshow(neighbor_count, origin='lower')
    plt.show()


if __name__ == "__main__":
    arr = np.zeros((30, 30))
    arr[10:13, 20:23] = 1
    better_pixel_shift(arr)
