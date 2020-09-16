"""
Mask function designed from specifications by LGM in the Sept 9, 2020 HELPSS
meeting.

The function takes in one 2D array, with float values.
The array represents a mask with three possible values.
The value 1 should represent valid pixels, "included" in some sense.
The value 0 should represent invalid pixels along whose boundaries with valid
pixels this function will operate.
The value NaN should should represent invalid pixels whose boundaries with
valid pixels are ignored.

The function also takes in a first_pad and a second_pad value, both specified
in number of iterations (mostly equivalent to pixels).
The function will dilate the 0 pixels <first_pad> iterations into 1 pixels
from the boundaries of 0 pixels (but not NaN pixels).
Then, it will dilate <second_pad> iterations into the 1 pixels from that new
boundary. That second round of dilation will be represented by 1s in a new mask
which will be returned. In the new mask, all other pixels will be 0s, except
for the NaNs which will remain NaN.

Created: September 15, 2020
"""
__author__ = "Ramsey Karim"

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits


def dilate(array, times=1, forbidden=None):
    """
    This is the old dilation routine from boolean_islands.py
    I want to benchmark this against a new routine I write.

    I modified it to make a "forbidden" set out of a bool array where 1s are
    forbidden
    """
    if times == 0:
        return np.copy(array)
    if forbidden is None:
        forbidden = set()
    else:
        forbidden = set(zip(*np.where(forbidden)))
    shape = array.shape
    def valid_point(i, j):
        within_bounds = not ((i < 0) or (j < 0) or (i >= shape[0]) or (j >= shape[1]))
        not_forbidden = (i, j) not in forbidden
        return within_bounds and not_forbidden
    dilating_ones = set(zip(*np.where(array)))
    new_ones = set()
    for iteration in range(times):
        totally_new_ones = set()
        for i0, j0 in dilating_ones:
            totally_new_ones.update(*(set((i0 + di, j0 + dj) for dj in range(-1, 2) if (di or dj) and valid_point(i0 + di, j0 + dj)) for di in range(-1, 2)))
        totally_new_ones.difference_update(dilating_ones)
        new_ones |= totally_new_ones
        dilating_ones = totally_new_ones
    new_array = np.copy(array).astype(bool)
    new_array[tuple(zip(*new_ones))] = True
    return new_array


def new_dilate(array, forbidden):
    """
    Both arguments should be boolean
    :param array: 1 where valid, 0 where invalid
    :param forbidden: 1 where we should never work (NaN equivalent)
        all array 1s should be 0 in forbidden
    """
    new_array = np.copy(array)
    # Shift right
    new_array[:, 1:] |= array[:, :-1]
    # Shift left
    new_array[:, :-1] |= array[:, 1:]
    # Shift up
    new_array[1:, :] |= array[:-1, :]
    # Shift down
    new_array[:-1, :] |= array[1:, :]

    # Shift up-right
    new_array[1:, 1:] |= array[:-1, :-1]
    # Shift down-right
    new_array[:-1, 1:] |= array[1:, :-1]
    # Shift down-left
    new_array[:-1, :-1] |= array[1:, 1:]
    # Shift up-left
    new_array[1:, :-1] |= array[:-1, 1:]

    new_array[forbidden] = False

    return new_array


def mask_from_old_dilation_routine(array, first_pad, second_pad):
    """
    Use the old dilation routine to make the mask as specified in this file's
    header.
    """
    nan_mask = np.isnan(array)
    valid_mask = (array == 1) # array == 0
    first_round_mask = dilate(valid_mask, times=first_pad, forbidden=nan_mask)
    second_round_mask = dilate(first_round_mask, times=second_pad, forbidden=nan_mask)
    result = (second_round_mask & ~first_round_mask).astype(float)
    result[nan_mask] = np.nan
    return result


def mask_from_new_dilation_routine(array, first_pad, second_pad):
    """
    Use the new array indexing shift dilation routine to make the mask as
    specified in this file's header.
    """
    nan_mask = np.isnan(array)
    first_round_mask = (array == 1) # array == 0
    for i in range(first_pad):
        first_round_mask = new_dilate(first_round_mask, nan_mask)
    second_round_mask = np.copy(first_round_mask)
    for i in range(second_pad):
        second_round_mask = new_dilate(second_round_mask, nan_mask)
    result = (second_round_mask & ~first_round_mask).astype(float)
    result[nan_mask] = np.nan
    return result


if __name__ == "__main__":
    arr = fits.getdata("../syp_mask_3layer.fits")
    # newmask = mask_from_old_dilation_routine(arr, 5, 20)
    newmask = mask_from_new_dilation_routine(arr, 5, 20)

    """
    The new method is vasly faster. Array shifting!
    """

    fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)
    axes[0].imshow(arr, origin='lower')
    axes[1].imshow(newmask, origin='lower')
    plt.show()
