"""
Mask function designed from specifications by LGM in the Sept 9, 2020 HELPSS
meeting.

The function takes an array with 3 possible values: 1, 0, and NaN

The function also takes in a first_pad and a second_pad value, both specified
in number of iterations (mostly equivalent to pixels).
The function will dilate the 1 pixels <first_pad> iterations into 0 pixels
from the boundaries of 1 pixels (but not NaN pixels).
Then, it will dilate <second_pad> iterations into the 0 pixels from that new
boundary. That second round of dilation will be represented by 1s in a new mask
which will be returned. In the new mask, all other pixels will be 0s, except
for the NaNs which will remain NaN.

Created: September 15, 2020
"""
__author__ = "Ramsey Karim"

import numpy as np


def padded_mask(array, first_pad, second_pad):
    """
    Masking function that selects a border around designated pixels.

    The argument array should contain 1s, 0s, and NaNs. The 1s represent
    pixels FROM which we dilate outwards. The 0s represent pixels INTO which
    we dilate. The NaNs represent pixels we should never bother with; those NaNs
    will be conserved in the resulting map, and no other NaNs will be created.

    The border between 1s and 0s is expanded into the 0s, pixel by pixel, for
    <first_pad> iterations. From there, it is expanded again in the same way
    for <second_pad> iterations. All pixels that were added in the second set
    of iterations are assigned 1 in the resulting mask array. All other pixels
    are assigned 0, except for NaNs, which remain NaN.

    :param array: np.ndarray of float
        The input mask, represented as a 2D array of floats containing only
        1s, 0s, and NaNs
    :param first_pad: int
        number of iterations (roughly equivalent to pixels) for the first
        round of dilation
    :param second_pad: int
        number of iterations for the second round of dilation
    :return: np.ndarray of float
        Same format as the input array; 1s represent included in mask, 0s
        represent excluded, and NaNs are conserved from the input.
        The 1s represent regions between the boundaries of the first and second
        dilations.
    """
    # Split the array into two boolean arrays
    # Make a NaN mask, 1 where NaN
    nan_mask = np.isnan(array)
    # Make a mask where array 1s are 1 and everything else (0 and NaN) is 0
    first_round_mask = (array == 1)
    # Dilate <first_pad> times
    for i in range(first_pad):
        first_round_mask = dilate(first_round_mask, nan_mask)
    # Save this state of the mask
    second_round_mask = np.copy(first_round_mask)
    # Dilate <second_pad> times
    for i in range(second_pad):
        second_round_mask = dilate(second_round_mask, nan_mask)
    # Compare the first and second round dilations
    result = (second_round_mask & ~first_round_mask).astype(float)
    # Conserve NaNs
    result[nan_mask] = np.nan
    return result


def dilate(mask, forbidden):
    """
    Dilate the 1s in mask into the 0s in mask. Never touch pixels that are
    1 in forbidden. Dilation is a one-pixel shift into all of the 8 adjacent
    available pixels.

    Helper function for padded_mask; this is NOT the main function.

    :param mask: np.ndarray of bool
        1 where we dilate FROM, 0 where we dilate INTO
    :param forbidden: np.ndarray of bool
        1 where we should never work (NaN equivalent)
        All 1s in mask should be 0 in forbidden
    """
    new_mask = np.copy(mask)
    # Shift right
    new_mask[:, 1:] |= mask[:, :-1]
    # Shift left
    new_mask[:, :-1] |= mask[:, 1:]
    # Shift up
    new_mask[1:, :] |= mask[:-1, :]
    # Shift down
    new_mask[:-1, :] |= mask[1:, :]

    # Shift up-right
    new_mask[1:, 1:] |= mask[:-1, :-1]
    # Shift down-right
    new_mask[:-1, 1:] |= mask[1:, :-1]
    # Shift down-left
    new_mask[:-1, :-1] |= mask[1:, 1:]
    # Shift up-left
    new_mask[1:, :-1] |= mask[:-1, 1:]

    new_mask[forbidden] = False

    return new_mask
