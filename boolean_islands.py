import numpy as np
if __name__ == "__main__":
    import matplotlib.pyplot as plt
from astropy.io.fits import getdata


"""
Author: Ramsey Karim

Imagine you have a 2D boolean array.
This array may have islands of connected 1s in it.

We want to create a mask of 1s based on these islands.
Specifically, the larger ones; small (~10 pixel) islands are probably noise.

"""

def hash_coord(i, j):
    # i should be the ROW number, j the COLUMN number
    # In other words, j should be the faster-switching index.
    # In numpy, array should be indexed array[i, j]
    return tuple(i, j)

class BooleanIslands:
    def __init__(self, array, neighbor_limit=1):
        self.set_neighbor_limit(neighbor_limit)
        self.unsorted_ones = set(zip(*np.where(array)))
        self.sorted_ones_sets = []
        self.sorted_ones_dict = {x:-1 for x in self.unsorted_ones}
        self._current_set_index = None
        self.find_islands()

    def find_islands(self):
        traversal_stack = set()
        traversal_holding_area = set()
        while self.unsorted_ones:
            traversal_stack.add(self.unsorted_ones.pop())
            self._current_set_index = len(self.sorted_ones_sets)
            while traversal_stack:
                current_point = traversal_stack.pop()
                # Good debug:
                # assert current_point in self.sorted_ones_sets[self.sorted_ones_dict[current_point]]
                i0, j0 = current_point
                adjacent_pairs = sum(tuple(tuple((i0 + di, j0 + dj) for dj in range(-1, 2) if (di or dj)) for di in range(-1, 2)), ())
                adjacency_count = 0
                for adjacent_point in adjacent_pairs:
                    if adjacent_point in self.sorted_ones_dict:
                        adjacency_count += 1
                        adjacent_index = self.sorted_ones_dict[adjacent_point]
                        ### If we want an immediate dilation
                        # if len(self.sorted_ones_sets) > self._current_set_index:
                        #     self.sorted_ones_sets[self._current_set_index].add(adjacent_point)
                        if adjacent_index == self._current_set_index:
                            # No action necessary! Already found this one.
                            pass
                        elif adjacent_index < 0:
                            # This one has not been sorted yet
                            self.tag(adjacent_point)
                            self.unsorted_ones.remove(adjacent_point)
                            traversal_holding_area.add(adjacent_point)
                        else:
                            # This is in another island, so we're not including it.
                            pass
                if adjacency_count >= self.neighbor_limit:
                    self.identify(current_point)
                    self.tag(current_point)
                    traversal_stack |= traversal_holding_area
                traversal_holding_area.clear()

    def tag(self, coordinate):
        self.sorted_ones_dict[coordinate] = self._current_set_index

    def identify(self, coordinate):
        if len(self.sorted_ones_sets) == self._current_set_index:
            self.sorted_ones_sets.append(set())
        self.sorted_ones_sets[self._current_set_index].add(coordinate)

    def set_neighbor_limit(self, neighbor_limit):
        if neighbor_limit < 0 or neighbor_limit > 8:
            raise ValueError("The neighbor_limit should be between 0 and 8, not {}".format(neighbor_limit))
        else:
            self.neighbor_limit = neighbor_limit


def dilate(array, times=1, forbidden=None):
    if times == 0:
        return np.copy(array)
    if forbidden is None:
        forbidden = set()
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


def get_islands(array, n):
    b = BooleanIslands(array, neighbor_limit=n)
    return b.sorted_ones_sets

def get_mask(array, n=8, min_size=800, dilation=1):
    sorted_sets = get_islands(array, n)
    sorted_sets.sort(key=len, reverse=True)
    a_copy = np.zeros(array.shape)
    for largest_set in (s for s in sorted_sets if len(s) > min_size):
        a_copy[tuple(zip(*largest_set))] = 1
    a_copy = a_copy.astype(bool)
    return dilate(a_copy, times=dilation)


def get_planck_mask():
    fn = "bool_mask_test.fits"
    array = getdata(fn)
    array = (array == 0)
    # return array
    return get_mask(array, dilation=0)


def fill_inwards(array, nanmask, min_size=800, n=8):
    # nanmask is false if NaN
    # sorted_ones = get_islands(array&nanmask, n)
    # sorted_ones.sort(key=len, reverse=True)
    # a_copy = np.full(array.shape, False)
    # for positive_set in sorted_ones:
    #     if len(positive_set) > min_size:
    #         a_copy[tuple(zip(*positive_set))] = True
    a_copy = array.copy()
    sorted_zeros = get_islands((~a_copy)&nanmask, 0)
    sorted_zeros.sort(key=len, reverse=True)
    outside_space = sorted_zeros.pop(0)
    for negative_set in sorted_zeros:
        if len(negative_set) < min_size:
            a_copy[tuple(zip(*negative_set))] = True
    return a_copy.astype(bool)



if __name__ == "__main__":
    import matplotlib.pyplot as plt
    fn = "bool_core.fits"
    array = getdata(fn).astype(bool)
    sorted_sets = get_islands(array, 4)
    sorted_sets.sort(key=len, reverse=True)
    def gen_subset(min_size):
        a_copy = np.zeros(array.shape)
        for largest_set in (s for s in sorted_sets if len(s) > min_size):
            a_copy[tuple(zip(*largest_set))] = 1
        return a_copy.astype(bool)
    def plot_both(a_copy):
        plt.subplot(131)
        plt.imshow(array, origin='lower')
        plt.subplot(132)
        plt.imshow(a_copy, origin='lower')
        plt.subplot(133)
        plt.imshow(dilate(a_copy, times=3), origin='lower')
        plt.show()
    plot_both(gen_subset(500))
