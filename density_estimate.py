plotting_remotely = False
import numpy as np
import matplotlib
if plotting_remotely:
	matplotlib.use('Agg')
SAVE_NAME = "Figure_X_current.png"
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from astropy.io import fits
from planck_mask import gen_hist_and_stats
import manticore_results as mtc
import sys

def show_plot():
	if plotting_remotely:
		plt.savefig(SAVE_NAME)
	else:
		plt.show()

# working_dir = "/n/sgraraid/filaments/data/TEST4/helpss_scratch_work/"
working_dir = "/home/ramsey/Documents/Research/Filaments/"
# soln = "T4-absdiff-Per1J-3param-plus046-full.fits"
soln = "T4-absdiff-Per1J-plus045-pow-1000-0.1-1.80.fits"

T = fits.getdata(working_dir + soln, 3)
# mask = mtc.get_filament_mask()
# T[~mask] = 0
# T = np.log(T)
# T -= np.nanmin(T)
# T += 1
# T = np.exp(T)
T[np.isnan(T)] = 0

Tflat = T.flatten(order='C')
Tsort_i = np.argsort(Tflat)
Tsort = Tflat[Tsort_i]

# plt.figure(figsize=(11, 8.5))
# plt.subplot(121)

# new_T = Tflat.reshape(T.shape)
cumulative_Tsort = np.cumsum(Tsort)
CDF_Tsort = cumulative_Tsort / cumulative_Tsort[-1]
PDF_Tsort = Tsort / cumulative_Tsort[-1]
## Plot the CDF
# plt.plot(Tsort, CDF_Tsort, '.')
# show_plot()
## Proof of rebuilding image from flattened array
# empty_map = np.empty(T.shape)
# for i, val in zip(Tsort_i, PDF_Tsort):
#	 empty_map[i // T.shape[1], i % T.shape[1]] = val
# plt.imshow(empty_map, origin='lower')
# plt.colorbar()


def get_position(i):
	# Assume the center of the pixel
	row = i // T.shape[1]
	col = i % T.shape[1]
	return np.array([col+0.5, row+0.5])

def sample_kernel(sigma=1.0):
	# Return two samples from a 2D gaussian kernel
	return np.random.normal(scale=sigma, size=2)

# plt.imshow(T, origin='lower')
# plt.subplot(122)
n_points = 7000
kernel_width = 3 # pixels
points = []
for i in range(n_points):
	rand_cdf_value = np.random.uniform()
	chosen_cdf_value = np.min(CDF_Tsort[CDF_Tsort > rand_cdf_value])
	rand_idx = np.where(chosen_cdf_value == CDF_Tsort)
	rand_i = Tsort_i[rand_idx][0]
	# x, y ordering is column, row
	xy = get_position(rand_i)
	xy += sample_kernel(sigma=kernel_width)
	points.append(xy)
points = np.array(points)
# plt.plot(points[:, 0], points[:, 1], ',', markersize=3)
# plt.xlim([0, T.shape[1]])
# plt.ylim([0, T.shape[0]])
# plt.gca().set_aspect('equal')
# show_plot()
working_dir = "/home/ramsey/Documents/Research/AlphaX/PyAlpha_drafting/test_data/"
np.savetxt("{:s}filament{:04d}_sampleNH2_betah1.80.dat".format(working_dir, n_points),
	 points)
print("done")
