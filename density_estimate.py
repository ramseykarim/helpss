plotting_remotely = False
import numpy as np
import matplotlib
if plotting_remotely:
	matplotlib.use('Agg')
SAVE_NAME = "Figure_X_current.png"
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

def show_plot():
	if plotting_remotely:
		plt.savefig(SAVE_NAME)
	else:
		plt.show()

# # Department directories
# working_dir = "/n/sgraraid/filaments/data/TEST4/helpss_scratch_work/"
# soln = "T4-absdiff-Per1J-3param-plus046-full.fits"
# # Local directories
working_dir = "/home/ramsey/Documents/Research/Filaments/"
soln = "T4-absdiff-Per1J-plus045-pow-1000-0.1-1.80.fits"

# Load image
img = fits.getdata(working_dir + soln, 3)
img[np.isnan(img)] = 0
# Rescale
img[img > 0] = np.sqrt(img[img > 0])
# img[img > 0] -= np.min(img[img > 0])
# Flatten
img_flat = img.flatten(order='C')
cumulative_img = np.cumsum(img_flat)
cdf = cumulative_img / cumulative_img[-1]
print("N0", np.min(img_flat[img_flat > 0]))
print("N1", np.max(img_flat[img_flat > 0]))
print("SNR N1/N0", np.max(img_flat[img_flat > 0])/np.min(img_flat[img_flat > 0]))
print("sum(N)", cumulative_img[-1])
print("ratio N0/sum(N)", np.min(img_flat[img_flat > 0])/cumulative_img[-1])
print("ratio N1/sum(N)", np.max(img_flat[img_flat > 0])/cumulative_img[-1])
print("ratio 5e22/sum(N)", np.sqrt(5e22)/cumulative_img[-1])
print("n", cdf.size)
best_n_points = int(cumulative_img[-1]*3*3/(9*np.sqrt(1e20)))
print("best n", best_n_points)
# pdf = img_flat / cumulative_img[-1]

def get_position(i):
	# Assume the center of the pixel
	row = i // img.shape[1]
	col = i % img.shape[1]
	return np.concatenate([col+0.5, row+0.5], axis=1)

def sample_kernel(size=1, sigma=1.0):
	# Return two samples from a 2D gaussian kernel
	return np.random.normal(scale=sigma, size=(size, 2))

n_points = 10000
n_points = best_n_points
kernel_width = 3 # pixels
rand_cdf_values = np.random.uniform(size=(n_points, 1))
pixel_indices = np.searchsorted(cdf, rand_cdf_values)
points = get_position(pixel_indices) + sample_kernel(size=n_points, sigma=kernel_width)

# plt.figure(figsize=(16, 9))
# plt.plot(points[:, 0], points[:, 1], ',', markersize=3, alpha=0.05)
# plt.xlim([0, img.shape[1]])
# plt.ylim([0, img.shape[0]])
# plt.gca().set_aspect('equal')
# show_plot()

working_dir = "/home/ramsey/Documents/Research/AlphaX/PyAlpha_drafting/test_data/"
np.savetxt("{:s}filament{:04d}_sampleNH2_betah1.80.dat".format(working_dir, n_points),
	 points)
print("done")
