import numpy as np
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io.fits import getdata
import manticore_results as mtc
import boolean_islands as bli
from planck_mask import get_spire_mask
import sys

SAVE_NAME = "Figure_X_current.png"
spire_mask = get_spire_mask(setting=2)


cold_fn = "bool_cold.fits"
core_fn = "bool_core.fits"
planck_fn = "bool_mask_test.fits"
cold_mask = getdata(cold_fn)
cold_mask = cold_mask.astype(bool)
# cold_mask = bli.get_mask(cold_mask, n=4, dilation=0)

# core_mask = getdata(core_fn)
# core_mask = core_mask.astype(bool)
# core_mask = bli.get_mask(core_mask, n=4, dilation=0)

planck_mask = getdata(planck_fn)
planck_mask = (planck_mask == 0)
# planck_mask = bli.get_mask(planck_mask, n=4, dilation=0, min_size=1500)

# r_proxy = planck_mask
# g_proxy = spire_mask
# b_proxy = cold_mask

# all_F = ~(r_proxy | g_proxy | b_proxy)
# all_T = r_proxy & g_proxy & b_proxy
# r = ((r_proxy | all_F) & ~all_T).astype(float)
# g = ((g_proxy | all_F) & ~all_T).astype(float)
# b = ((b_proxy | all_F) & ~all_T).astype(float)

# r = (r_proxy).astype(float)
# g = (g_proxy).astype(float)
# b = (b_proxy).astype(float)

# both_mask = np.stack([r, g, b], axis=2)

all_mask = (spire_mask.astype(int) + cold_mask.astype(int) + planck_mask.astype(int))
# all_mask = (all_mask >= 3).astype(bool)
min_size = spire_mask.shape[0]*spire_mask.shape[1]//1000
print("1pct of the image is {}".format(min_size))
# all_mask = bli.get_mask(all_mask, n=5, min_size=min_size, dilation=3)

plt.figure(figsize=(14, 8))
plt.subplot(111)
plt.imshow(all_mask, origin='lower')
plt.tight_layout()
plt.show()
