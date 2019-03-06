import sys
from importlib import import_module

dir_stub, src_stub = [], []
configs = []
len_count = 1
while len(sys.argv) > len_count:
    if sys.argv[len_count] == '-c':
        try:
            target = sys.argv[len_count+1]
            if target[:2] == "./":
                target = target[2:]
            target = target[:-3].replace('/', '.')
            config = import_module(target)
            configs.append(config)
            dir_stub.append(config.directory)
            src_stub.append(config.stub)
        except AttributeError:
            print("check config file")
            sys.exit()
        except ValueError:
            print("value error?")
    else:
        print("what sys args are you using?")
    len_count += 2
if len_count == 1:
    print("need config file")
    sys.exit()


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, PowerNorm
from astropy.io import fits

NH2lim = (19.5, 25.)
BINS = 2**7

gen_fn_stub = lambda s: "T4-absdiff-"+s+"-3param.fits"
def gen_fn(index):
    directory = dir_stub[index]
    fn_stub = gen_fn_stub(src_stub[index])
    return directory + fn_stub

def hist(img, log=True, c='k', lims=None, ax=None, alpha=1, label=None, density=False):
    imgr = img.ravel()
    imgr = imgr[~np.isnan(imgr)]
    if lims is not None and type(lims) is not tuple:
        imgr = imgr[imgr < lims]
        lims = None
    plt.hist(imgr, bins=BINS, log=log, histtype='step',
             fill=False, color=c, range=lims, alpha=alpha, label=label,
             normed=density)

def hist2d(imgs, fs=None, log=True, cmap=None, lims=None, ax=None, normed=True):
    c = 'k'
    alpha = 1
    imgr_list = []
    nanmask = np.isnan(imgs[0]) | np.isnan(imgs[1])
    for img in imgs:
        img = img[~nanmask]
        imgr = img.ravel()
        imgr_list.append(imgr)
    if fs is not None:
        imgr_list = list(fs(*tuple(imgr_list)))
    imgrx, imgry = imgr_list
#    plt.figure()
#    plt.hist(imgrx, bins=BINS, log=True, histtype='step',
#             fill=False, color=c, range=lims_list[0],
#             normed=False)
#    plt.figure()
#    plt.hist(imgry, bins=BINS, log=True, histtype='step',
#             fill=False, color=c, range=lims_list[1],
#             normed=False)
    plt.hist2d(imgrx, imgry, bins=BINS, range=lims, norm=LogNorm())

Nlimhi = 25

def modify_N(x):
    sign_arr = np.sign(x)
    abs_arr = np.abs(x)
    abs_arr = np.log10(abs_arr)
    return sign_arr * abs_arr


def cold_temp_hist(index):
    with fits.open(gen_fn(index)) as hdul:
        Tc = hdul[1].data
        #Nc = hdul[3].data
        #Nh = hdul[7].data
        #chisq = hdul[10].data
    hist(Tc, lims=20)
    plt.show()

def cold_gas_hist2d(index):
    with fits.open(gen_fn(index)) as hdul:
        Tc = hdul[1].data
        Nc = hdul[3].data
        #Nh = hdul[7].data
        #chisq = hdul[10].data
    f = lambda x, y: (modify_N(x), y)
    hist2d((Nc, Tc), fs=f, lims=([-Nlimhi, Nlimhi], [-1, 20]))
    plt.xlabel("log10(Nc)")
    plt.ylabel("Tc")
    plt.show()

def cold_gas_hist2d_REARRANGE(index):
    with fits.open(gen_fn(index)) as hdul:
        Tc = hdul[1].data
        Nc = hdul[3].data
        #Nh = hdul[7].data
        #chisq = hdul[10].data
    def f(x, y):
        x_mod = modify_N(x)
        y_mod = y*np.sign(x_mod)
        x_mod = np.abs(x_mod)
        return x_mod, y_mod
    hist2d((Nc, Tc), fs=f, lims=(list(NH2lim), [-20, 20]))
    plt.xlabel("log10(Nc), kinda")
    plt.ylabel("Tc, kinda")
    plt.show()

def warm_v_cold_hist2d(index):
    with fits.open(gen_fn(index)) as hdul:
        Nc = hdul[3].data
        Nh = hdul[7].data
    f = lambda x, y: (modify_N(x), modify_N(y))

    plt.figure()
    for i in range(2):
        for j in range(2):
            lims = [[-1*x if not i else x for x in NH2lim], [-1*x if j else x for x in NH2lim]]
            lims = list(map(sorted, lims))
            # false j: y axis is positive
            # true i: x axis is positive
            # i: toggles right side
            # j: toggles bottom row
            plt.subplot(2, 2, j*2 + i + 1)
            if not i:
                if not j:
                    plt.ylabel("positive Nh")
                if j:
                    plt.ylabel("negative Nh")
            if j:
                if not i:
                    plt.xlabel("negative Nc")
                if i:
                    plt.xlabel("positive Nc")
#            for l in lims:
#                print(l)
            if not (j==1 and i==0):
                hist2d((Nc, Nh), fs=f, lims=lims)
            else:
                plt.xlim(lims[0]), plt.ylim(lims[1])
            print()
    plt.show()


def plot_with_mask(index, plot=None, masks=None, lims=None, f=None):
    instructions = """
plot=PLOT_NUMBER (1-4)
masks=(4 element tuple of functions of 1 numeric variable)
lims=optional tuple (lo, hi) of color scale limits
f=optional function to apply to the plotted data (e.g. log10 for column density, etc)

Mask functions should all evaluate to True if you want to see that pixel
If function is None it is assumed to be all True
Even if all masks are None, still need a 4-element tuple
1: Tc
2: Nc
3: Nh
4: chisq
"""
    if plot is None or masks is None:
        print(instructions)
        return
    with fits.open(gen_fn(index)) as hdul:
        Tc = hdul[1].data
        Nc = hdul[3].data
        Nh = hdul[7].data
        chisq = hdul[10].data
    frames = [Tc, Nc, Nh, chisq]
    names = ['Tc', 'Nc', 'Nh', 'chi sq']
    arr_masks = []
    for frame, mask in zip(frames, masks):
        if mask is not None:
            arr_masks.append(mask(frame))
    master_mask = np.all(np.array(arr_masks), axis=0) if len(arr_masks) > 1 else None
    plot -= 1
    img = frames[plot]
    if f is not None:
        img = f(img)
    # TODO imshow with mask (in fix160.py)
    if lims is not None:
        lo, hi = lims
        img[img < lo] = lo
        img[img > hi] = hi
    plt.imshow(img, cmap='hot', origin='lower')
    plt.colorbar()
    if master_mask is not None:
        mask = np.ma.masked_where(master_mask, np.ones(master_mask.shape))
        plt.imshow(mask, cmap='jet', origin='lower')
    plt.title(names[plot])
    plt.show()

def debug(index):
    with fits.open(gen_fn(index)) as hdul:
        Tc = hdul[1].data
        Nc = hdul[3].data
        #Nh = hdul[7].data
        #chisq = hdul[10].data
    def f(x, y):
        x_mod = modify_N(x)
        y_mod = y*np.sign(x_mod)
        x_mod = np.abs(x_mod)
        return x_mod, y_mod
    Nc_m, Tc_m = f(Nc, Tc)
    plt.imshow(Nc_m)
    plt.show()
    

print("""
cold_temp_hist(index)
cold_gas_hist2d(index)
cold_gas_hist2d_REARRANGE(index)
warm_v_cold_hist2d(index)
plot_with_mask(index, plot=None, masks=None, lims=None)
""")
