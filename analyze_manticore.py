#!/astromake/opt/python/anaconda3/bin/python
import sys
from importlib import import_module


def error_exit(msg, flag, whatever):
    """
    Helper function for quick error reporting while reading command line args
    MSG is intended to be a string, somewhat of a description
    FLAG is the flag "-f", "-c", etc, where the hangup occured
    WHATEVER is intended to be the offending object, maybe sys.argv[1:]
    """
    print(msg)
    print("Fix your "+flag+" flag usage. Problem with: ", whatever)
    sys.exit()


def set_globals(imported_config):
    global fits_fn
    fits_fn = imported_config.filename
    global src_stub
    src_stub = imported_config.stub + "_"
    global limits
    limits = imported_config.limits
    global color
    color = imported_config.color
    try:
        global dir_stub
        dir_stub = imported_config.directory
    except AttributeError:
        print("Using current directory ./")
    try:
        global gradientSearch
        gradientSearch = imported_config.gradientSearch
    except AttributeError:
        print("no gradient search, ", end="")
    try:
        global scatMaskHelp
        scatMaskHelp = imported_config.scatMaskHelp
    except AttributeError:                
        print("no scatter mask, ", end="")
    try:
        global NH2cutoff
        NH2cutoff = imported_config.NH2cutoff
    except AttributeError:
        print("default NH2cutoff, ", end="")
    print("NH2 cutoff = %.2f" % NH2cutoff)
    

"""
These defaults will not allow the program to proceed
"""
dir_stub = "./"
fits_fn = ""
src_stub = ""
limits = {}
color = {}
median_T = 15.5
Tlim = (8, 45)
NH2lim = (19.5, 23.)
NH2cutoff = 21.25
gradientSearch = None
scatMaskHelp = None
if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == '-f':
            try:
                fits_fn = sys.argv[2]
            except IndexError:
                error_exit("-f", sys.argv[1:])
        elif sys.argv[1] == "-c":
            try:
                target = sys.argv[2]
                if target[:2] == "./":
                    target = target[2:]
                target = target[:-3].replace('/', '.')
                config = import_module(target)
                set_globals(config)
            except IndexError:
                error_exit("", "-c", sys.argv[1:])
            except ModuleNotFoundError:
                error_exit("Config file not found.", "-c", target)
            except AttributeError as e:
                error_exit("Check config file variables: "+sys.argv[2], "-c", str(e))
    else:
        msg = "Need fits file or config file.\n"
        msg += "Use -f flag to specify fits FILENAME. This usage does not support boxed region analysis.\n"
        msg += "\tThe fits file must be the results of a MANTICORE run.\n"
        msg += "Use -c flag to specify config FILENAME. Config files require 4 variable definitions:\n"
        msg += "\t1) filename: string name of MANTICORE result fits file.\n"
        msg += "\t2) limits: string->list(list(int)) dictionary of boxed region limits. \"region_name\": [[x0, x1], [y0, y1]].\n"
        msg += "\t3) color: string->string dictionary of associated region colors across all plots. Keys must MATCH limits keys.\n"
        msg += "\t4) stub: string name for insertion into PNG filenames. Should be concise & appropriate for filenames."
        error_exit(msg, "-f or -c", sys.argv)


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.ticker import NullFormatter
nullfmt = NullFormatter()
BINS = 128

"""
############################################################# REGULAR PLOTS
"""
def hist(T, NH2, fname="plot_hists.png", axes=None, c='k', log=True):
    """
    Plot histograms for T and NH2 images
    If AXES is specified, fname is disregarded
    AXES, if specified, is a tuple: (T_ax, NH2_ax)
    """
    if axes is None:
        plt.figure()
        plt.subplot(211)
    else:
        plt.sca(axes[0])
    prep_arr = lambda x: x[~np.isnan(x)].ravel()
    nT, bins, patches = plt.hist(prep_arr(T), bins=BINS, log=log,
                                 histtype='step', fill=False, color=c, range=Tlim,
                                 orientation='horizontal')
    bin_centers = (bins[:-1]+bins[1:])/2
    medT = bin_centers[nT == np.max(nT)]
    if axes is None:
        plt.subplot(212)
    else:
        plt.sca(axes[1])
    nN, bins, patches = plt.hist(prep_arr(np.log10(NH2)), bins=BINS, log=log,
                                 histtype='step', fill=False, color=c, range=NH2lim,
                                 orientation='vertical')
    bin_centers = (bins[:-1]+bins[1:])/2
    medN = bin_centers[nN == np.max(nN)]
    if axes is None:
        plt.savefig(fname)
    return medT, medN  # returns medians of histograms


def format_hist(scax=None, axes=None):
    if axes is None:
        return
    else:
        plt.sca(axes[0])
        plt.ylim(scax.get_ylim())
        plt.sca(axes[1])
        plt.xlim(scax.get_xlim())


def setup_avgT_Nbins(fig_T_Nbins, figsize=(9, 10.5)):
    plt.figure(fig_T_Nbins, figsize=figsize)
    histAx, avgTAx, stdTAx = plt.subplot(311), plt.subplot(312), plt.subplot(313)
    return histAx, avgTAx, stdTAx


def avgT_Nbins(T, NH2, fname=None, axes=None, c='k'):
    if axes is None:
        fig = plt.figure()
        histAx, avgTAx, stdTAx = setup_avgT_Nbins(fig.number)
    else:
        histAx, avgTAx, stdTAx = axes
    notnanmask = ~(np.isnan(T) | np.isnan(NH2))
    T, NH2 = T[notnanmask], NH2[notnanmask]
    logN = np.log10(NH2)
    Nhist, Nedges = np.histogram(logN.ravel(), bins=BINS, range=NH2lim)
    prep_arr = lambda a, b: np.array([a, b]).T.flatten()
    histx, histy = prep_arr(Nedges[:-1], Nedges[1:]), prep_arr(Nhist, Nhist)
    histAx.plot(histx, histy, '-', c=c)
    bin_centers = np.empty(Nhist.shape)
    Tavgs, Tstds = np.empty(Nhist.shape), np.empty(Nhist.shape)
    for i, lims in enumerate(zip(Nedges[:-1], Nedges[1:])):
        bin_centers[i] = np.mean(lims)
        T_inBin = T[(logN >= lims[0]) & (logN < lims[1])]
        if T_inBin.size == 0:
            Tavgs[i], Tstds[i] = np.NaN, np.NaN
        else:
            Tavgs[i], Tstds[i] = np.mean(T_inBin), np.std(T_inBin)
    avgTAx.plot(bin_centers, Tavgs, '.', color=c)
    stdTAx.plot(bin_centers, Tstds, '.', color=c)
    if axes is None:
        if fname is None:
            fname = name_constr(src_stub, "TavgNbins", dir_stub=dir_stub)
        format_avgT_Nbins((histAx, avgTAx, stdTAx))
        plt.savefig(fname)


def format_avgT_Nbins(axes=None):
    if axes is None:
        return
    histAx, avgTAx, stdTAx = axes
    histAx.set_yscale('log')
    histAx.set_xlim(list(NH2lim))
    histAx.set_ylabel("N(H2) bin counts")
    avgTAx.set_xlim(list(NH2lim))
    avgTAx.set_ylabel("AVG(T)")
    stdTAx.set_xlim(list(NH2lim))
    stdTAx.set_ylabel("STD(T)")
    stdTAx.set_xlabel("log(N(H2)) bin")
    histAx.xaxis.set_major_formatter(nullfmt)
    avgTAx.xaxis.set_major_formatter(nullfmt)


def imbox(lims, fig=None, c='y'):
    """
    Plot boxes inside image given x,y limits
    LIMS must be [[x0, x1], [y0, y1]]
    And keep in mind that Python does x, y -> row, col
      so x might be vertical.. trial and error!
    FIG must be specified or else the program will exit.
    """
    if fig is None:
        return
    else:
        plt.figure(fig)
    ylim, xlim = lims
    for i in range(2):
        plt.plot([xlim[i], xlim[i]], [ylim[0], ylim[1]], '-', color=c)
        plt.plot([xlim[0], xlim[1]], [ylim[i], ylim[i]], '-', color=c)


def imaging(T, NH2, fname="plot_images.png", fig=None, noplot=False):
    """
    Plot RGB images given T and NH2
    
    RED: temperature above 16 K; normalized (16, 45)
    GREEN: NH2 column density; normalized 10^(21, 23)
    BLUE: temperature below 16 K; normalized (16, 10)

    If fig is specified, fname is disregarded
    """
    lo_T, hi_T = Tlim
    Nlo, Nhi = NH2lim

    Thi, Tlo = hi_T, median_T
    Thot = T.copy()
    Thot[Thot < Tlo] = Tlo
    Thot[Thot > Thi] = Thi
    Thot -= Tlo
    Thot /= Thi-Tlo
    
    Thi, Tlo = median_T, lo_T
    Tcold = T.copy()
    Tcold[Tcold < Tlo] = Tlo
    Tcold[Tcold > Thi] = Thi
    Tcold -= Tlo
    Tcold /= Thi-Tlo
    
    Tcold = 1 - Tcold
    
    NH2log = np.log10(NH2)
    NH2log[NH2log < Nlo] = Nlo
    NH2log[NH2log > Nhi] = Nhi
    NH2log = (NH2log - Nlo)/(Nhi - Nlo)
    
    g = np.zeros(Thot.shape)
    rgb = np.stack([Thot, NH2log, Tcold], axis=2)
    if noplot:
        return rgb
    if fig is None:
        plt.figure()
    else:
        plt.figure(fig)
    plt.imshow(rgb)
    plt.gca().invert_yaxis()
    if fig is None:
        plt.savefig(fname)
    return rgb
    

def scatter(T, NH2, fname="plot_scatter.png", ax=None, c='k', alpha=0.05, label=None, skip=1):
    """
    Plot T vs log10(NH2) scatter for given field

    If fig is specified, fname is disregarded
    """
    if ax is None:
        plt.figure()
    else:
        plt.sca(ax)
    notnanmask = ~(np.isnan(T) | np.isnan(NH2))
    plt.plot(np.log10(NH2[notnanmask].ravel()[::skip]), T[notnanmask].ravel()[::skip], ',', color=c, alpha=alpha, label=label)
    if ax is None:
        plt.savefig(fname)


def setup_scaHist(fig_scaHist, figsize=(15, 15)):
    # Plot axes setup; from matplotlib example scatter_hist.html
    anchor, width = 0.1, 0.62
    anchor_h = anchor + width + 0.03
    rect_scatter = [anchor, anchor, width, width]
    rect_histNH2 = [anchor, anchor_h, width, 0.22]
    rect_histT = [anchor_h, anchor, 0.2, width]
    plt.figure(fig_scaHist, figsize=figsize)
    axSctr = plt.axes(rect_scatter)
    axHistNH2, axHistT = plt.axes(rect_histNH2), plt.axes(rect_histT)
    axHistNH2.xaxis.set_major_formatter(nullfmt)
    axHistT.yaxis.set_major_formatter(nullfmt)
    return axSctr, axHistNH2, axHistT


def setup_image(fig_img, figsize=(10, 10)):
    plt.figure(fig_img, figsize=figsize)

def format_scatter(ax=None):
    if ax is None:
        return
    else:
        plt.sca(ax)
        plt.xlabel("N(H$_{2}$ (cm$^{-2}$)")
        plt.ylabel("T (K)")
        plt.xlim(list(NH2lim))
        plt.ylim(list(Tlim))

"""
############################################################# MASKS
"""
# This runs the low density mask and creates the hist/scatter plots
def redblue_mask(T, NH2, fname=None, fig=None, noplot=False):
    mask = np.log10(NH2) < NH2cutoff
    if noplot:
        return mask
    Tmask = T.copy()
    Tlo, Thi = 15, 18
    Tmask[Tmask < Tlo] = Tlo
    Tmask[Tmask > Thi] = Thi
    Tmask -= Tlo
    Tmask /= Thi-Tlo
    dens = np.log10(NH2)
    Nlo, Nhi = 20.9, 21.5
#    dens[dens < NH2lim[0]] = NH2lim[0]
#    dens[dens > NH2lim[1]] = NH2lim[1]
#    dens -= NH2lim[0]
#    dens /= NH2lim[1]-NH2lim[0]
    dens = np.zeros(T.shape)
    dens[~mask] = 1
    rgb = np.stack([Tmask, np.zeros(T.shape), dens], axis=2)
    if fig is None:
        plt.figure()
    else:
        plt.figure(fig)
    plt.imshow(rgb)
    plt.gca().invert_yaxis()
    if fig is None:
        if fname is None:
            fname = name_constr(src_stub, "mask", dir_stub=dir_stub)
        plt.savefig(fname)

# Functionality for checking out regions of the scatter plot in the image
def scatterMask(T, NH2, fname=None, fig=None):
    if fig is None:
        f = plt.figure()
    rgb = imaging(T, NH2, noplot=True)
    masks, clrs = scatMaskHelp(T, np.log10(NH2))
    clrs = tuple(map(lambda x: np.array(x)/255., clrs))
    for msk, new_rgb in zip(masks, clrs):
        rgb[msk] = new_rgb
    plt.imshow(rgb)
    plt.gca().invert_yaxis()
    if fname is None:
        fname = name_constr(src_stub, "search", dir_stub=dir_stub)
    plt.savefig(fname)


"""
############################################################# CREATE PLOTS
"""
np.warnings.filterwarnings('ignore')
if __name__ == "__main__":
    fig_scaHist = 1
    fig_img = 2
    fig_T_Nbins = 3
    axSctr, axHistNH2, axHistT = setup_scaHist(fig_scaHist)
    setup_image(fig_img)
    histAx, avgTAx, stdTAx = setup_avgT_Nbins(fig_T_Nbins)
    
    T, NH2 = None, None
    
    with fits.open(fits_fn) as hdul:
        print(hdul.info())
        T = hdul[1].data
        NH2 = hdul[3].data

    plots = {
        fig_scaHist: "scHst",
        fig_img: "image",
        fig_T_Nbins: "Tstat",
        }

    plot_stub = "plot_"
    png_stub = ".png"


    def name_constr(field_stub, plot_type, dir_stub='./'):
        return dir_stub + plot_stub + field_stub + plot_type + png_stub


    # All-pixel plots
    scatter(T, NH2, ax=axSctr)
    hist(T, NH2, axes=(axHistT, axHistNH2))
    rgb = imaging(T, NH2, fig=fig_img)
    avgT_Nbins(T, NH2, axes=(histAx, avgTAx, stdTAx))

    # Select larger areas using masks set with gradientSearch
    if gradientSearch is not None:
        masks, clrs, lines = gradientSearch(T.shape)
        overlay_mask = ~np.isnan(T)#redblue_mask(T, NH2, noplot=True)
        plt.figure(fig_img)
        rgb[~overlay_mask] = np.array([0, 0, 1])
        plt.imshow(rgb)
        plt.gca().invert_yaxis()
        for msk, clr in zip(masks, clrs):
            msk = msk & overlay_mask
            Tbox = T[msk]
            NH2box = NH2[msk]
            scatter(Tbox, NH2box, ax=axSctr, c=clr)
            hist(Tbox, NH2box, axes=(axHistT, axHistNH2), c=clr)
            avgT_Nbins(Tbox, NH2box, axes=(histAx, avgTAx, stdTAx), c=clr)
        plt.figure(fig_img)
        for l in lines:
            plt.plot(l[1], l[0], '--', c='y')

    # Boxed regions set with limits dict
            
    for field in limits:
        xlim, ylim = limits[field]
        Tbox = T[xlim[0]:xlim[1], ylim[0]:ylim[1]]
        NH2box = NH2[xlim[0]:xlim[1], ylim[0]:ylim[1]]
        scatter(Tbox, NH2box, ax=axSctr, c=color[field])
        hist(Tbox, NH2box, axes=(axHistT, axHistNH2), c=color[field])
        imbox(limits[field], fig=fig_img, c=color[field])
        avgT_Nbins(Tbox, NH2box, axes=(histAx, avgTAx, stdTAx), c=color[field])
        
    # This part writes out the multicolored plots
    format_scatter(ax=axSctr)
    format_hist(scax=axSctr, axes=(axHistT, axHistNH2))
    format_avgT_Nbins(axes=(histAx, avgTAx, stdTAx))
    for fig in (fig_scaHist, fig_img, fig_T_Nbins):
        plt.figure(fig)
        plt.savefig(name_constr(src_stub, plots[fig], dir_stub=dir_stub))

    # N-masked, rescaled temperature image
    redblue_mask(T, NH2)
