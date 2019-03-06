#!/astromake/opt/python/anaconda3/bin/python -i
import sys
from importlib import import_module
import analyze_manticore as am


def error_exit(msg):
    print(msg)
    sys.exit(22)

configs = []
NH2lim = (19.5, 23.)
Ncutoff = []
dir_stub, src_stub = [], []
global_color, syms = [], []
pacs_offset, pacs_mods = [], []
len_count = 1
# Read in multiple config files following -c flags
while len(sys.argv) > len_count:
    if sys.argv[len_count] == "-c":
        try:
            target = sys.argv[len_count+1]
            if target[:2] == "./":
                target = target[2:]
            target = target[:-3].replace('/', '.')
            config = import_module(target)
            configs.append(config)
            dir_stub.append(config.directory)
            src_stub.append(config.stub)
            Ncutoff.append(config.NH2cutoff)
            global_color.append(config.global_color)
            syms.append(config.marker)
            pacs_offset.append(config.pacs_offset)
            pacs_mods.append(config.pacs_mods)
        except AttributeError:
            error_exit("Check your config file <%s>, something fell through" % sys.argv[len_count+1])
        except ValueError:
            print(sys.argv[len_count+1][:-3].split('/'))
    else:
        error_exit("What sys argument are you using? <%s>" % sys.argv[len_count])
    len_count += 2
if len_count == 1:
    error_exit("Need config file")

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
#plt.rc('text', usetex=True)
from astropy.io import fits
from astropy.stats import median_absolute_deviation as mad
from scipy.interpolate import UnivariateSpline

rcc = "-remapped-conv.fits"
rccl = "-remapped-conv-clean.fits"
img, err = "image", "error"
pacs160 = "PACS160um-"
spire250 = "SPIRE250um-"
spire350 = "SPIRE350um-"
spire500 = "SPIRE500um-"
bands = {
    pacs160: 'DarkViolet',
    spire250: 'b',
    spire350: 'g',
    spire500: 'r'
}
BINS = 128
diffLIM = (-10, 250)
XsMAX, XsMIN = 1.e0, 1.e-3
XsMOD = 2
XsMASK = lambda x: (x < XsMAX) & (x > XsMIN)
XsPLOTLIM = (-5, 3.1)
SKIP = 4

def plotIMG(img, ax=None):
    if ax is not None:
        plt.sca(ax)
    plt.imshow(img)
    plt.gca().invert_yaxis()
    plt.show()


def hist(img, log=True, c='k', lims=None, ax=None, alpha=1, label=None, density=False):
    imgr = img.ravel()
    if lims is not None and type(lims) is not tuple:
        imgr = imgr[imgr < lims]
    if ax is not None:
        plt.sca(ax)
    plt.hist(imgr[~np.isnan(imgr)], bins=BINS, log=log,
             histtype='step', fill=False, color=c, range=lims,
             alpha=alpha, label=label, normed=density)


def fn_constr(b, index, im=True):
    """
    This retrieves the original flux & error fits files
    """
    directory = dir_stub[index]
    ie = img if im else err
    end_stub = rccl if src_stub[index][-1] == 'F' else rcc
    return directory + b + ie + end_stub


def result_constr(plus, div, index):
    """
    This retrieves the manticore results
    """
    if plus == -99:
        return configs[index].filename
    stub = dir_stub[index] + "T4-absdiff-" + src_stub[index]
    if plus == -1:
        stub += "-4bandLErr"
    elif plus == "G":
        stub += "-plusGRD"
    elif plus != 0 or div != 1:
        stub += "-plus%03d" % plus
        if div != 1:
            stub += "div%02d" % div
    return stub + ".fits"


def pixel_histograms(index, axes=None):
    with fits.open(result_constr(-1, 1, index)) as hdul:
        N = np.log10(hdul[3].data)
    mask = N < Ncutoff[index]
    if axes is None:
        plt.figure()
        regAx, maskAx = plt.subplot(211), plt.subplot(212)    
    else:
        regAx, maskAx = axes
    l = 1
    legend = []
    for b in bands:
        c = bands[b]
        legend.append(c)
        with fits.open(fn_constr(b, index)) as hdul:
            flux = hdul[0].data
            #if b == bands[0]:
            #    flux += 35
            hist(flux, log=l, c=c, ax=regAx)
            hist(flux[mask], log=l, c=c, ax=maskAx)
    plt.legend(legend)
    regAx.set_title("unfiltered, below 100MJy/sr")
    maskAx.set_title("low density mask, below 100MJy/sr")


def modify_HDU(modifier, index, im=True):
    b = pacs160
    data, header = fits.getdata(fn_constr(b, index, im=im), header=True)
    if im:
        data += modifier
        comment = "added %d to these 160 um fluxes" % modifier
        new_stub = "plus%03d-" % modifier
        ie = img
    else:
        data /= modifier
        comment = "divided this 160 um error map by factor of %d" % modifier
        new_stub = "div%02d-" % modifier
        ie = err
    header['HISTORY'] = comment
    end_stub = rccl if src_stub[index][-1] == 'F' else rcc
    new_fn = dir_stub[index] + b + new_stub + ie + end_stub
    fits.writeto(new_fn, data, header, overwrite=True)
    print(new_fn)


def modify_HDU_slope(p1, p2, val1, val2, index, constant_offset=0, custom=False):
    """
    Apply a gradient offset to PACS160um image.
    Supply two points and two diff160 values;
     the function will extrapolate a linear gradient parallel to
     the line passing between those two points and will assume
     value val1 at point p1 and val2 at p2.
    p1 and p2 should be in the form (x, y), and val1&2 should be
     numerical (int, float, etc).
    The x value of p1 MUST be less than the x value of p1!!!!
     so adjust your gradient accordingly (maybe reverse).
    """
    m = (p2[1]-p1[1])/(p2[0]-p1[0])  # Slope of parallel-to-gradient line
    b = (p1[1] - p1[0]*m)  # Zero-offset of parallel-to-gradient line
    b1 = (p1[1] + p1[0]/m)  # Zero-offset of perpendicular-to-gradient line
    # Distance from perpendicular-gradient line passing through p1=(x1, y1)
    dist = lambda i, j: (-1 if m<0 else 1)*((i/m)+j-b1)/np.sqrt((1/m**2)+1)
    b = pacs160
    data, header = fits.getdata(fn_constr(b, index), header=True)
    comment = "added "
    # Apply constant, if applicable (if we read gradient off already-offset result)
    if constant_offset != 0:
        data += constant_offset
        comment += "constant %d and "
    # Apply gradient
    xx, yy = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]), sparse=False)
    grad_offset = (dist(xx, yy)/dist(p2[0], p2[1]))*(val2-val1) + val1
    grad_offset[np.isnan(data)] = np.nan
    if custom:  # This is throwaway, don't keep this
        p3, p4 = (587, 1192), (964, 464)
        m2 = (p4[1]-p3[1])/(p4[0]-p3[0])
        b2 = (p3[1] - p3[0]*m2)
        grad_offset[(yy < (m2*xx + b2)) & (~np.isnan(data))] = 0
#    return grad_offset, data
    data += grad_offset
    # Save
    args = (p1[0], p1[1], val1, p2[0], p2[1], val2)
    comment += "slope: (%d, %d):%.2f -> (%d, %d):%.2f" % args
    header['HISTORY'] = comment
    new_stub = "plusGR2-"
    end_stub = rccl if src_stub[index][-1] == 'F' else rcc
    new_fn = dir_stub[index] + b + new_stub + img + end_stub
    fits.writeto(new_fn, data, header)
    print(new_fn)
    plt.imshow(grad_offset, cmap='hot', origin='lower')
    plt.show()


def chisq_by_N(index, ax=None, errorbars=False):
    """
    Chi Squared for logN bins.
    Column density N is taken from the large 160um error run.
    Chi Squared is taken from various runs where 160um maps
    have been modified in such a way to see the effect on chi squared.

    This currently only works for Perseus; it could be expanded out to
    others if more modified 160um runs are performed.
    """
    if ax is None:
        plt.figure()
        axAvg = plt.subplot(111)
    else:
        axAvg = ax
    all_modifiers, all_formatters = [], []
    plus_mods, colors, div_mods, symbols = pacs_mods[index]
    print(plus_mods, div_mods)
    try:
        with fits.open(result_constr(-1, 1, index)) as hdul:
            N = np.log10(hdul[3].data)
    except FileNotFoundError:
        with fits.open(result_constr(-99, 1, index)) as hdul:
            N = np.log10(hdul[3].data)
    Nhist, Nedges = np.histogram(N.ravel(), bins=BINS, range=NH2lim)
    prep_arr = lambda a, b: np.array([a, b]).T.flatten()
    histx, histy = prep_arr(Nedges[:-1], Nedges[1:]), prep_arr(Nhist, Nhist)
    bin_centers = np.empty(Nhist.shape)
    for i, lims in enumerate(zip(Nedges[:-1], Nedges[1:])):
        bin_centers[i] = np.mean(lims)
    for pm in plus_mods:
        for dm in div_mods:
            all_modifiers.append((pm, dm))
    for pm in colors:
        for dm in symbols:
            all_formatters.append((pm, dm))
    prep_chsq = np.log10 # lambda x: x
    for modifier, fmt in zip(all_modifiers, all_formatters):
        c, s = fmt
        try:
            with fits.open(result_constr(modifier[0], modifier[1], index)) as hdul:
                chisq = prep_chsq(hdul[5].data)
            print("found: ", modifier)
        except FileNotFoundError:
            if modifier[0]==0 and modifier[1]==1:
                print("Trying default instead of", modifier)
                try:
                    with fits.open(result_constr(-99, 1, index)) as hdul:
                        chisq = prep_chsq(hdul[5].data)
                except FileNotFoundError:
                    continue
            else:
                print("did not find:", modifier)
                continue
        ch_avgs, ch_stds = np.empty(Nhist.shape), np.empty(Nhist.shape)
        ch_meds, ch_mads = np.empty(Nhist.shape), np.empty(Nhist.shape)
        for i, lims in enumerate(zip(Nedges[:-1], Nedges[1:])):
            ch_inBin = chisq[(N >= lims[0]) & (N < lims[1])]
            if ch_inBin.size == 0:
                ch_avgs[i] = ch_stds[i] = ch_meds[i] = np.NaN
            else:
                ch_avgs[i], ch_stds[i] = np.mean(ch_inBin), np.std(ch_inBin)
                ch_meds[i], ch_mads[i] = np.median(ch_inBin), mad(ch_inBin)
        # Averages:
        #axAvg.errorbar(bin_centers, ch_avgs, fmt=s, c=c, yerr=ch_stds, capsize=2)
        # Medians:
        label = "F+%02d, $\sigma$/%02d " % (modifier[0], modifier[1])
        if errorbars:
            axAvg.errorbar(bin_centers, ch_meds, fmt=s, label=label,
                           mfc=c, mec='k', mew=1, ms=5, ecolor=c,
                           yerr=ch_mads, capsize=2, alpha=0.75)
        else:
            axAvg.plot(bin_centers, ch_meds, s, label=label,
                       mfc=c, mec='k', mew=1, ms=5, alpha=0.75)
    axAvg.legend()
    axAvg.set_title("Median $\chi^{2}$ per log(N) for %s" % src_stub[index])
    axAvg.set_ylabel("Median $\chi^{2}$")
    axAvg.set_xlabel("[Column density N]")
    axHist = axAvg.twinx()
    axHist.plot(histx, histy, '-', color='k', alpha=0.5)
    axHist.set_ylabel("N bin pixel count")

def hist_160diff(ax=None, alpha=1):
    """
    Histograms of the full diff160 maps
    """
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    for i, cloud in enumerate(src_stub):
        c = global_color[i]
        stub = result_constr(-1, 0, i)
        with fits.open(stub) as hdul:
            diff160 = hdul[6].data
        hist(diff160, log=False, lims=diffLIM, c=c, ax=ax, alpha=alpha, label=cloud)
    ax.legend()
    ax.set_xlabel("Diff 160 Flux")
    ax.set_ylabel("Pixel count")
    ax.set_title("Full-image diff160 histogram")

def chisq_mask(index, ax=None, plot=True):
    """
    Chi Squared and diff160 are obtained from the fit
    during which 160um was given a massive error.
    The diff160 thus shows what the SPIRE fit wants PACS to have.
    Low Chi Squared regions should indicate places where PACS
    should fit well to the single temp fit by SPIRE.
    """
    if ax is None and plot:
        plt.figure()
        ax = plt.subplot(111)
    with fits.open(result_constr(-1, 0, index)) as hdul:
        chisq = hdul[5].data * XsMOD
        diff160 = hdul[6].data
    mask = XsMASK(chisq) & (~np.isnan(diff160))
    dhist, dedges = np.histogram(diff160[mask].ravel(), bins=BINS, range=diffLIM)
    prep_arr = lambda a, b: np.array([a, b]).T.flatten()
    histx, histy = prep_arr(dedges[:-1], dedges[1:]), prep_arr(dhist, dhist)
    bin_centers = (dedges[:-1]+dedges[1:])/2
    peak_val = np.max(dhist)
    med = bin_centers[dhist == peak_val]  # was: dedges[dhist == peak_val]; why?
    try:
        spline = UnivariateSpline(bin_centers, dhist - peak_val/2, s=0)
        r1, r2 = spline.roots()
        fwhm = np.abs(r1 - r2)
        sigma = fwhm/2.355
    except ValueError:
        fwhm, sigma = np.nan, np.nan
    if plot:
        ax.plot(histx, histy, '-', color=global_color[index])
        ax.text(0.1, 0.9, "Mode: %.1f\nSTD: %.2f"
                % (med, sigma),
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes)
        ax.set_xlabel("Diff 160 Flux")
        ax.set_ylabel("Pixel count")
        ax.set_title("%s diff160 with $%.1E < \chi^{2} < %.1E$ mask" % (src_stub[index], XsMIN, XsMAX),
                     fontsize=10)
    else:
        print("%s: Mode: %.1f, STD: %.2f" % (src_stub[index], med, sigma))

def noiseFraction(index, ax=None):
    evf_c = 25000
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    legend = []
    for b in bands:
        c = bands[b]
        i_stub, e_stub = fn_constr(b, index, im=True), fn_constr(b, index, im=False)
        flux_map = fits.getdata(i_stub)
        if b == pacs160:
            flux_map += pacs_offset[index]
        flux_map = np.abs(flux_map)
        error_map = np.abs(fits.getdata(e_stub))
        evf = error_map/flux_map
        evf = evf.flatten()
        flux_map = flux_map.flatten()
        mask = ~(np.isnan(evf) | np.isnan(flux_map))
        evf = evf[mask]
        flux_map = flux_map[mask]
        mask = (evf < evf_c) & (evf > -evf_c)
        evf = evf[mask]
        flux_map = flux_map[mask]
        ax.plot(np.log10(flux_map)[::20], np.log10(evf)[::20], ',', c=c, alpha=0.005)
        legend.append(mpatches.Patch(color=c, label=b[:-1]))
    ax.legend(handles=legend)
    ax.set_xlabel("Flux (MJy/sr)")
    ax.set_ylabel("Error (% of flux)")
    ax.set_title("Noise Fraction for %s (with %d MJy/sr added to 160um)"
                 % (src_stub[index], pacs_offset[index]))


def noiseHist(index, b=pacs160, ax=None):
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    error_map = fits.getdata(fn_constr(b, index, im=False))
    med = np.nanmedian(error_map)
    label = "%s: Median $\sigma$ = %.2f" % (src_stub[index], med)
    hist(error_map, log=False, c=global_color[index], ax=ax, lims=(0, 35), label=label)
    ax.set_xlabel("160um error bins")
    ax.set_ylabel("Pixel count")
    ax.set_title("Stated 160um error map histogram")

def chisqHist(index, ax=None, fourband=False):
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    with fits.open(result_constr(-1+fourband, 1, index)) as hdul:
        chisq = hdul[5].data * (1 if fourband else XsMOD)
    chisq = np.log10(chisq)
    label = src_stub[index]
    hist(chisq, log=False, c=global_color[index], ax=ax, label=label,
         lims=XsPLOTLIM, density=True)
    ax.set_xlabel("[$\chi^{2}$]")
    ax.set_ylabel("Pixel count (normalized)")
    ax.set_title("$\chi^{2}$ Histogram")

def scaHist_XsMask(index, fign=(1, 2), fourband=False):
    fig_SH, fig_I = fign
    am.set_globals(configs[index])
    axSctr, axHistN, axHistT = am.setup_scaHist(fig_SH)
    am.setup_image(fig_I)
    if not fourband:
        with fits.open(result_constr(-1, 1, index)) as hdul:
            T = hdul[1].data
            N = hdul[3].data
            chisq = hdul[5].data * XsMOD
    else:
        with fits.open(result_constr(-1, 1, index)) as hdul:
            T = hdul[1].data
            N = hdul[3].data
        with fits.open(result_constr(-1+fourband, 1, index)) as hdul:
            chisq = hdul[5].data
    am.scatter(T, N, ax=axSctr)
    am.hist(T, N, axes=(axHistT, axHistN), log=False)
    rgb = am.imaging(T, N, fig=fig_I)
    mask = XsMASK(chisq) #& (np.log10(N) < Ncutoff[index])
    rgb[~mask] = np.array([0, 0, 1])
    plt.figure(fig_I)
    plt.imshow(rgb)
    plt.gca().invert_yaxis()
    Tbox, Nbox = T[mask], N[mask]
    clr = 'r'
    am.scatter(Tbox, Nbox, ax=axSctr, c=clr)
    medT, medN = am.hist(Tbox, Nbox, axes=(axHistT, axHistN), c=clr, log=False)
    am.format_scatter(ax=axSctr)
    am.format_hist(scax=axSctr, axes=(axHistT, axHistN))
    print("%s Median T: %.2f K" % (src_stub[index], medT))


def scaHist_diffXs(index, legend, axes=None, fourband=False):
    Tlim, NH2lim = am.Tlim, am.NH2lim
    am.Tlim, am.NH2lim = diffLIM, XsPLOTLIM  # hacky but easy. not actually a T or N lim.
    axSctr, axHistXs, axHistd = axes
    with fits.open(result_constr(-1+fourband, 1, index)) as hdul:
        chisq = hdul[5].data * (1 if fourband else XsMOD)
        diff160 = hdul[6].data
    clr = global_color[index]
    am.scatter(diff160, chisq, ax=axSctr, c=clr, skip=SKIP)
    am.hist(diff160, chisq, axes=(axHistd, axHistXs), log=False, c=clr)
    legend.append(mpatches.Patch(color=clr, label=src_stub[index]))
    plt.sca(axSctr)
    plt.xlabel("[$\chi^{2}$]")
    plt.ylabel("diff160")
    plt.ylim(diffLIM)
    plt.xlim(XsPLOTLIM)
    am.format_hist(scax=axSctr, axes=(axHistd, axHistXs))
    am.Tlim, am.NH2lim = Tlim, NH2lim  # clean up the mess


def scaXsTN(index, legend=([], []), ax=None, fourband=False):
    if ax is None:
        plt.figure()
        ax = (plt.subplot(121), plt.subplot(122))
    with fits.open(result_constr(-1+fourband, 1, index)) as hdul:
        T = hdul[1].data
        N = np.log10(hdul[3].data)
        chisq = hdul[5].data * (1 if fourband else XsMOD)
    def gen_fixarr(m):
        return lambda x: x[m].ravel()[::SKIP]
    alpha = 0.1
    l1, l2 = legend
    legend = l1
    ax1, ax2 = ax
    ax = ax1
    maxXs = 200
    mask = ~(np.isnan(T) | np.isnan(chisq)) & (chisq < 10**XsPLOTLIM[1])
    fixarr = gen_fixarr(mask)
    ax.plot(np.log10(fixarr(chisq)), fixarr(T), ',', color=global_color[index], alpha=0.1)
    legend.append(mpatches.Patch(color=global_color[index], label=src_stub[index]))
    ax.set_xlabel("$[\chi^{2}]$")
    ax.set_ylabel("T (K)")
    ax.set_title("T vs $\chi^{2}$")
#    ax.set_xlim([0, maxXs])
    ax.set_xlim(XsPLOTLIM)
    ax.set_ylim(list(am.Tlim))
    legend = l2
    ax = ax2
    mask = ~(np.isnan(N) | np.isnan(chisq)) & (chisq < 10**XsPLOTLIM[1])
    fixarr = gen_fixarr(mask)
    ax.plot(np.log10(fixarr(chisq)), fixarr(N), ',', color=global_color[index], alpha=0.1)
    legend.append(mpatches.Patch(color=global_color[index], label=src_stub[index]))
    ax.set_xlabel("$[\chi^{2}]$")
    ax.set_ylabel("[N]")
    ax.set_title("log(N) vs $\chi^{2}$")
#    ax.set_xlim([0, maxXs])
    ax.set_xlim(XsPLOTLIM)
    ax.set_ylim(list(am.NH2lim))
    

def temp_XsMask(index, fig=None, Tlim=(15, 18), fourband=False):
    plt.figure(fig)
    with fits.open(result_constr(-1+fourband, 1, index)) as hdul:
        T = hdul[1].data
        chisq = hdul[5].data * (1 if fourband else XsMOD)
    mask = np.ma.masked_where(XsMASK(chisq), np.ones(chisq.shape))
    Tlo, Thi = Tlim
    T[T < Tlo] = Tlo
    T[T > Thi] = Thi
    plt.imshow(T, cmap='hot', origin='lower')
    plt.colorbar()
    plt.imshow(mask, cmap='jet', origin='lower')


def diffXsMask(index, fig=None, difflim=(0, 150), fourband=False):
    plt.figure(fig)
    with fits.open(result_constr(-1+fourband, 1, index)) as hdul:
        T = hdul[1].data
        N = np.log10(hdul[3].data)
        chisq = hdul[5].data * (1 if fourband else XsMOD)
        diff160 = hdul[6].data
    dlo, dhi = difflim
    diffMASK = lambda x: ~np.isnan(x)#(x < dhi) & (x > dlo)
    mask = np.ma.masked_where(XsMASK(chisq) & diffMASK(diff160), np.ones(chisq.shape))
    diff160[diff160 < dlo] = dlo
    diff160[diff160 > dhi] = dhi
    plt.imshow(diff160, cmap='hot', origin='lower')
    plt.colorbar()
    plt.imshow(mask, cmap='jet', origin='lower')


def scaHist_diffFlux(index, legend=[], axes=None, dlim=diffLIM, flim=(-100, 100), f_option=0):
    # f_option = 0: real flux vs diff
    # f_option = 1: real flux + calc offset vs diff
    # f_option = 2: predicted flux vs diff
    # f_option = 3: predicted flux + calc offset vs diff
    log = False
    axSctr, axHistf, axHistd = axes
    div = (f_option%2==1)+1 if src_stub[:-1] == 'F' else 1
    with fits.open(result_constr(pacs_mods[index][0][f_option%2==1], div, index)) as hdul:
        chisq = hdul[5].data
        diff160 = hdul[6].data
        band160 = hdul[10].data
    clr = global_color[index]
    notnanmask = ~(np.isnan(diff160) | np.isnan(band160))# & XsMASK(chisq)
    prep_arr = lambda x: x[notnanmask].ravel()[::SKIP]
    band160, diff160 = prep_arr(band160 + diff160) if f_option>1 else prep_arr(band160), prep_arr(diff160)
    axSctr.plot(band160, diff160, ',', color=clr, alpha=0.2)
    for ax, data, r, o in zip((axHistf, axHistd), (band160, diff160), (flim, dlim), ('vertical', 'horizontal')):
        ax.hist(data, bins=am.BINS, log=log, histtype='step',
                fill=False, color=clr, range=r, orientation=o)
    legend.append(mpatches.Patch(color=clr, label=src_stub[index]))


def scaHist_Xstbfb(index, legend=[], axes=None):
    log = False
    axSctr, axHistx, axHisty = axes
    with fits.open(result_constr(-1, 1, index)) as hdul:
        chisq3 = hdul[5].data
    with fits.open(result_constr(pacs_offset[index], 1, index)) as hdul:
        chisq4 = hdul[5].data
        diff160 = hdul[6].data # for NaN mask where no PACS
    clr = global_color[index]
    notnanmask = ~(np.isnan(chisq3) | np.isnan(chisq4) | np.isnan(diff160))
    prep_arr = lambda x: np.log10(x[notnanmask].ravel()[::SKIP])
    chisq3, chisq4 = prep_arr(chisq3), prep_arr(chisq4)
    axSctr.plot(chisq3, chisq4, ',', color=clr, alpha=0.2)
    for ax, data, o in zip((axHistx, axHisty), (chisq3, chisq4), ('vertical', 'horizontal')):
        ax.hist(data, bins=am.BINS, log=log, histtype='step',
                fill=False, color=clr, range=XsPLOTLIM, orientation=o)
    legend.append(mpatches.Patch(color=clr, label=src_stub[index]))


def temp_XsMasktbfb(index, fig=None, Tlim=(15, 18), chisq4lim=XsPLOTLIM):
    plt.figure(fig)
    with fits.open(result_constr(-1, 1, index)) as hdul:
        T = hdul[1].data
        chisq3 = np.log10(hdul[5].data * XsMOD)
    with fits.open(result_constr(pacs_offset[index], 1, index)) as hdul:
        chisq4 = np.log10(hdul[5].data)
        diff160 = hdul[6].data # for NaN mask where no PACS
    xslo, xshi = chisq4lim
    mask = (chisq4 < xshi) & (chisq4 > xslo)
    mask = np.ma.masked_where(XsMASK(chisq3) & ~np.isnan(diff160) & mask, np.ones(chisq3.shape))
    Tlo, Thi = Tlim
    T[T < Tlo] = Tlo
    T[T > Thi] = Thi
    plt.imshow(T, cmap='hot', origin='lower')
    plt.colorbar()
    plt.imshow(mask, cmap='jet', origin='lower')



###################################################
######### Quick-Run Macro Functions ###############
###################################################

def run_NF():
    """
    Noise fraction scatter plot
    Plots all regions (up to 6) and noise of each of 4 bands in log space
    """
    plt.figure()
    for i in range(min(len(src_stub), 6)):
        noiseFraction(i, plt.subplot(3, 2, i+1))
    plt.show()

def run_XsN():
    """
    Median chi squared as function of logN
    Across several modified 160um maps
    Plots all regions (up to 6)
    """
    plt.figure()
    for i in range(min(len(src_stub), 6)):
        ax = plt.subplot(3, 2, i+1)
        chisq_by_N(i, ax=ax, errorbars=1)
    plt.show()

def run_XsMask():
    """
    The dispersion of 160um 'preferred' diffs masked with XsMASK
    Also includes median error in stated error map
    Plots all regions (up to 5)
    Uses chi squared from SPIRE-only exclusively
    """
    plt.figure()
    ax1 = plt.subplot2grid((3, 6), (0, 0), colspan=2)
    ax2 = plt.subplot2grid((3, 6), (0, 2), colspan=2)
    ax3 = plt.subplot2grid((3, 6), (0, 4), colspan=2)
    ax4 = plt.subplot2grid((3, 6), (1, 1), colspan=2)
    ax5 = plt.subplot2grid((3, 6), (1, 3), colspan=2)
    ax6 = plt.subplot2grid((3, 6), (2, 0), colspan=3)
    ax7 = plt.subplot2grid((3, 6), (2, 3), colspan=3)
    axes = [ax1, ax2, ax3, ax4, ax5]
    for i in range(min(len(src_stub), 5)):
        chisq_mask(i, axes[i])
        noiseHist(i, ax=ax7)
    if len(src_stub) > 5:
        chisq_mask(5, axes[src_stub.index('Per2')])
        chisq_mask(6, axes[src_stub.index('Per2')])
        axes[src_stub.index('Per2')].legend(['Per2', src_stub[5], src_stub[6]])
    if len(src_stub) > 7:
        chisq_mask(7, axes[src_stub.index('Ser2')])
        axes[src_stub.index('Ser2')].legend(['Ser2', src_stub[7]])
    hist_160diff(ax=ax6)
    ax7.legend(loc='upper left')
    plt.show()

def run_TNscat(n, fb=False):
    """
    The scatter-histogram of T vs N that we were looking at during TEST3
    Includes a second overlay of T vs N given the XsMASK
    Also plots up the density image masked with XsMASK in second figure
    Uses chi squared, T, & N from SPIRE-only by default, but supports
     fb=True for fourband chi squared. T & N will remain from SPIRE-only!
    Use SKIP global variable to reduce scatter points.
    """
    scaHist_XsMask(n, fourband=fb)
    plt.show()

def run_dXscat(n, fb=False):
    """
    Scatter-histogram of diff160 vs chi squared.
    SPIRE-only by default, but supports fb=True for four band diff160 & chi squared
    """
    axSctr, axHistXs, axHistd = am.setup_scaHist(None)
    legend = []
    if type(n) == int:
        n = [n]
    for i in n:
        scaHist_diffXs(i, legend, axes=(axSctr, axHistXs, axHistd), fourband=fb)
    axSctr.legend(handles=legend)
    plt.show()

def run_dfscat(n, offset=0):
    """
    Scatter-histogram of diff160 versus 160 band flux.
    Supports several flux options in the offset=# keyword:
    - offset = 0: real flux
    - offset = 1: real flux + calc offset
    - offset = 2: predicted flux
    - offset = 3: predicted flux + calc offset
    Reduce scatter points with SKIP global variable.
    """
    axSctr, axHistf, axHistd = am.setup_scaHist(None)
    legend = []
    diffLim = tuple([x - 100 for x in diffLIM]) if offset%2==1 else diffLIM
    fluxLim = (-100, 300) 
    if type(n) == int:
        n = [n]
    for i in n:
        scaHist_diffFlux(i, legend, axes=(axSctr, axHistf, axHistd), dlim=diffLim, flim=fluxLim, f_option=offset)
    axSctr.set_xlabel("160um flux")
    axSctr.set_ylabel("diff160")
    axSctr.set_ylim(diffLim)
    axSctr.set_xlim(fluxLim)
    axSctr.legend(handles=legend)
    am.format_hist(scax=axSctr, axes=(axHistd, axHistf))
    plt.show()

def run_Xs34Bscat(n):
    """
    Scatter-histogram; compares the chi squared from the
     SPIRE-only fit to the chi squared in the same pixels
     from the four band fit, including the PACS offset.
    Can be used to deduce what will happen to the chi squared
     after offset in the four band region within the original
     SPIRE-motivated chi squared mask.
    """
    axSctr, axHistx, axHisty = am.setup_scaHist(None)
    legend = []
    if type(n) == int:
        n = [n]
    for i in n:
        scaHist_Xstbfb(i, legend, axes=(axSctr, axHistx, axHisty))
    axSctr.set_xlabel("SPIRE-only $[\chi^{2}]$")
    axSctr.set_ylabel("4-band(PACS+offset) $[\chi^{2}]$")
    axSctr.set_xlim(XsPLOTLIM)
    axSctr.set_ylim(XsPLOTLIM)
    am.format_hist(scax=axSctr, axes=(axHisty, axHistx))
    plt.show()

def run_XsHist(fb=False):
    """
    Histogram (normalized) of log chi squared for all regions present.
    SPIRE-only by default, but fb=True sets it to four band chi squared.
    """
    plt.figure()
    ax = plt.subplot(111)
    for i in range(len(src_stub)):
        chisqHist(i, ax=ax, fourband=fb)
    ax.legend(loc='upper left')
    plt.show()

def run_Tmap(n, Tlim=(15, 18), fb=False):
    """
    Shows T image masked with XsMASK and scaled with Tlim.
    SPIRE-only by default, but fb=True will switch to four band T & chi squared.
    """
    temp_XsMask(n, Tlim=Tlim, fourband=fb)
    plt.show()

def run_Tmap34B(n, Tlim=(15, 18), X4lim=XsPLOTLIM):
    """
    Shows T image scaled with Tlim and masked using BOTH SPIRE-only chi squared
     as well as four band PACS+offset chi squared.
    Control SPIRE-only mask with XsMIN, XsMAX (LINEAR)
    Control four band mask with the function keyword X4lim=(lo, hi) (LOG10)
    """
    temp_XsMasktbfb(n, Tlim=Tlim, chisq4lim=X4lim)
    plt.show()


def run_diffXsmap(n, difflim=diffLIM, fb=False):
    """
    Shows diff160 image masked with XsMASK and scaled with difflim
    SPIRE-only by default, but fb=True will switch to four band diff160 & chi squared.
    """
    diffXsMask(n, difflim=difflim, fourband=fb)
    plt.show()

def run_XsTN(n, fb=False):
    """
    Chi squared versus T in one frame, chi squared versus N in a second frame.
    SPIRE-only by default, but fb=True switches T, N, & chi squared all over to four band.
    Reduce scatter points with SKIP global variable.
    """
    if type(n) == int:
        n = [n]
    plt.figure()
    ax = (plt.subplot(121), plt.subplot(122))
    legend = ([], [])
    for i in n:
        scaXsTN(i, legend=legend, ax=ax, fourband=fb)
    [a.legend(handles=l) for a, l in zip(ax, legend)]
    plt.show()

def run_XsMask_fast():
    """
    Really fast printout of the mode diff160s under XsMASK.
    These can theoretically be used as the offset.
    Same numbers as printed in the corner of run_XsMask().
    """
    for i in range(len(src_stub)):
        chisq_mask(i, plot=False)

modify_HDU(46, 0, im=True)
#run_NF()
#run_XsN()
#run_XsMask()
print("""Plot functions:
-- run_NF()   |   run_XsN()   |   run_XsMask()   |   run_XsMask_fast(),

-- run_TNscat(n, fb=1/0)   |   run_XsHist(fb=1/0)   |   run_Tmap(n, Tlim=(lo, hi), fb=1/0),

-- run_XsTN(n, fb=1/0)   |   run_diffXsmap(n, difflim=(lo, hi), fb=1/0)   |   run_dXscat(n, fb=1/0),

-- run_dfscat(n, offset=0-3)   |   run_Xs34Bscat(n)

-- run_Tmap34B(n, Tlim=(lo, hi), X4lim=(lo, hi))

-- Modify PACS functions:
-- modify_HDU_(##, n, im=True)
-- modify_HDU_slope((x1,y1), (x2,y2), val1, val2, n, constant_offset=##):

-- help(function) for help""")
print("n can be:")
for e in enumerate(src_stub):
    print("  %d: %s" % e)
#plt.show()
