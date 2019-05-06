import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import sys
import boolean_islands as boolis
from datetime import datetime, timezone
from matplotlib.font_manager import FontProperties

per1_dir = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
power_stub = "pow-1000-0.1-"
power_run_stub = "T4-absdiff-Per1J-plus045-"
power_run_3p_stub = "T4-absdiff-Per1J-3param-plus045-"
single_comp_dir = "single_comp_beta_grid/"
two_comp_dir_stub = "two_comp_beta_grid_betah"
fits_stub = ".fits"

GNILC_T_fn = per1_dir + "dustModel_Per1_SPIRE500umgrid_Temperature.fits"
multi_beta_fn = per1_dir + "single_component_multi_beta.fits"

def gen_power_fn(beta):
	return "{:s}{:s}{:s}{:s}{:s}{:s}".format(
		per1_dir, single_comp_dir, power_run_stub,
		power_stub, beta, fits_stub
	)

def gen_3p_power_fn(beta_c, beta_h):
	return "{:s}{:s}{:s}/{:s}c{:s}{:s}h{:s}{:s}{:s}".format(
		per1_dir, two_comp_dir_stub, beta_h,
		power_run_3p_stub,
		power_stub, beta_c,
		power_stub, beta_h,
		fits_stub
	)

filament_mask_fn = per1_dir+"filament_mask_syp.fits"

local_test = per1_dir+"T4-absdiff-Per1J-3param-plus046-full.fits"
b21n1e20 = per1_dir+"T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-bcreinit-Nh1E20.fits"
b21b16 = per1_dir+"T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.10hpow-1000-0.1-1.60-bcreinit-Th16.0-Nh5E19,1E22.fits"
b21n0 = per1_dir+"T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-2.10hpow-1000-0.1-1.80-bcreinit-Nh0.0.fits"
b15n1e20 = per1_dir+"T4-absdiff-Per1J-3param-plus045-cpow-1000-0.1-1.50hpow-1000-0.1-1.80-bcreinit-Nh1E20.fits"

mask = fits.getdata(filament_mask_fn).astype(bool)

target_3p_fn = b21n1e20
target_2p_fn = gen_power_fn("1.80")

font_mono = FontProperties()
font_mono.set_family('monospace')
font_mono.set_variant('normal')
font_serif_large = FontProperties()
font_serif_large.set_family('serif')
font_serif_large.set_variant('normal')
font_serif_large.set_size('x-large')
font_serif_small = font_serif_large.copy()
font_serif_small.set_size('medium')


def make_2p_plot():
    frames_to_plot = (1, 3)
    frame2p_data = []
    with fits.open(target_2p_fn) as hdul:
        hdr = hdul[1].header
        for frame in frames_to_plot:
            frame2p_data.append(hdul[frame].data)
    frames_to_plot = frame2p_data
    frames_to_plot[1] = np.log10(frames_to_plot[1])
    limits = [
        (10.5, 17.5), # T
        (20.3, 22.7), # N
    ]
    titles = [
        r"Single-component T",
        r"Single-component N(H$_{2}$)",
    ]
    c_labels = [
        r"Temperature (K)",
        r"log$_{10}$[Column density ($H_{2}$/cm$^2$)]",
    ]
    ticks = [
        (11, 12, 13, 14, 15, 16, 17),
        (20.5, 21, 21.5, 22, 22.5),
    ]
    w = WCS(hdr)
    size = np.array((10, 17.5))*0.6
    plt.figure(figsize=size)
    for i in range(2):
        ax = plt.subplot(211+i, projection=w)
        lo, hi = limits[i]
        plt.imshow(frames_to_plot[i], origin='lower',
            vmin=lo, vmax=hi,
            cmap='inferno')
        if i > 0:
            plt.xlabel("Right Ascension", fontproperties=font_serif_small)
        plt.ylabel("Declination", fontproperties=font_serif_small)
        plt.title(titles[i], fontproperties=font_serif_large)
        c = plt.colorbar(pad=0.0, ticks=list(ticks[i]))
        c.ax.set_yticklabels([str(x) for x in ticks[i]])
        for t in c.ax.get_yticklabels():
            t.set_fontproperties(font_serif_small)
            t.set_rotation(-90.)
        c.set_label(c_labels[i],
            fontproperties=font_serif_small, rotation=-90.,
            labelpad=18)
    plt.subplots_adjust(
        top=0.96,
        bottom=0.051,
        left=0.08,
        right=0.97,
        hspace=0.18,
        wspace=0.125
    )
    plt.savefig("/home/rkarim/Desktop/syp_figure_2p.png")

def make_larger2p_plot():
    frames_to_plot = (1, 2, 3, 4, 5)
    frame2p_data = []
    with fits.open(target_2p_fn) as hdul:
        hdr = hdul[1].header
        for frame in frames_to_plot:
            frame2p_data.append(hdul[frame].data)
    frames_to_plot = frame2p_data
    frames_to_plot[2] = np.log10(frames_to_plot[2])
    frames_to_plot[3] = np.log10(frames_to_plot[3])
    limits = [
        (10.5, 17.5), # T
        (0, 1), # dT
        (20.3, 22.7), # N
        (19.5, 21.5), # dN
        (0, 5), # Xs
    ]
    titles = [
        r"Single-component T",
        r"Single-component T uncertainty",
        r"Single-component N(H$_{2}$)",
        r"Single-component N(H$_{2}$) uncertainty",
        r"Single-component reduced $\chi^{2}$",
    ]
    c_labels = [
        r"Temperature (K)",
        r"Temperature (K)",
        r"log$_{10}$[Column density ($H_{2}$/cm$^2$)]",
        r"log$_{10}$[Column density ($H_{2}$/cm$^2$)]",
        r"$\chi^{2}$/DoF",
    ]
    ticks = [
        (11, 12, 13, 14, 15, 16, 17),
        (0, 0.25, 0.5, .75, 1),
        (20.5, 21, 21.5, 22, 22.5),
        (19.5, 20, 20.5, 21, 21.5),
        (0, 0.25, 0.5, 0.75),
    ]
    w = WCS(hdr)
    plt.figure(figsize=(10, 15))
    for i in range(5):
        ax = plt.subplot(321+i, projection=w)
        lo, hi = limits[i]
        plt.imshow(frames_to_plot[i], origin='lower',
            vmin=lo, vmax=hi,
            cmap='inferno')
        plt.xlabel("Right Ascension", fontproperties=font_serif_small)
        plt.ylabel("Declination", fontproperties=font_serif_small)
        plt.title(titles[i], fontproperties=font_serif_large)
        c = plt.colorbar(pad=0.0, ticks=list(ticks[i]))
        c.ax.set_yticklabels([str(x) for x in ticks[i]])
        for t in c.ax.get_yticklabels():
            t.set_fontproperties(font_serif_small)
            t.set_rotation(-90.)
        c.set_label(c_labels[i],
            fontproperties=font_serif_small, rotation=-90.,
            labelpad=18)
    # plt.subplots_adjust(
    #     top=0.905,
    #     bottom=0.085,
    #     left=0.08,
    #     right=0.935,
    #     hspace=0.2,
    #     wspace=0.125
    # )
    plt.show()

def make_T2p_plot():
    frames_to_plot = (1, 2)
    frame2p_data = []
    with fits.open(target_2p_fn) as hdul:
        hdr = hdul[1].header
        for frame in frames_to_plot:
            frame2p_data.append(hdul[frame].data)
    frames_to_plot = frame2p_data
    limits = [
        (10.5, 17.5), # T
        (0, 1), # dT
    ]
    titles = [
        r"Single-component T",
        r"Single-component T uncertainty",
    ]
    c_labels = [
        r"Temperature (K)",
        r"Temperature uncertainty (K)",
    ]
    ticks = [
        (11, 12, 13, 14, 15, 16, 17),
        (0, 0.25, 0.5, .75, 1),
    ]
    w = WCS(hdr)
    plt.figure(figsize=(9.7, 17.5))
    for i in range(2):
        ax = plt.subplot(211+i, projection=w)
        lo, hi = limits[i]
        plt.imshow(frames_to_plot[i], origin='lower',
            vmin=lo, vmax=hi,
            cmap='inferno')
        plt.xlabel("Right Ascension", fontproperties=font_serif_small)
        plt.ylabel("Declination", fontproperties=font_serif_small)
        plt.title(titles[i], fontproperties=font_serif_large)
        c = plt.colorbar(pad=0.0, ticks=list(ticks[i]))
        c.ax.set_yticklabels([str(x) for x in ticks[i]])
        for t in c.ax.get_yticklabels():
            t.set_fontproperties(font_serif_small)
            t.set_rotation(-90.)
        c.set_label(c_labels[i],
            fontproperties=font_serif_small, rotation=-90.,
            labelpad=18)
    plt.subplots_adjust(
        top=0.96,
        bottom=0.04,
        left=0.08,
        right=0.97,
        hspace=0.1,
        wspace=0.125
    )
    plt.savefig("/home/rkarim/Desktop/syp_figure_T2p.png")


def make_N2p_plot():
    frames_to_plot = (3, 4)
    frame2p_data = []
    with fits.open(target_2p_fn) as hdul:
        hdr = hdul[1].header
        for frame in frames_to_plot:
            frame2p_data.append(hdul[frame].data)
    frames_to_plot = frame2p_data
    frames_to_plot[0] = np.log10(frames_to_plot[0])
    frames_to_plot[1] = np.log10(frames_to_plot[1])
    limits = [
        (20.3, 22.7), # N
        (19.5, 21.2), # dN
    ]
    titles = [
        r"Single-component N(H$_{2}$)",
        r"Uncertainty in the single-component N(H$_{2}$)",
    ]
    c_labels = [
        r"log$_{10}$[Column density ($H_{2}$/cm$^2$)]",
        r"log$_{10}$[Column density ($H_{2}$/cm$^2$)]",
    ]
    ticks = [
        (20.5, 21, 21.5, 22, 22.5),
        (19.5, 20, 20.5, 21),
    ]
    w = WCS(hdr)
    plt.figure(figsize=(9.7, 17.5))
    for i in range(2):
        ax = plt.subplot(211+i, projection=w)
        lo, hi = limits[i]
        plt.imshow(frames_to_plot[i], origin='lower',
            vmin=lo, vmax=hi,
            cmap='inferno')
        plt.xlabel("Right Ascension", fontproperties=font_serif_small)
        plt.ylabel("Declination", fontproperties=font_serif_small)
        plt.title(titles[i], fontproperties=font_serif_large)
        c = plt.colorbar(pad=0.0, ticks=list(ticks[i]))
        c.ax.set_yticklabels([str(x) for x in ticks[i]])
        for t in c.ax.get_yticklabels():
            t.set_fontproperties(font_serif_small)
            t.set_rotation(-90.)
        c.set_label(c_labels[i],
            fontproperties=font_serif_small, rotation=-90.,
            labelpad=18)
    plt.subplots_adjust(
        top=0.96,
        bottom=0.04,
        left=0.08,
        right=0.97,
        hspace=0.1,
        wspace=0.125
    )
    plt.savefig("/home/rkarim/Desktop/syp_figure_N2p.png")

def make_Xs2p_plot():
    frames_to_plot = (5,)
    frame2p_data = []
    with fits.open(target_2p_fn) as hdul:
        hdr = hdul[1].header
        for frame in frames_to_plot:
            frame2p_data.append(hdul[frame].data)
    frames_to_plot = frame2p_data
    frames_to_plot[0] = np.log10(frames_to_plot[0])
    limits = [
        (np.log10(0.5), np.log10(10)), # Xs
    ]
    titles = [
        r"Single-component reduced $\chi^{2}$",
    ]
    c_labels = [
        r"$\chi^{2}$/DoF",
    ]
    ticks = [
        (np.log10(0.5), np.log10(1), np.log10(2), np.log10(5), np.log10(10)),
    ]
    tick_labels = [
        (0.5, 1, 2, 5, 10)
    ]
    w = WCS(hdr)
    plt.figure(figsize=(10, 10))
    for i in range(1):
        ax = plt.subplot(111+i, projection=w)
        lo, hi = limits[i]
        plt.imshow(frames_to_plot[i], origin='lower',
            vmin=lo, vmax=hi,
            cmap='inferno')
        plt.xlabel("Right Ascension", fontproperties=font_serif_small)
        plt.ylabel("Declination", fontproperties=font_serif_small)
        plt.title(titles[i], fontproperties=font_serif_large)
        c = plt.colorbar(pad=0.0, ticks=list(ticks[i]))
        c.ax.set_yticklabels([str(x) for x in tick_labels[i]])
        for t in c.ax.get_yticklabels():
            t.set_fontproperties(font_serif_small)
            t.set_rotation(-90.)
        c.set_label(c_labels[i],
            fontproperties=font_serif_small, rotation=-90.,
            labelpad=18)
    plt.subplots_adjust(
        top=0.91,
        bottom=0.08,
        left=0.08,
        right=1.02,
        hspace=0.0,
        wspace=0.0
    )
    plt.savefig("/home/rkarim/Desktop/syp_figure_Xs2p.png")


def make_3p_plot():
    frames_to_plot = (1, 3, 5, 7)
    frame3p_data = []
    with fits.open(target_3p_fn) as hdul:
        hdr = hdul[1].header
        for frame in frames_to_plot:
            frame3p_data.append(hdul[frame].data)

    frames_to_plot = (1, 3)
    frame2p_data = []
    with fits.open(target_2p_fn) as hdul:
        for frame in frames_to_plot:
            frame2p_data.append(hdul[frame].data)
    Th = frame2p_data[0].copy()
    Th[mask] = frame3p_data[2][mask]
    Nh = frame2p_data[1].copy()
    Nh[mask] = frame3p_data[3][mask]
    Nh = np.log10(Nh)
    Tc = frame3p_data[0].copy()
    Tc[~mask] = np.nan
    Nc = frame3p_data[1].copy()
    Nc[~mask] = np.nan
    Nc = np.log10(Nc)

    frames_to_plot = [Th, Nh, Tc, Nc]

    limits = [
        (13.5, 17.5), # Th
        (20.3, 22.2), # Nh
        (8.5, 12.5), # Tc
        (20.8, 22.7), # Nc
    ]

    ticks = [
        (14, 15, 16, 17),
        (20.5, 21, 21.5, 22),
        (9, 10, 11, 12),
        (21, 21.5, 22, 22.5),
    ]

    titles = [
        r"Hot component T",
        r"Hot component N(H$_{2}$)",
        r"Cold component T",
        r"Cold component N(H$_{2}$)",
    ]

    c_labels = [
        r"Temperature (K)",
        r"log$_{10}$[Column density ($H_{2}$/cm$^2$)]",
        r"Temperature (K)",
        r"log$_{10}$[Column density ($H_{2}$/cm$^2$)]",
    ]

    w = WCS(hdr)
    size = np.array((14.5, 12))*0.6
    plt.figure(figsize=size)
    for i in range(4):
        ax = plt.subplot(221+i, projection=w)
        lo, hi = limits[i]
        plt.imshow(frames_to_plot[i], origin='lower',
            vmin=lo, vmax=hi,
            cmap='inferno')
        if i >= 2:
            plt.xlabel("Right Ascension", fontproperties=font_serif_small)
        if i % 2 == 0:
            plt.ylabel("Declination", fontproperties=font_serif_small)
        if i == 2 or i == 3:
            plt.xlim([75, 775])
            plt.ylim([150, 850])
        plt.title(titles[i], fontproperties=font_serif_large)
        c = plt.colorbar(pad=0.0, ticks=list(ticks[i]))
        c.ax.set_yticklabels([str(x) for x in ticks[i]])
        for t in c.ax.get_yticklabels():
            t.set_fontproperties(font_serif_small)
            t.set_rotation(-90.)
        c.set_label(c_labels[i],
            fontproperties=font_serif_small, rotation=-90.,
            labelpad=18)
    plt.subplots_adjust(
        top=0.95,
        bottom=0.09,
        left=0.1,
        right=.95,
        hspace=0.2,
        wspace=0.2
    )
    plt.savefig("/home/rkarim/Desktop/syp_figure_3p.png")

def make_T3p_plot():
    frames_to_plot = (1, 2, 5)
    frame3p_data = []
    with fits.open(target_3p_fn) as hdul:
        hdr = hdul[1].header
        for frame in frames_to_plot:
            frame3p_data.append(hdul[frame].data)

    frames_to_plot = (1,)
    frame2p_data = []
    with fits.open(target_2p_fn) as hdul:
        for frame in frames_to_plot:
            frame2p_data.append(hdul[frame].data)
    Th = frame2p_data[0].copy()
    Th[mask] = frame3p_data[2][mask]
    Tc = frame3p_data[0].copy()
    Tc[~mask] = np.nan
    dTc = frame3p_data[1].copy()
    dTc[~mask] = np.nan

    frames_to_plot = [Th, Tc, dTc]

    limits = [
        (13.5, 17.5), # Th
        (8.5, 12.5), # Tc
        (0, 1), # dTc
    ]

    ticks = [
        (14, 15, 16, 17),
        (9, 10, 11, 12),
        (0, 1),
    ]

    titles = [
        r"Hot component T",
        r"Cold component T",
        r"Cold component T uncertainty",
    ]

    c_labels = [
        r"Temperature (K)",
        r"Temperature (K)",
        r"Temperature (K)",
    ]

    w = WCS(hdr)
    plt.figure(figsize=(16, 13))
    for i in range(3):
        ax = plt.subplot(221+i, projection=w)
        lo, hi = limits[i]
        plt.imshow(frames_to_plot[i], origin='lower',
            vmin=lo, vmax=hi,
            cmap='inferno')
        plt.xlabel("Right Ascension", fontproperties=font_serif_small)
        plt.ylabel("Declination", fontproperties=font_serif_small)
        if i > 0:
            plt.xlim([75, 775])
            plt.ylim([150, 850])
        plt.title(titles[i], fontproperties=font_serif_large)
        c = plt.colorbar(pad=0.0, ticks=list(ticks[i]))
        c.ax.set_yticklabels([str(x) for x in ticks[i]])
        for t in c.ax.get_yticklabels():
            t.set_fontproperties(font_serif_small)
            t.set_rotation(-90.)
        c.set_label(c_labels[i],
            fontproperties=font_serif_small, rotation=-90.,
            labelpad=18)
    # plt.subplots_adjust(
    #     top=0.92,
    #     bottom=0.08,
    #     left=0.0,
    #     right=1,
    #     hspace=0.2,
    #     wspace=0.0
    # )
    plt.subplots_adjust(
        top=0.97,
        bottom=0.05,
        left=0.0,
        right=1,
        hspace=0.15,
        wspace=-0.0
    )
    plt.savefig("/home/rkarim/Desktop/syp_figure_T3p.png")

def make_N3p_plot():
    frames_to_plot = (3, 4, 7, 8)
    frame3p_data = []
    with fits.open(target_3p_fn) as hdul:
        hdr = hdul[1].header
        for frame in frames_to_plot:
            frame3p_data.append(hdul[frame].data)

    frames_to_plot = (3,)
    frame2p_data = []
    with fits.open(target_2p_fn) as hdul:
        for frame in frames_to_plot:
            frame2p_data.append(hdul[frame].data)
    Nh = frame2p_data[0].copy()
    Nh[mask] = frame3p_data[2][mask]
    Nh = np.log10(Nh)
    Nc = frame3p_data[0].copy()
    Nc[~mask] = np.nan
    Nc = np.log10(Nc)
    dNc = frame3p_data[1].copy()
    dNc[~mask] = np.nan
    dNc = np.log10(dNc)
    dNh = frame3p_data[3].copy()
    dNh[~mask] = np.nan
    dNh = np.log10(dNh)
    frames_to_plot = [Nh, dNh, Nc, dNc]

    limits = [
        (20.3, 22.), # Nh
        (19.5, 21.), # dNh
        (20.8, 22.7), # Nc
        (19.5, 22.), # dNc
    ]

    ticks = [
        (20.5, 21, 21.5, 22),
        (19.5, 20, 20.5, 21),
        (21, 21.5, 22, 22.5),
        (19.5, 20, 20.5, 21, 21.5, 22),
    ]

    titles = [
        r"Hot component N(H$_{2}$)",
        r"Uncertainty in the hot component N(H$_{2}$)",
        r"Cold component N(H$_{2}$)",
        r"Uncertainty in the cold component N(H$_{2}$)",
    ]

    c_labels = [
        r"log$_{10}$[Column density ($H_{2}$/cm$^2$)]",
        r"log$_{10}$[Column density ($H_{2}$/cm$^2$)]",
        r"log$_{10}$[Column density ($H_{2}$/cm$^2$)]",
        r"log$_{10}$[Column density ($H_{2}$/cm$^2$)]",
    ]

    w = WCS(hdr)
    plt.figure(figsize=(16, 13))
    for i in range(4):
        ax = plt.subplot(221+i, projection=w)
        lo, hi = limits[i]
        plt.imshow(frames_to_plot[i], origin='lower',
            vmin=lo, vmax=hi,
            cmap='inferno')
        plt.xlabel("Right Ascension", fontproperties=font_serif_small)
        plt.ylabel("Declination", fontproperties=font_serif_small)
        if i > 0:
            plt.xlim([75, 775])
            plt.ylim([150, 850])
        plt.title(titles[i], fontproperties=font_serif_large)
        c = plt.colorbar(pad=0.0, ticks=list(ticks[i]))
        c.ax.set_yticklabels([str(x) for x in ticks[i]])
        for t in c.ax.get_yticklabels():
            t.set_fontproperties(font_serif_small)
            t.set_rotation(-90.)
        c.set_label(c_labels[i],
            fontproperties=font_serif_small, rotation=-90.,
            labelpad=18)
    plt.subplots_adjust(
        top=0.97,
        bottom=0.05,
        left=0.0,
        right=1,
        hspace=0.15,
        wspace=-0.0
    )
    plt.savefig("/home/rkarim/Desktop/syp_figure_N3p.png")


def make_Xs_plot():
    Xs_3p, hdr = fits.getdata(target_3p_fn, 9, header=True)
    Xs_2p = fits.getdata(target_2p_fn, 5, header=False)
    Xs_3p[~mask] = np.nan
    Xs_2p[~mask] = np.nan
    lin_diff = -(Xs_3p - Xs_2p)
    diff = np.log10(lin_diff)
    diff[lin_diff < 0] = -5
    tick_positions = [0.5, 1, 2, 4, 8]
    ticks = [np.log10(x) for x in  tick_positions]
    plt.figure(figsize=(10.5, 10))
    plt.subplot(projection=WCS(hdr))
    plt.imshow(diff,
        origin='lower', vmin=np.log10(0.5), vmax=np.log10(8),
        cmap='inferno')
    plt.xlabel("Right Ascension", fontproperties=font_serif_small)
    plt.ylabel("Declination", fontproperties=font_serif_small)
    plt.title(r"Reduced $\chi^{2}$ improvement through filaments",
        fontproperties=font_serif_large)
    plt.xlim([75, 775])
    plt.ylim([150, 850])
    c = plt.colorbar(pad=0.0, ticks=list(ticks))
    c.ax.set_yticklabels([str(x) for x in tick_positions])
    for t in c.ax.get_yticklabels():
        t.set_fontproperties(font_serif_small)
        t.set_rotation(-90.)
    c.set_label(r"Goodness of fit improvement ($\chi^{2}$/DoF)",
        fontproperties=font_serif_small, rotation=-90.,
        labelpad=18)
    plt.subplots_adjust(
        top=0.91,
        bottom=0.08,
        left=0.08,
        right=1.02,
        hspace=0.0,
        wspace=0.0
    )
    plt.savefig("/home/rkarim/Desktop/syp_figure_Xs3pdiff.png")

def make_diff_plot():
    frames_to_plot = (1, 3)
    frame3p_data = []
    with fits.open(target_3p_fn) as hdul:
        hdr = hdul[1].header
        for frame in frames_to_plot:
            frame3p_data.append(hdul[frame].data)

    frames_to_plot = (1, 3)
    frame2p_data = []
    with fits.open(target_2p_fn) as hdul:
        for frame in frames_to_plot:
            frame2p_data.append(hdul[frame].data)
    diffT = frame3p_data[0] - frame2p_data[0]
    diffT[~mask] = np.nan
    diffN = np.log10(frame3p_data[1]) - np.log10(frame2p_data[1])
    diffN[~mask] = np.nan
    diffs = [diffT, diffN]
    limits = [
        (-4.5, 0),
        (-0.2, 0.2),
    ]
    w = WCS(hdr)
    plt.figure(figsize=(15, 6))
    for i in range(2):
        plt.subplot(121+i, projection=w)
        lo, hi = limits[i]
        plt.imshow(diffs[i], origin='lower', vmin=lo, vmax=hi)
        plt.colorbar()
    plt.show()


def make_herschel_plot():
    frames_to_plot = (15, 17, 19, 21)
    frame3p_data = []
    with fits.open(target_3p_fn) as hdul:
        hdr = hdul[1].header
        for frame in frames_to_plot:
            frame3p_data.append(hdul[frame].data)
    frames_to_plot = frame3p_data
    limits = [
        (0, 180),
        (0, 180),
        (0, 100),
        (0, 50),
    ]
    titles = ["160", "250", "350", "500"]
    units = "MJy/sr"
    ticks = [
        (0, 100, 180),
        (0, 100, 180),
        (0, 100),
        (0, 25),
    ]
    w = WCS(hdr)
    plt.figure(figsize=(16, 12.5))
    for i in range(4):
        ax = plt.subplot(221+i, projection=w)
        lo, hi = limits[i]
        plt.imshow(frames_to_plot[i], origin='lower',
            vmin=lo, vmax=hi,
            cmap='gray')
        if i > 1:
            plt.xlabel("Right Ascension", fontproperties=font_serif_small)
        if i % 2 == 0:
            plt.ylabel("Declination", fontproperties=font_serif_small)
        plt.title("{} micron".format(titles[i]), fontproperties=font_serif_large)
        c = plt.colorbar(pad=0.0, ticks=list(ticks[i]))
        c.ax.set_yticklabels([str(x) for x in ticks[i]])
        for t in c.ax.get_yticklabels():
            t.set_fontproperties(font_serif_small)
            t.set_rotation(-90.)
        if i % 2 == 1:
            c.set_label("Observed emission ({})".format(units),
                fontproperties=font_serif_small, rotation=-90.,
                labelpad=18)
    plt.subplots_adjust(
        top=0.97,
        bottom=0.05,
        left=0.0,
        right=1,
        hspace=0.15,
        wspace=-0.0
    )
    plt.savefig("/home/rkarim/Desktop/syp_figure_4bands.png")


make_2p_plot()

### GENERATING AND WRITING MASK TO FITS
def write_mask():
    data = fits.getdata(local_test, 1)
    nanmask = ~np.isnan(data)
    data, hdr = fits.getdata(b15n1e20, 1, header=True)
    mask = (data>11)
    mask = boolis.get_mask(mask, n=2, min_size=100, dilation=0)
    mask = boolis.fill_inwards(mask, nanmask)
    w = WCS(hdr)

    h['COMMENT'] = "BEST FILAMENT MASK TO DATE (April 23 2019)"
    h = fits.Header()
    h['CREATOR'] = ("Ramsey: {}".format(str(__file__)), "FITS file creator")
    h['OBJECT'] = ("per_04_nom", "Target name")
    h['DATE'] = (datetime.now(timezone.utc).astimezone().isoformat(),
    	"File creation date")
    h.update(w.to_header())
    fits.writeto("../filament_mask_syp.fits", mask.astype(int), header=h)
    print("WRITTEN")
