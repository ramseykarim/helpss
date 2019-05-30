import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cst
import mpy_utils as mpu
from Dust import Dust
from Greybody import Greybody
from Instrument import Instrument, get_Herschel
import mantipyfit as mpfit
import corner
import pickle

"""
This is where I will run science code on my laptop!
"""

def main():
    mtest_grid_to_single_pixel()


def desktop_main():
    mtest_3dgrid_to_single_pixel()

"""
Scripts below here
"""
from computer_config import manticore_soln_2p, manticore_soln_3p, mask_fn

def mtest_2pixel_scatter():
    pij1 = 478-1, 376-1 # looks ok
    pij2 = 479-1, 376-1 # low Nh
    pix_dict = {"ok": pij1, "lowNh": pij2}
    pix_infos = {}
    for k in pix_dict:
        pix_infos[k] = [
            mpu.get_manticore_info(x, *pix_dict[k])
            for x in
            [manticore_soln_2p, manticore_soln_3p]
        ]
    plt.figure()
    linestyles = ['-', '--']
    markerstyles = ['>', 'x']
    dusts = [Dust(beta=2.10), Dust(beta=1.80)]
    nu_range = np.exp(np.linspace(np.log(cst.c/(1500*1e-6)), np.log(cst.c/(50*1e-6)), 100))
    colors = iter(['blue', 'orange', 'red'])
    herschel = get_Herschel()
    for k in pix_infos:
        color = next(colors)
        for i, info_dict in enumerate(pix_infos[k]):
            Ts = [info_dict[x] for x in ("Th", "Tc") if x in info_dict]
            Ns = [info_dict[x] for x in ("Nh", "Nc") if x in info_dict]
            ds = dusts[1-i:]
            cloud = Greybody(Ts, Ns, ds)
            plt.plot(mpu.f_hz_micron(nu_range), cloud.radiate(nu_range),
                linestyle=linestyles[i], color=color,
                label=f"{str(cloud).upper()} ({k})")
            plt.plot([x*1.01 for x in mpu.H_WL],
                [d.detect(cloud) for d in herschel],
                markerstyles[i], label=f"{str(cloud).upper()} ({k})",
                color=color, markersize=8)
            if i == 0:
                obs, err = mpu.get_obs(info_dict), mpu.get_err(info_dict)
                plt.errorbar([x*0.99 for x in mpu.H_WL], obs, yerr=err,
                    fmt='.', color=color,
                    capsize=5, markersize=8)
    plt.legend()
    plt.xscale('log')
    plt.xlabel("wavelength (micron)")
    plt.ylabel("flux (MJy/sr)")
    plt.show()

def mtest_mask_chisq():
    m = mpu.fits.getdata(mask_fn).astype(bool)
    soln = mpu.fits.open(manticore_soln_3p)
    nanmask = ~np.isnan(soln[1].data)
    plt.figure()
    plt.subplot(121)
    plt.imshow((m & nanmask), origin='lower')
    plt.subplot(122)
    plt.imshow((~m & nanmask), origin='lower')
    plt.figure()
    herschel = get_Herschel()
    dusts = [Dust(beta=2.10), Dust(beta=1.80)]
    for mask in (m, ~m):
        info_dict = mpu.get_manticore_info(soln, mask & nanmask)
        Ts = tuple(zip(*[tuple(info_dict[x]) for x in ("Th", "Tc") if x in info_dict]))
        Ns = tuple(zip(*[tuple(info_dict[x]) for x in ("Nh", "Nc") if x in info_dict]))
        observations = tuple(zip(*[tuple(x) for x in mpu.get_obs(info_dict)]))
        errors = tuple(zip(*[tuple(x) for x in mpu.get_err(info_dict)]))
        residuals = mpu.deque()
        for i in range(0, len(Ts), len(Ts)//1000):
            cloud = Greybody(Ts[i], Ns[i], dusts)
            residual = sum((d.detect(cloud) - o)**2 / (e*e) for d, o, e in zip(herschel, observations[i], errors[i]))/1
            residuals.append(residual)
        print('done', end=", ")
        residuals = np.array(residuals)
        hist = mpu.histogram(residuals, x_lim=(0, 3))
        plt.plot(*hist, '-')
    soln.close()
    print()
    plt.show()

def mtest_mask_fit_2p_manticore_agreement():
    dust = Dust(beta=1.80)
    herschel = get_Herschel()
    m = mpu.fits.getdata(mask_fn).astype(bool)
    soln = mpu.fits.open(manticore_soln_2p)
    nanmask = ~np.isnan(soln[1].data)
    info_dict = mpu.get_manticore_info(soln, m&nanmask)
    soln.close()
    Ts = tuple(zip(*[tuple(info_dict[x]) for x in ("Th", "Tc") if x in info_dict]))
    Ns = tuple(zip(*[tuple(info_dict[x]) for x in ("Nh", "Nc") if x in info_dict]))
    observations = tuple(zip(*[tuple(x) for x in mpu.get_obs(info_dict)]))
    errors = tuple(zip(*[tuple(x) for x in mpu.get_err(info_dict)]))
    residuals = mpu.deque()
    values = mpu.deque()
    l = 500
    for i in range(0, len(Ts), len(Ts)//l):
        result = mpfit.fit_source_2p(observations[i], errors[i], herschel, dust)
        residuals.append([abs(mr - pr)/mr for mr, pr in zip([Ts[i][0], np.log10(Ns[i][0])], result)])
        values.append(list(result))
    residuals = np.array(residuals)
    values = np.array(values)
    plt.figure()
    plt.subplot(121)
    plt.plot(values[:, 0], residuals[:, 0], '.')
    plt.title("T")
    plt.yscale('log')
    plt.subplot(122)
    plt.plot(values[:, 1], residuals[:, 1], '.')
    plt.title("N")
    plt.yscale('log')
    plt.show()
    print("{:.1f} calls per fit".format(mpfit.ITER['a']/l))

def mtest_mask_fit_3p_manticore_agreement():
    dust = [Dust(beta=2.10), Dust(beta=1.80)]
    herschel = get_Herschel()
    m = mpu.fits.getdata(mask_fn).astype(bool)
    with mpu.fits.open(manticore_soln_3p) as soln:
        nanmask = ~np.isnan(soln[1].data)
        info_dict = mpu.get_manticore_info(soln, m&nanmask)
    Ts = tuple(zip(*[tuple(info_dict[x]) for x in ("Th", "Tc") if x in info_dict]))
    Ns = tuple(zip(*[tuple(info_dict[x]) for x in ("Nh", "Nc") if x in info_dict]))
    observations = tuple(zip(*[tuple(x) for x in mpu.get_obs(info_dict)]))
    errors = tuple(zip(*[tuple(x) for x in mpu.get_err(info_dict)]))
    residuals = mpu.deque()
    values = mpu.deque()
    l = 10
    for i in range(0, len(Ts), len(Ts)//l):
        result = mpfit.fit_source_3p(observations[i], errors[i], herschel, dust, Th=15.95)
        nres = [abs(np.log10(mr) - pr) for pr, mr in zip(result[1:], Ns[i])]
        residuals.append([abs(Ts[i][1] - result[0])]+nres)
        values.append(list(result))
        print(".", end="")
    print()
    print("{:.1f} calls per fit".format(mpfit.ITER['a']/l))
    residuals = np.array(residuals)
    values = np.array(values)
    plt.figure()
    plt.subplot(131)
    plt.plot(values[:, 0], residuals[:, 0], '.')
    plt.title("Tc")
    plt.subplot(132)
    plt.plot(values[:, 1], 10**residuals[:, 1], '.')
    plt.title("Nh")
    plt.yscale('log')
    plt.subplot(133)
    plt.plot(values[:, 2], 10**residuals[:, 2], '.')
    plt.title("Nc")
    plt.yscale('log')
    plt.show()


def mtest_selected_pixels_error():
    good_list, name_list, coords, color_list = zip(*mpu.PIXELS_OF_INTEREST)
    info_dict = mpu.get_manticore_info(manticore_soln_3p, coords)
    nu_range = np.exp(np.linspace(np.log(cst.c/(1500*1e-6)), np.log(cst.c/(50*1e-6)), 100))
    wl_range = 1e6*cst.c/nu_range
    herschel = get_Herschel()
    dusts = [Dust(beta=1.80), Dust(beta=2.10)]
    err_names = ('dTc', 'dNh', 'dNc')
    plt.figure(figsize=(8, 4.5))
    plt.subplot(111)
    for i, name in enumerate(name_list):
        if i != 6:
            continue
        # plt.subplot(331 + i)
        label = name + " (" + ("good" if good_list[i] else "bad") + ")"
        print(label)
        fmt = "o" if good_list[i] else "x"
        # plot points+error, original manticore fits
        obs = [x[i] for x in mpu.get_obs(info_dict)]
        err = [x[i] for x in mpu.get_err(info_dict)]
        # err[0] = err[0]/2
        plt.errorbar(mpu.H_WL, obs, yerr=err, fmt=fmt, capsize=6,
            color=color_list[i], label=label)
        cloud = Greybody([info_dict[x][i] for x in ("Th", "Tc")],
            [info_dict[x][i] for x in ('Nh', 'Nc')], dusts
        )
        plt.plot(wl_range, cloud.radiate(nu_range), '-', color=color_list[i],
            label=str(cloud))
        Tcf, Nhf, Ncf = mpfit.fit_source_3p(obs, err,
            herschel, dusts, Th=info_dict['Th'][i])
        Ncf, Nhf = (10**x for x in (Ncf, Nhf))
        cloudf = Greybody([info_dict['Th'][i], Tcf],
            [Nhf, Ncf], dusts)
        plt.plot(wl_range, cloudf.radiate(nu_range), '--', color=color_list[i],
            label=str(cloudf))
        p_sets, p_errs = mpfit.bootstrap_errors(obs, err, herschel,
            dusts, niter=5, fit_f=mpfit.fit_source_3p,
            dof=1, Th=info_dict['Th'][i])
        manticore_errors = tuple(info_dict[x][i] for x in err_names)
        for x in zip(err_names, manticore_errors, p_errs):
            print("{}: manticore({:.2E}), python({:.2E})".format(*x))
        title = "dTc({0:.2f}|{3:.2f}) / dNh({1:.2E}|{4:.2E}) / dNc({2:.2E}|{5:.2E})".format(
            *manticore_errors, *p_errs
        )
        # plt.title(title)
        nominal = [info_dict[x][i] for x in ("Tc", "Nh", "Nc")]
        print("nominal>> Tc:{:.2f}, Nh:{:.2E}, Nc:{:.2E}".format(*nominal))
        print("fitted >> Tc:{:.2f}, Nh:{:.2E}, Nc:{:.2E}".format(Tcf, Nhf, Ncf))
        for s in p_sets:
            cloudf = Greybody([info_dict['Th'][i], s[0]], s[1:], dusts)
            plt.plot(wl_range, cloudf.radiate(nu_range), '-', alpha=0.15,
                color='grey')
            print(">>> Tc:{:5.2f}, Nh:{:.2E}, Nc:{:.2E}".format(*s))
        plt.xscale('log')
        plt.legend()
        print()
    plt.show()

def mtest_corner_boostrap_single_pixel():
    pixel_index = 6
    good, name, coords, color = mpu.PIXELS_OF_INTEREST[pixel_index]
    info_dict = mpu.get_manticore_info(manticore_soln_3p, *coords)

    herschel = get_Herschel()
    dusts = [Dust(beta=1.80), Dust(beta=2.10)]

    obs, err = mpu.get_obs(info_dict), mpu.get_err(info_dict)
    err[0] = err[0]/2
    nominal = [info_dict[x] for x in ("Tc", "Nh", "Nc") if x in info_dict]
    Tcf, Nhf, Ncf = mpfit.fit_source_3p(obs, err, herschel, dusts, Th=info_dict['Th'])
    print("nominal>> Tc:{:5.2f}, Nh:{:.2E}, Nc:{:.2E}".format(*nominal))
    print("fitted >> Tc:{:5.2f}, Nh:{:.2E}, Nc:{:.2E}".format(Tcf, 10**Nhf, 10**Ncf))
    p_sets, p_errs = mpfit.bootstrap_errors(obs, err, herschel,
        dusts, niter=10, fit_f=mpfit.fit_source_3p, dof=1)
    for s in p_sets:
        print(">>> Tc:{:5.2f}, Nh:{:.2E}, Nc:{:.2E}".format(*s))
    return
    params = list(zip(*p_sets))
    Tcs, Nhs, Ncs = map(np.array, params)
    Nhs, Ncs = np.log10(Nhs), np.log10(Ncs)
    nominal[1] = np.log10(nominal[1])
    nominal[2] = np.log10(nominal[2])
    print("log:")
    print("nominal>> Tc:{:5.2f}, Nh:{:5.2f}, Nc:{:5.2f}".format(*nominal))
    print("fitted >> Tc:{:5.2f}, Nh:{:5.2f}, Nc:{:5.2f}".format(Tcf, Nhf, Ncf))
    labels = ['Tc', 'log(Nh)', 'log(Nc)']
    params = np.stack([Tcs, Nhs, Ncs], axis=1)
    # fig = corner.corner(params, labels=labels, truths=nominal,)
#        range=[(0, 15), (18, 22), (18, 23.5)])
    plt.show()
    return


def mtest_grid_to_single_pixel():
    pixel_index = 6
    good, name, coords, color = mpu.PIXELS_OF_INTEREST[pixel_index]
    info_dict = mpu.get_manticore_info(manticore_soln_2p, *coords)
    # nu_range = np.exp(np.linspace(np.log(cst.c/(1500*1e-6)), np.log(cst.c/(50*1e-6)), 100))
    # wl_range = 1e6*cst.c/nu_range
    herschel = get_Herschel()
    dof = 2.
    dusts = Dust(beta=1.80) # [Dust(beta=1.80), Dust(beta=2.10)]
    obs, err = mpu.get_obs(info_dict), mpu.get_err(info_dict)
    nominal = [info_dict[x] for x in ('Tc', 'Nc')]
    nominal[1] = np.log10(nominal[1])
    print("T:{:5.2f}, N:{:5.2f}".format(*nominal))
    Tclim, Nclim = (14, 17), (21.2, 21.8)
    ex = (*Tclim, *Nclim)
    aspect = (Tclim[1] - Tclim[0]) / (Nclim[1] - Nclim[0])
    Tcrange = np.linspace(*Tclim, 40)
    Ncrange = np.linspace(*Nclim, 40)
    Tcgrid, Ncgrid = np.meshgrid(Tcrange, Ncrange, indexing='xy')
    gofgrid = np.empty(Tcgrid.size)
    for i, pvec in enumerate(zip(Tcgrid.ravel(), Ncgrid.ravel())):
        gof = mpfit.goodness_of_fit_f_2p(pvec, dusts, obs, err, herschel, dof)
        gofgrid[i] = gof
    gofgrid = np.log10(gofgrid.reshape(Tcgrid.shape))
    plt.imshow(gofgrid, origin='lower', extent=ex, aspect=aspect)
    plt.show()

def mtest_3dgrid_to_single_pixel():
    pixel_index = 6
    good, name, coords, color = mpu.PIXELS_OF_INTEREST[pixel_index]
    info_dict = mpu.get_manticore_info(manticore_soln_3p, *coords)
    # nu_range = np.exp(np.linspace(np.log(cst.c/(1500*1e-6)), np.log(cst.c/(50*1e-6)), 100))
    # wl_range = 1e6*cst.c/nu_range
    herschel = get_Herschel()
    dof = 1.
    dusts = [Dust(beta=1.80), Dust(beta=2.10)]
    obs, err = mpu.get_obs(info_dict), mpu.get_err(info_dict)
    nominal = [info_dict[x] for x in ('Tc', 'Nh', 'Nc')]
    Th = info_dict['Th']
    nominal[1] = np.log10(nominal[1])
    nominal[2] = np.log10(nominal[2])
    print("Tc:{:5.2f}, Nh:{:5.2f}, Nc:{:5.2f}".format(*nominal))
    Tclim, Nclim = (10, 14), (20., 23.)
    Nhlim = (20., 22.)
    # ex = (*Tclim, *Nclim)
    # aspect = (Tclim[1] - Tclim[0]) / (Nclim[1] - Nclim[0])
    Tcrange = np.arange(*Tclim, 0.03)
    Nhrange = np.arange(*Nhlim, 0.03)
    Ncrange = np.arange(*Nclim, 0.03)
    Tcgrid, Nhgrid, Ncgrid = np.meshgrid(Tcrange, Nhrange, Ncrange, indexing='ij')
    print(Tcgrid.shape, Tcgrid.size)
    if int(input("Ready? Type '1'")) != 1:
        print("quitting...")
        return
    gofgrid = np.empty(Tcgrid.size)
    for i, pvec in enumerate(zip(Tcgrid.ravel(), Nhgrid.ravel(), Ncgrid.ravel())):
        gof = mpfit.goodness_of_fit_f_3p(pvec, dusts, obs, err, herschel, Th, dof)
        gofgrid[i] = gof
    gofgrid = np.log10(gofgrid.reshape(Tcgrid.shape))
    # plt.imshow(gofgrid, origin='lower', extent=ex, aspect=aspect)
    # plt.show()
    with open("grid_file.pkl", 'wb') as pfl:
        pickle.dump(gofgrid, pfl)
    print("Written and finished")
    return

if __name__ == "__main__":
    main()