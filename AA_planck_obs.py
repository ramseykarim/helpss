import numpy as np
from scipy.constants import k, h, c
from astropy.io import fits
import AA_regrid_healpy as rhp
import sys


def sizeof_fmt(num, suffix='B'):
    # Taken directly from stackoverflow
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)

FN_T = "COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.00.fits"
FN_tau = "COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.00.fits"
FN_beta = "COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.00.fits"
# Above filenames are obsolete: this work is done on my laptop (For now..)

component_stubs = ("Temperature", "Opacity", "SpectralIndex")
def gen_component_fn(field_stub, band_stub, component_stub):
    return "dustModel_%s_%sgrid_%s.fits" % (field_stub, band_stub, component_stub)

REF_F = 353*1e9  # 353 GHz is the reference frequency for Planck opacity

angstroms_to_Hz = lambda x: c / (x * 1e-10)

hb_dr = "/n/sgraraid/filaments/data/TEST4/regridding_stuff/Herschel_bands/"
p_RIMO = "/n/sgraraid/filaments/data/TEST4/regridding_stuff/HFI_stuff/HFI_RIMO_R3.00.fits"
spire_stub = lambda x: hb_dr+"Herschel_SPIRE.P"+x+"W_ext.dat"
kevin_bp_fn = lambda b: "{}_fromManticore.dat".format(b)
bandpass_files = {
    "PACS160um": hb_dr+kevin_bp_fn("PACS160um"),
    "SPIRE250um": hb_dr+kevin_bp_fn("SPIRE250um"),
    "SPIRE350um": hb_dr+kevin_bp_fn("SPIRE350um"),
    "SPIRE500um": hb_dr+kevin_bp_fn("SPIRE500um"),
    "F100": 3,
    "F143": 4,
    "F217": 5,
    "F353": 6,
    "F545": 7,
    "F857": 8,
}

bandpass_frequencies = {
    # Input in Angstroms, converted immediately
    "PACS160um": angstroms_to_Hz(1539451.3),
    "SPIRE250um": angstroms_to_Hz(2471245.1),
    "SPIRE350um": angstroms_to_Hz(3467180.4),
    "SPIRE500um": angstroms_to_Hz(4961067.7),
    "F100": 100*1e9,
    "F143": 143*1e9,
    "F217": 217*1e9,
    "F353": 353*1e9,
    "F545": 545*1e9,
    "F857": 857*1e9,
}

pfreq_lims = (1e10, 2e12) # Hz
limit_planck_frequency = lambda arr, f_arr: arr[(f_arr > pfreq_lims[0]) & (f_arr < pfreq_lims[1])]

def get_bandpass_data(stub):
    """
    Get photometry weighting function for either Herschel or Planck HFI
    Returns tuple(frequency array in Hz, weight array)
    """
    if stub[0] == "F":
        with fits.open(p_RIMO) as hdul:
            i = bandpass_files[stub]
            assert stub == hdul[i].header['EXTNAME'][-4:]
            frequency_SI = hdul[i].data['WAVENUMBER'] * c * 1e2
            weight = hdul[i].data['TRANSMISSION']
        # Limit the frequency range of the outputs to avoid large/small floats
        weight = limit_planck_frequency(weight, frequency_SI)
        frequency_SI = frequency_SI[(frequency_SI > 1e10) & (frequency_SI < 2e12)] #limit_planck_frequency(frequency_SI, frequency_SI)
    else:
        fn = bandpass_files[stub]
        # FIXME do we need the hb_dr filename in here???
        # fixed: nah, we're good. delete this comment eventually
        bp_data = np.loadtxt(fn)
        frequency_SI, weight = bp_data[:, 0], bp_data[:, 1]
    return frequency_SI, weight


def B(T, v):
    # Planck function, SI units
    front = 2*h*(v**3) / (c*c)
    exponent = h*v / (k*T)
    return front/(np.exp(exponent) - 1)


class DustModel:
    def __init__(self, field_stub, band_stub, target_fits_data, target_fits_header, directory):
        fn_T, fn_tau, fn_beta = FN_T, FN_tau, FN_beta
        self.projection_wizard = rhp.ProjectionWrapper(target_fits_data, target_fits_header)
        self.T = self.regrid(directory + fn_T)
        self.tau = self.regrid(directory + fn_tau)
        self.beta = self.regrid(directory + fn_beta)
        self.frequency = None
        self.weight = None
        self.load_bandpass_data(band_stub)

    def regrid(self, filename):
        # Specifically for temperature, opacity, and spectral index
        # Reworked to use healpy. Option of either regrid_healpix_to_fits or using ProjectionWrapper class
        # Read regrid_healpy documentation for more information
        # Open a HEALPix with open_healpix(filename, nest=nest)
        print("Loading in %s...", end="")
        component_image = self.projection_wizard.generate_image(rhp.open_healpix(filename))
        print("done.")
        return component_image[np.newaxis, :]

    def load_bandpass_data(self, band_stub):
        self.frequency, self.weight = get_bandpass_data(band_stub)
        # Incorporate the d_frequency integral term into weight
        df = np.empty(self.frequency.shape)
        df[:-1] = self.frequency[1:] - self.frequency[:-1]
        # Duplicate last df term to account for missing "end" frequency
        df[-1] = df[-2]
        self.weight *= df
        # Calculate the normalization for the bandpass integration. Includes color correction.
        # CC: spectral index -1
        # self.normalization = np.sum(self.weight)
        # CC: spectral index 0
        self.normalization = np.sum(self.weight * bandpass_frequencies[band_stub] / self.frequency)
        print("Band centered at %.1f GHz" % (bandpass_frequencies[band_stub]/1e9))
        # Expand axes so these are flat 3D arrays
        self.weight = self.weight[:, np.newaxis, np.newaxis]
        self.frequency = self.frequency[:, np.newaxis, np.newaxis]

    def flux_helper(self, n_start, n_end):
        # For flat (though multi-D) frequency array (SI units--Hertz)
        subset_retval = np.empty((self.frequency.shape[0], n_end-n_start, self.T.shape[2]))
        crop = lambda arr: arr[:, n_start:n_end, :]
        subset_retval[:] = B(crop(self.T), self.frequency)
        subset_retval *= crop(self.tau)
        subset_retval *= (self.frequency/REF_F)**crop(self.beta)
        subset_retval *= self.weight
        subset_retval = np.sum(subset_retval, axis=0) / self.normalization
        # Convert from SI to MJy/sr
        return subset_retval * 1e20

    def observe_planck(self, band_stub, islarge=False, row_step=None):
        """
        This is the function the user should call to create a Herschel-like Planck image
        """
        if islarge:
            # Will calculate several rows at a time to conserve memory
            if row_step is None:
                row_step = 1024
            n_rows, n_cols = self.T.shape[1], self.T.shape[2]
            n_steps = n_rows // row_step
            n_rows_leftover = n_rows % row_step
            return_value = np.zeros((n_rows, n_cols))
            print("ROWS: %d. LEFTOVER: %d. Need %d steps." % (n_rows, n_rows_leftover, n_steps))
            print("Total size: %d x %d x %d x 64 bits = %s" % (self.frequency.shape[0], n_rows, n_cols, sizeof_fmt(self.frequency.shape[0] * n_rows * n_cols * 8)))
            print("Step size: %d x %d x %d x 64 bits = %s" % (self.frequency.shape[0], row_step, n_cols, sizeof_fmt(self.frequency.shape[0] * row_step * n_cols * 8)))
            for i in range(n_steps):
                n_start, n_end = i*row_step, (i+1)*row_step
                print("-> array[0:%d, %d:%d, 0:%d]\r" % (self.frequency.shape[0], n_start, n_end, n_cols), end="")
                return_value[n_start:n_end, :] = self.flux_helper(n_start, n_end)
            n_start, n_end = n_steps*row_step, n_steps*row_step + n_rows_leftover
            print("-> array[0:%d, %d:%d, 0:%d]" % (self.frequency.shape[0], n_start, n_end, n_cols))
            return_value[n_start:n_end, :] = self.flux_helper(n_start, n_end)
        else:
            print("Calculating sky at %s band..." % band_stub, end="")
            return_value = self.flux_helper(0, self.T.shape[1])
            print("done.")
        print("Interpolating back to Herschel grid...", end="")
        return_value = self.projection_wizard.interpolate_to_target(return_value)
        print("done.")
        return return_value


if __name__ == '__main__':
    ## COMPARE Kevin's bandpass with SVOC: using Kevin's now
    band = "PACS160um"
    import matplotlib.pyplot as plt
    plt.figure(figsize=(15, 15))
    axes = (plt.subplot(221 + i) for i in range(4))
    for band in bandpass_frequencies:
        if band[0] == 'F':
            continue
        plt.sca(next(axes))
        frequency1, weight1 = get_bandpass_data(band)
        manticore_dat = np.loadtxt(hb_dr+kevin_bp_fn(band))
        weight1 = weight1 * np.max(manticore_dat[:, 1]) / np.max(weight1)
        plt.plot(manticore_dat[:, 0], manticore_dat[:, 1], '.', label='manticore')
        plt.plot(frequency1, weight1, '.', label='SVOC')
        plt.title(band)
        plt.legend()
    plt.show()
    sys.exit()

    if len(sys.argv) <= 1:
        raise InputError("Need filename of Herschel image. Will assume dust model parameters are in the same directory.")
    reference_filename = sys.argv[1]
    field_stub = sys.argv[2]
    band_stub = None
    for bs in bandpass_files:
        if bs in reference_filename:
            band_stub = bs
    if band_stub is None:
        raise KeyError("Couldn't find %s in band stubs" % reference_filename)
    region_directory = "/n/sgraraid/filaments/data/TEST4/Per/testregion1342190326JS/"
    reference_filename = region_directory + reference_filename
    with fits.open(reference_filename) as hdul:
        h_data = hdul[0].data
        h_header = hdul[0].header
    # The filename below is correct, even after the directory move
    dust_directory = "/n/sgraraid/filaments/data/TEST4/regridding_stuff/"
    sky = DustModel(field_stub, band_stub, h_data, h_header, dust_directory)
    planck350 = sky.observe_planck(band_stub)
    new_filename = region_directory+"%s-image-PLANCKfabricated.fits" % band_stub
    h_header['HISTORY'] = "THIS IS NOT FROM HERSCHEL -- it is Herschel-like, made from the Planck dust model (from Ramsey 3/4/19)"
    print("can save at \n\t%s\nusing command:\n > fits.writeto(new_filename, planck350, h_header)" % new_filename)
    fits.writeto(new_filename, planck350, h_header)
    import matplotlib.pyplot as plt
    plt.imshow(planck350, origin='lower')
    plt.show()
