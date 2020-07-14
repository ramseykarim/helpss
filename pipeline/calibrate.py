#!/usr/bin/env python
"""
Wrapper script for the PACS calibration code in scripts/
More automatic version of example.py; uses command line arguments.

To use this script from any directory, run this shell command from the
directory containing this file:
~$ PATH=$(pwd):$PATH
And then you can leave the directory and run commands like
~$ calibrate.py ./processed

Created: July 9, 2020
"""
__author__ = "Ramsey Karim"

import sys
import os
import argparse

from scripts import calc_offset, modify_fits
# calc_offset, modify_fits = None, None


def get_data_path():
    """
    Pull the data directory from the first command line argument
    :returns: absolute path to the PACS data folder for reading and writing,
        dictionary with other useful keys (like TEST, True if we shouldn't
        do any I/O)
    """
    parser = argparse.ArgumentParser(description="Command line tool to zero-point calibrate the PACS data and add systematic uncertainties the error maps.")
    parser.add_argument('directory', type=str, nargs='?', default='./', help="directory containing the Herschel PACS/SPIRE data (default: <current directory> ).")
    parser.add_argument('--test', action='store_true', help="print out the actions that would be taken, but do not execute any of them. No I/O at all.")
    args = parser.parse_args()
    data_path = args.directory

    if not os.path.isabs(data_path):
        data_path = os.path.abspath(data_path)
    if not os.path.isdir(data_path):
        raise RuntimeError(f"Invalid directory: {data_path}")

    other_kwargs = {'test': args.test}
    return data_path, other_kwargs


uncertainties = {"PACS70um": 6., "PACS160um": 8., "SPIRE250um": 5.5, "SPIRE350um": 5.5, "SPIRE500um": 5.5}

if __name__ == "__main__":
    data_path, other_kwargs = get_data_path()
    # Can list 70 and 160 here if you wanted
    band_stubs = ["PACS160um"]
    modified_flux_files = {}
    for band_stub in band_stubs:
        pacs_flux_filename = f"{data_path}/{band_stub}-image-remapped.fits"
        if not other_kwargs['test']:
            model = calc_offset.GNILCModel(pacs_flux_filename, target_bandpass=band_stub)
            derived_offset = model.get_offset(full_diagnostic=True)
        else:
            print("--- model calculation ---")
            derived_offset = -99.99
        print(f"Working on {pacs_flux_filename}")
        assigned_offset = int(input("Derived offset is {:.2f}. Assign: ".format(derived_offset)))
        pacs_flux_filename = f"{data_path}/{band_stub}-image-remapped-conv.fits"
        calibrated_pacs_flux_filename = modify_fits.add_offset(int(round(assigned_offset)), pacs_flux_filename, savename=data_path, test=other_kwargs['test'])
        modified_flux_files[band_stub] = calibrated_pacs_flux_filename
        print("written")

    for band_stub in uncertainties:
        # Consistency between image-conv and error-conv
        if band_stub in modified_flux_files:
            flux_filename = modified_flux_files[band_stub]
        else:
            flux_filename = f"{data_path}/{band_stub}-image-remapped-conv.fits"
        error_filename = f"{data_path}/{band_stub}-error-remapped-conv.fits"
        if not os.path.isfile(flux_filename):
            print(f"{band_stub} not found, skipping")
            continue
        """
        SPIRE from https://www.cosmos.esa.int/documents/12133/1035800/QUICK-START+GUIDE+TO+HERSCHEL-SPIRE
        PACS from https://www.cosmos.esa.int/documents/12133/996891/PACS+Photometer+Quick+Start+Guide
        """
        if not other_kwargs['test']:
            modify_fits.add_systematic_error(uncertainties[band_stub]*0.01, error_filename, flux_filename,
                savename=data_path)
        else:
            print(f"--- writing {uncertainties[band_stub]} % error to {error_filename} based on {flux_filename} ---")
