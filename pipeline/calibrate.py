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

bands = {70: "PACS70um", 160: "PACS160um"}
uncertainties = {"PACS70um": 6., "PACS160um": 8., "SPIRE250um": 5.5, "SPIRE350um": 5.5, "SPIRE500um": 5.5}


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
    parser.add_argument('--band', type=int, nargs='*', default=[160], help="select the PACS bands to be calibrated. Use integer wavelengths in microns.")
    parser.add_argument('--assign', action='store_true', help="skip straight to the offset assignment. No calculations or diagnostic figures will be produced.")
    args = parser.parse_args()
    data_path = args.directory

    if not os.path.isabs(data_path):
        data_path = os.path.abspath(data_path)
    if not os.path.isdir(data_path):
        raise RuntimeError(f"Invalid directory: {data_path}")

    other_kwargs = {'test': args.test, 'bands': [bands[k] for k in args.band], 'assign': args.assign}
    return data_path, other_kwargs


def safe_cast_integer(x):
    """
    Check if this is an integer, and if so, return the integer
    :param x: the thing to check
    :returns: False if cannot cast to int, or the resulting integer if it can
    """
    try:
        return int(x)
    except:
        return False


def get_assigned_offset(band_stub, derived_offset):
    """
    Manage input loop and return an integer offset from user
    :param band_stub: the string bandpass name to reference in the prompt
    :param derived_offset: the float offset derived by the automatic procedure
    :returns: the int assigned offset. Or, None if the file writing should
        be skipped.
    """
    first_msg = f"Derived {band_stub} offset is {derived_offset:.2f}. Assign: "
    next_msg = lambda s: f"Having trouble interpreting your response ({s}) as an integer. Enter 'q' to quit. Assign: "
    quit_commands = ['exit', '', 'q', 'x', 'quit']
    response = input(first_msg)
    while (response.lower() not in quit_commands) and (safe_cast_integer(response) is False):
        # This response is NEITHER a quit command NOR an integer. Try again
        response = input(next_msg(response))
    # Passed through the while loop under one of the conditions.
    # Find out which one
    if response.lower() in quit_commands:
        # Return None, indicating that we should not assign anything
        return None
    else:
        # Must be an integer
        return safe_cast_integer(response)


if __name__ == "__main__":
    STOP = False
    data_path, other_kwargs = get_data_path()
    # Can list 70 and 160 here if you wanted
    band_stubs = other_kwargs['bands']
    modified_flux_files = {}
    for band_stub in band_stubs:
        # Line up all the filenames
        # Use the NATIVE RESOLUTION (but regridded) PACS map
        pacs_flux_filename = f"{data_path}/{band_stub}-image-remapped.fits"
        # SPIRE flux map filenames for the masks in calc_offset
        spire_filenames = {
            "spire250_filename": f"{data_path}/SPIRE250um-image-remapped.fits",
            "spire500_filename": f"{data_path}/SPIRE500um-image-remapped.fits",
        }
        print(f"Working on {pacs_flux_filename}")
        # Handle the possibility of debug runs that don't need calculation
        if other_kwargs['test']:
            print("--- model calculation ---")
            derived_offset = -99.99
        elif other_kwargs['assign']:
            print("Skipping to assignment. No calculations will be made.")
            derived_offset = -99.99
        else:
            model = calc_offset.GNILCModel(pacs_flux_filename, target_bandpass=band_stub, **spire_filenames)
            derived_offset = model.get_offset(full_diagnostic=True)
        # Use the input-handling function to get the assigned offset
        assigned_offset = get_assigned_offset(band_stub, derived_offset)
        if assigned_offset is None:
            # Do other bands, if present, but do not write errors
            STOP = True
            continue
        # Write out the CONVOLVED flux map plus assigned offset
        # The offset should NOT have a resolution dependence
        pacs_flux_filename = f"{data_path}/{band_stub}-image-remapped-conv.fits"
        calibrated_pacs_flux_filename = modify_fits.add_offset(int(round(assigned_offset)), pacs_flux_filename, savename=data_path, test=other_kwargs['test'])
        modified_flux_files[band_stub] = calibrated_pacs_flux_filename
        print("written")

    if STOP:
        sys.exit()

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
