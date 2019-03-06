import sys
from astropy.io import fits

"""
This quickly edits FITS headers for the remap/conv/crop images/errors
It adds OBJECT (appears to be a problem)

First cmd line arg is the filename up to image/error (it'll do both)
Second is the new OBJECT value

It will not overwrite if OBJECT already exists, but will instead tell you the value.
"""

rcc = "-remapped-conv.fits"
img = "image"
err = "error"

# GET FILE
try:
    data, header = fits.getdata(sys.argv[1]+img+rcc, header=True)
except IndexError:
    print("need fits name")
    sys.exit()
except FileNotFoundError:
    print("File not found: %s" % sys.argv[1]+img+rcc)
    sys.exit()

# FIX OBJECT
try:
    print("HEADER[OBJECT] ALREADY CONTAINS: ", header['OBJECT'])
    try:
        data2, header2 = fits.getdata(sys.argv[1]+err+rcc, header=True)
        print("error has it too: ", header2['OBJECT'])
    except KeyError:
        print("image is ok but error is not??")
except KeyError:
    try:
        i = sys.argv.index('-o')
        header['OBJECT'] = sys.argv[i+1]
        fits.writeto(sys.argv[1]+img+rcc, data, header, overwrite=True)
        data2, header2 = fits.getdata(sys.argv[1]+err+rcc, header=True)
        header2['OBJECT'] = sys.argv[i+1]
        fits.writeto(sys.argv[1]+err+rcc, data2, header2, overwrite=True)
    except ValueError:
        print("Does not have OBJECT header! Put desired value after -o flag")
    except IndexError:
        print("whoops, did you use flags wrong?")
        sys.exit()

# FIX HISTORY
try:
    print("HEADER[HISTORY] ALREADY CONTAINS: ", header['HISTORY'])
    try:
        print("error has it too: ", header2['HISTORY'])
    except NameError:
        data2, header2 = fits.getdata(sys.argv[1]+err+rcc, header=True)
        print("error has it too: ", header2['HISTORY'])
    except KeyError:
        print("image is ok but error is not??")
except KeyError:
    try:
        i = sys.argv.index('-h')
        header['HISTORY'] = sys.argv[i+1]
        fits.writeto(sys.argv[1]+img+rcc, data, header, overwrite=True)
        try:
            header2['HISTORY'] = sys.argv[i+1]
        except NameError:
            data2, header2 = fits.getdata(sys.argv[1]+err+rcc, header=True)
            header2['HISTORY'] = sys.argv[i+1]
        fits.writeto(sys.argv[1]+err+rcc, data2, header2, overwrite=True)
    except ValueError:
        print("Does not have HISTORY header! Put desired value after -h flag")
    except IndexError:
        print("whoops, did you use flags wrong?")
        sys.exit()
