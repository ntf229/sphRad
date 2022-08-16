#!/bin/env python

import os
import sys
import argparse
import subprocess
import numpy as np
import astropy.io.ascii as ascii
import fitsio
from astropy.io import fits
import scipy.ndimage

# Call this code with
#
#  mkprof -g [galaxy] -b [band]
#
# This code depends on autoprof:
#  https://autoprof.readthedocs.io/en/latest/index.html
#
# To install autoprof:
#
#  (1) pip install photutils
#  (2) Follow install instructions for autoprof
#  (3) But instead of using an alias for autoprof.py, add its directory
#      into your unix path
#  (4) Run the tests in the "test" directory; make sure it actually
#      outputs a whole bunch of files, and look at them. It claims it
#      will find errors in the output logs, but it turns out to miss
#      some.
# 
# The code below should then work when run in a directory with each
# galaxy as a subdirectory.
#
# It will output (in the subdirectory) a csv file called:
#   {galaxy}-{band}-ba.csv
# which has the band name, the b/a estimated, and the position angle
# (in degrees from the + Y-axis toward the - X-axis, i.e. it would be
# the usual astronomical definition for a North-up, East-left aligned
# image).
#
# The b/a and the PA are estimated as the average of the ellipses between
# the 50% and 90% light radii. This choice is to avoid the inner bulge or
# bar regions.


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Measure profile of NIHAO galaxy')
    parser.add_argument('-g', '--galaxy', dest='galaxy',
                        type=str, help='name of galaxy')
    parser.add_argument('-n', '--name', dest='name',
                        type=str, help='inc and az values')
    parser.add_argument('-p', '--path', dest='path',
                        type=str, help='path to galaxy')
    parser.add_argument('-s', '--save', dest='save',
                        type=str, help='path to save results')
    args = parser.parse_args()
    galaxy = args.galaxy
    path = args.path
    save = args.save
    name = args.name
    # Make a FITS version of the band's image
    images_file = fits.open(path+'/sph_'+name+'_total.fits')
    images = np.asarray(images_file[0].data) # cube of broadbands in MJy/sr
    image = images[0, :, :]
    randScale = np.mean(image[np.nonzero(image)]) / 10
    skyLevel = randScale * 2
    randArray = np.random.randint(1, high=100, size=image.shape) / 100 # random samples between 0 and 1 (excluding 0)
    image = image + skyLevel + (randScale * randArray) # add background and noise to image for FFT computation
    outfile = save+'galaxy-r.fits'
    fitsio.write(outfile, image, clobber=True)
    # Write autoprof's config file
    cfgfile = save+'prof_config_r.py'
    fp = open(cfgfile, 'w')
    fp.write("""
ap_process_mode = "image"

ap_image_file = "{s}galaxy-{b}.fits"
ap_name = "galaxy-{b}"
ap_pixscale = 1.0
ap_zeropoint = 22.5
ap_samplegeometricscale = 0.1
ap_doplot = True
ap_isoclip = False
ap_isofit_perturbscale_ellip = 0.03
ap_isofit_perturbscale_pa = 0.06
ap_isofit_iterlimitmax = 1000
ap_isofit_iterlimitmin = 100
ap_set_psf = 1
ap_fit_limit = 9
ap_isofit_iterstopnochange = 3
ap_centeringring = 50
ap_isoband_fixed = True
ap_regularize_scale = 0.1
""".format(s=save,b='r'))
    fp.close()
    # Run autoprof
    origDir = os.getcwd()
    os.chdir(save)
    subprocess.run(['autoprof.py', cfgfile]) 
    os.chdir(origDir)
    # Read autoprof output
    prof = ascii.read(save+'galaxy-r.prof', comment='^#.*')
    # Find 50% and 90% light thresholds
    mag = prof['totmag'].min()
    mag90 = mag - 2.5 * np.log10(0.9)
    mag50 = mag - 2.5 * np.log10(0.5)
    # Now find average ellipticity in that range
    inrange = ((prof['totmag'] > mag50) & (prof['totmag'] < mag90))
    ellip = prof['ellip'].mean()
    pa = prof['pa'].mean()
    # un-documented definition determined by test with fake data 
    ba = 1. - ellip
    # Output
    fp = open(save+galaxy+'-r-ba.csv', "w")
    fp.write("""band, ba, pa
{b}, {ba}, {pa}
""".format(b='r', ba=ba, pa=pa))
    fp.close()
