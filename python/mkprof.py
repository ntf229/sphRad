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


band_names = np.array(['FUV', 'NUV',
                       'u', 'g', 'r', 'i', 'z',
                       '2MASS_J', '2MASS_H', '2MASS_KS',
                       'W1', 'W2', 'W3', 'W4',
                       'PACS70', 'PACS100', 'PACS160',
                       'SPIRE250', 'SPIRE350', 'SPIRE500'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Measure profile of NIHAO galaxy')
    parser.add_argument('-g', '--galaxy', dest='galaxy',
                        type=str, help='name of galaxy')
    parser.add_argument('-b', '--band', dest='band',
                        type=str, help='name of band')
    parser.add_argument('-p', '--path', dest='path',
                        type=str, help='path to galaxy')
    parser.add_argument('-s', '--save', dest='save',
                        type=str, help='path to save results')
    args = parser.parse_args()
    band = args.band
    galaxy = args.galaxy
    path = args.path
    save = args.save
    # Use this to get to the right galaxy directory
    #os.chdir(galaxy)
    #os.chdir(path)
    # Find index of band
    iband = np.where(band_names == band)[0][0]
    # Make a FITS version of the band's image
    #images = np.load('phot.npy')
    images_file = fits.open(path+'/sph_broadband_total.fits')
    images = np.asarray(images_file[0].data) # cube of broadbands in MJy/sr
    image = images[iband, :, :]
    #minFlux = np.min(image[np.nonzero(image)])
    #minFlux = np.percentile(image[np.nonzero(image)], 10) # 10th percentile 
    randScale = np.mean(image[np.nonzero(image)]) / 10
    skyLevel = randScale * 2
    randArray = np.random.randint(1, high=100, size=image.shape) / 100 # random samples between 0 and 1 (excluding 0)
    image = image + skyLevel + (randScale * randArray) # add background and noise to image for FFT computation
    #image = image + (np.min(image[np.nonzero(image)])/1000) # add small flux to all pixels to prevent errors from flux=0 pixels
    #image = scipy.ndimage.gaussian_filter(image, 2) # apply gaussian filter with std=3 
    #outfile = 'galaxy-{b}.fits'.format(b=band)
    outfile = save+'galaxy-'+band+'.fits'
    fitsio.write(outfile, image, clobber=True)
    # Write autoprof's config file
    #cfgfile = 'prof_config_{b}.py'.format(b=band)
    cfgfile = save+'prof_config_'+band+'.py'
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
ap_set_psf = 2.5
ap_fit_limit = 6
ap_isofit_iterstopnochange = 5
ap_centeringring = 50
""".format(s=save,b=band))
    fp.close()
    # Run autoprof
    origDir = os.getcwd()
    os.chdir(save)
    subprocess.run(['autoprof.py', cfgfile]) 
    os.chdir(origDir)
    # Read autoprof output
    prof = ascii.read(save+'galaxy-'+band+'.prof', comment='^#.*')
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
    fp = open(save+galaxy+'-'+band+'-ba.csv', "w")
    fp.write("""band, ba, pa
{b}, {ba}, {pa}
""".format(b=band, ba=ba, pa=pa))
    fp.close()
