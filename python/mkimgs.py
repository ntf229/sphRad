#!/bin/env python

import os
import sys
import argparse
import numpy as np
import scipy.signal as signal
import matplotlib.image
import astropy.visualization as viz


band_names = np.array(['FUV', 'NUV',
                       'u', 'g', 'r', 'i', 'z',
                       '2MASS_J', '2MASS_H', '2MASS_KS',
                       'W1', 'W2', 'W3', 'W4',
                       'PACS70', 'PACS100', 'PACS160',
                       'SPIRE250', 'SPIRE350', 'SPIRE500'])

band_indx = dict()
for iband, band in enumerate(band_names):
    band_indx[band] = iband

i = dict()

i['irg'] = {'R': {'i': 1.},
            'G': {'r': 1.},
            'B': {'g': 0.7},
            'Rscale': 50.,
            'Gscale': 50.,
            'Bscale': 50.,
            'stretch': 1.,
            'Q': 10.,
            'minimum': 0.}

i['rgu'] = {'R': {'r': 1.},
            'G': {'g': 1.},
            'B': {'u': 1.},
            'Rscale': 30.,
            'Gscale': 30.,
            'Bscale': 30.,
            'stretch': 1.,
            'Q': 10.,
            'minimum': 0.}

i['W3W2W1'] = {'R': {'W3': 1.},
               'G': {'W2': 1.},
               'B': {'W1': 1.},
               'Rscale': 20.,
               'Gscale': 50.,
               'Bscale': 30.,
               'stretch': 1.,
               'Q': 10.,
               'minimum': 0.}

i['S2P1W3'] = {'R': {'SPIRE250': 1.},
               'G': {'PACS100': 1.},
               'B': {'W3': 1.},
               'Rscale': 0.01,
               'Gscale': 0.005,
               'Bscale': 0.005,
               'stretch': 1.,
               'Q': 30.,
               'minimum': 0.}

i['rNF'] = {'R': {'r': 1.},
            'G': {'NUV': 1.},
            'B': {'FUV': 1.},
            'Rscale': 5.,
            'Gscale': 40.,
            'Bscale': 40.,
            'stretch': 30.,
            'Q': 30.0,
            'minimum': 0.}

i['irgF'] = {'R': {'i': 1., 'PACS100': 0.006},
             'G': {'r': 1.},
             'B': {'g': 0.7, 'FUV': 8.},
             'Rscale': 50.,
             'Gscale': 50.,
             'Bscale': 50.,
             'stretch': 1.,
             'Q': 10.0,
             'minimum': 0.}

i['rguF'] = {'R': {'r': 1., 'PACS100': 0.006},
             'G': {'g': 1.},
             'B': {'u': 1., 'FUV': 5.},
             'Rscale': 30.,
             'Gscale': 30.,
             'Bscale': 30.,
             'stretch': 1.,
             'Q': 10.0,
             'minimum': 0.}


def rgb_image(rimage=None, gimage=None, bimage=None,
              minimum=0., stretch=1., Q=10., filename=None):

    uint8max = np.float32(np.iinfo(np.uint8).max)

    intensity = (rimage + gimage + bimage) / 3.
    intensity = intensity + 1.e-9 * intensity.max()
    sintensity = np.arcsinh(stretch * Q * (intensity - minimum)) / Q
    R = np.where(sintensity > 0, sintensity * rimage / intensity, 0.)
    G = np.where(sintensity > 0, sintensity * gimage / intensity, 0.)
    B = np.where(sintensity > 0, sintensity * bimage / intensity, 0.)

    RGB = np.stack([R, G, B], axis=2)
    maxRGB = RGB.max(axis=2)
    scale = np.where(maxRGB > 1., 1. / (maxRGB + 1.e-16), 1.)

    RGB[:, :, 0] *= scale * uint8max
    RGB[:, :, 1] *= scale * uint8max
    RGB[:, :, 2] *= scale * uint8max

    RGB = np.uint8(RGB)

    if(filename is not None):
        matplotlib.image.imsave(filename, RGB, origin='lower')

    return(RGB)


def binimage(image=None, binsize=4):
    npsf = binsize * 10 + 1
    d = np.arange(npsf, dtype=np.float32) - np.float32(npsf // 2)
    xpsf = np.outer(d, np.ones(npsf, dtype=np.float32))
    ypsf = np.outer(np.ones(npsf, dtype=np.float32), d)
    r2psf = (xpsf**2 + ypsf**2) / np.float32(binsize)**2
    psf = np.exp(- 0.5 * r2psf) / (2. * np.pi * np.float32(binsize)**2)
    image2 = signal.fftconvolve(image, psf, mode='same')
    image_out = np.zeros((image2.shape[0] // binsize,
                          image2.shape[1] // binsize), dtype=np.float32)
    for i in np.arange(binsize, dtype=int):
        for j in np.arange(binsize, dtype=int):
            image_out = image_out + image2[i::binsize, j::binsize]
    return(image_out)


def image(images=None, imagetype='irg', filename=None):

    if(imagetype not in i):
        print("No such image type {i}".format(i=imagetype))
        return

    itype = i[imagetype]

    rimage = None
    for band in itype['R']:
        scale = itype['R'][band]
        tmp_image = scale * binimage(images[band_indx[band], :, :])
        if(rimage is None):
            rimage = tmp_image
        else:
            rimage = rimage + tmp_image
    rimage = rimage * itype['Rscale']

    gimage = None
    for band in itype['G']:
        scale = itype['G'][band]
        tmp_image = scale * binimage(images[band_indx[band], :, :])
        if(gimage is None):
            gimage = tmp_image
        else:
            gimage = gimage + tmp_image
    gimage = gimage * itype['Gscale']

    bimage = None
    for band in itype['B']:
        scale = itype['B'][band]
        tmp_image = scale * binimage(images[band_indx[band], :, :])
        if(bimage is None):
            bimage = tmp_image
        else:
            bimage = bimage + tmp_image
    bimage = bimage * itype['Bscale']

    stretch = itype['stretch']
    Q = itype['Q']
    minimum = itype['minimum']
    rgb = rgb_image(rimage=rimage, gimage=gimage, bimage=bimage,
                    stretch=stretch, Q=Q, minimum=minimum,
                    filename=filename)
    return(rgb)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description='Make images of NIHAO')

    parser.add_argument('-g', '--galaxy', dest='galaxy',
                        type=str, help='name of galaxy')
    parser.add_argument('-i', '--image-type', dest='imagetype',
                        type=str, help='name of image type')

    args = parser.parse_args()
    imagetype = args.imagetype
    galaxy = args.galaxy

    # Use this to get to the right galaxy directory
    os.chdir(galaxy)

    # Get images
    images = np.load('phot.npy')

    if(imagetype != 'all'):
        imagetypes = [imagetype]
    else:
        imagetypes = i.keys()

    for cimagetype in imagetypes:
        filename = "{g}-{i}.png".format(g=galaxy, i=cimagetype)
        print(filename)
        image(images=images, imagetype=cimagetype, filename=filename)
