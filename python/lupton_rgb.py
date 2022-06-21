import numpy as np
import pts.simulation as sm
from sedpy import observate
import astropy.io.fits as fits
import os
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb

plotPath = '/scratch/ntf229/RT_fit/resources/SKIRT/lupton_images/'
os.system('mkdir -p '+plotPath)

#bands = ["sdss_z0","sdss_r0","sdss_g0"] # corresponding to RGB channels
bands = ["herschel_pacs_100","sdss_g0","galex_FUV"]

global_scaling = False # if True, colors are scaled such that the average color of the population of galaxies is white
#                        if False, colors are scaled such that the average color of each galaxy is white

# All SecondBatch galaxies (69)
#galaxies = ['g1.12e12','g1.92e12','g2.39e11','g2.79e12','g3.23e11','g3.49e11','g3.59e11','g3.61e11','g5.31e11',
#            'g5.36e11','g5.38e11','g5.55e11','g7.08e11','g7.44e11','g8.26e11','g8.28e11','g1.05e11','g1.23e10',
#            'g1.50e10','g1.88e10','g1.95e10','g2.37e10','g2.57e11','g2.83e10','g3.21e11','g4.36e09','g4.99e09',
#            'g5.84e09','g7.12e10','g9.59e10','g1.08e11','g1.37e11','g1.52e11','g1.89e10','g2.04e11','g2.39e10',
#            'g2.63e10','g2.94e10','g3.44e10','g4.48e10','g5.22e09','g6.12e10','g7.34e09','g1.09e10','g1.44e10',
#            'g1.57e10','g1.90e10','g2.09e10','g2.41e11','g2.64e10','g3.06e11','g3.67e10','g4.86e10','g5.41e09',
#            'g6.96e10','g8.89e10','g1.18e10','g1.47e10','g1.64e11','g1.92e10','g2.34e10','g2.54e11','g2.80e10',
#            'g3.19e10','g3.93e10','g4.94e10','g5.59e09','g7.05e09','g9.26e09']

# Final sample (66)
galaxies = ['g1.88e10','g1.89e10','g1.90e10','g2.34e10','g2.63e10','g2.64e10','g2.80e10','g2.83e10','g2.94e10',
            'g3.19e10','g3.44e10','g3.67e10','g3.93e10','g4.27e10','g4.48e10','g4.86e10','g4.94e10','g4.99e10',
            'g5.05e10','g6.12e10','g6.37e10','g6.77e10','g6.91e10','g7.12e10','g8.89e10','g9.59e10','g1.05e11',
            'g1.08e11','g1.37e11','g1.52e11','g1.57e11','g1.59e11','g1.64e11','g2.04e11','g2.19e11','g2.39e11',
            'g2.41e11','g2.42e11','g2.54e11','g3.06e11','g3.21e11','g3.23e11','g3.49e11','g3.55e11','g3.59e11',
            'g3.71e11','g4.90e11','g5.02e11','g5.31e11','g5.36e11','g5.38e11','g5.46e11','g5.55e11','g6.96e11',
            'g7.08e11','g7.44e11','g7.55e11','g7.66e11','g8.06e11','g8.13e11','g8.26e11','g8.28e11','g1.12e12',
            'g1.77e12','g1.92e12','g2.79e12']


# We want each RGB channel to have the same value when averaged over the entire population of galaxies

avg_r = np.zeros(len(galaxies)) # scaled between 0 and 1 
avg_g = np.zeros(len(galaxies))
avg_b = np.zeros(len(galaxies))

if global_scaling:
    for i in range(len(galaxies)):
        filePath = '/scratch/ntf229/RT_fit/resources/SKIRT/'+galaxies[i]+'/maxLevel13/wavelengths601/numPhotons1e9/inc0/dust/dustFraction0.2/maxTemp16000/'
        r_grid = np.load(filePath+'grids/'+bands[0]+'.npy')
        g_grid = np.load(filePath+'grids/'+bands[1]+'.npy')
        b_grid = np.load(filePath+'grids/'+bands[2]+'.npy')
        tot_max = np.amax([np.amax(r_grid),np.amax(g_grid),np.amax(b_grid)]) 
        avg_r[i] = np.average(r_grid / tot_max)
        avg_g[i] = np.average(g_grid / tot_max)
        avg_b[i] = np.average(b_grid / tot_max)
    global_avg_r = np.average(avg_r)
    global_avg_g = np.average(avg_g)
    global_avg_b = np.average(avg_b)

for i in range(len(galaxies)):
    filePath = '/scratch/ntf229/RT_fit/resources/SKIRT/'+galaxies[i]+'/maxLevel13/wavelengths601/numPhotons1e9/inc0/dust/dustFraction0.2/maxTemp16000/'
    r_grid = np.load(filePath+'grids/'+bands[0]+'.npy')
    g_grid = np.load(filePath+'grids/'+bands[1]+'.npy')
    b_grid = np.load(filePath+'grids/'+bands[2]+'.npy') 
    if global_scaling:
        r_grid = r_grid / global_avg_r
        g_grid = g_grid / global_avg_g
        b_grid = b_grid / global_avg_b
        tot_max = np.amax([np.amax(r_grid),np.amax(g_grid),np.amax(b_grid)])
        r_grid = 255 * r_grid / tot_max
        g_grid = 255 * g_grid / tot_max
        b_grid = 255 * b_grid / tot_max
    else:
        r_grid = 255 * r_grid / np.amax(r_grid)
        g_grid = 255 * g_grid / np.amax(g_grid)
        b_grid = 255 * b_grid / np.amax(b_grid)
    plt.figure(figsize=(2000/300, 2000/300), dpi=300)
    image = make_lupton_rgb(r_grid, g_grid, b_grid, stretch=0.2, Q=15)
    plt.axis('off')
    plt.imshow(image,interpolation='none')
    plt.savefig(plotPath+galaxies[i]+'_'+bands[0]+'_'+bands[1]+'_'+bands[2]+'.png',dpi=300)
    plt.close()




