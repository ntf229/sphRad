import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import datetime
import matplotlib.patches as mpatches
import fitsio

def directoryStructure(tauClear, dustFraction):
	particlePath = resultPath+'resources/NIHAO/Particles/'
	SKIRTPath = resultPath+'resources/bestParamsRT/'
	noDustSKIRTPath = resultPath+'resources/bestParamsRT/'
	plotPath = resultPath+'resources/bestParamsPlots/'
	if eval(args.ageSmooth):
	    particlePath += 'ageSmooth/'
	    SKIRTPath += 'ageSmooth/'
	    noDustSKIRTPath += 'ageSmooth/'
	    plotPath += 'ageSmooth/'
	else:
	    particlePath += 'noAgeSmooth/'
	    SKIRTPath += 'noAgeSmooth/'
	    noDustSKIRTPath += 'noAgeSmooth/'
	    plotPath += 'noAgeSmooth/'
	if eval(args.SF):
	    particlePath += 'SF/tauClear'+np.format_float_positional(tauClear, trim='-')+'/'
	    SKIRTPath += 'SF/tauClear'+np.format_float_positional(tauClear, trim='-')+'/'
	    noDustSKIRTPath += 'SF/tauClear'+np.format_float_positional(tauClear, trim='-')+'/'
	    plotPath += 'SF/tauClear'+np.format_float_positional(tauClear, trim='-')+'/'
	else:
	    particlePath += 'noSF/'
	    SKIRTPath += 'noSF/'
	    noDustSKIRTPath += 'noSF/'
	    plotPath += 'noSF/'
	noDustPlotPath = plotPath
	SKIRTPath += 'dust/dustFraction'+str(dustFraction)+'/maxTemp'+args.maxTemp+'/'
	plotPath += 'dust/dustFraction'+str(dustFraction)+'/maxTemp'+args.maxTemp+'/'
	noDustSKIRTPath += 'noDust/'
	noDustPlotPath += 'noDust/'
	if eval(args.clumps):
	    particlePath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
	    SKIRTPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
	    plotPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
	else:
	    particlePath += 'noClumps/'
	    SKIRTPath += 'noClumps/'
	    plotPath += 'noClumps/'
	SKIRTPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
	plotPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
	noDustSKIRTPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
	noDustPlotPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
	return SKIRTPath, plotPath, noDustSKIRTPath, noDustPlotPath, particlePath

def makeImages(galaxy):
    nameMask = names == galaxy
    os.system('mkdir -p '+plotPath+'Images/composite/'+galaxy+'/')
    for i in range(len(flux[nameMask][:,0,0,0])): # loop over orientations
        #for j in range(numOrientations):
        # RGB = z, r, g:
        r_grid = flux[nameMask][i,6,:,:] # sdss_z band (see band_names for indicies)
        g_grid = flux[nameMask][i,4,:,:] # sdss_r band
        b_grid = flux[nameMask][i,3,:,:] # sdss_g band
        fuv_grid = flux[nameMask][i,0,:,:]
        w4_grid = flux[nameMask][i,13,:,:]
        # set brightest pixels to value of 50th brightest pixel
        numVisable = len(r_grid[r_grid > np.amax(r_grid)*1e-5])
        numBrightest = int(numVisable*0.001)
        print('numBrightest:',numBrightest)
        max_r = r_grid.flatten()[np.argsort(r_grid.flatten())][-numBrightest]
        max_g = g_grid.flatten()[np.argsort(g_grid.flatten())][-numBrightest]
        max_b = b_grid.flatten()[np.argsort(b_grid.flatten())][-numBrightest]
        max_fuv = fuv_grid.flatten()[np.argsort(fuv_grid.flatten())][-numBrightest]
        max_w4 = w4_grid.flatten()[np.argsort(w4_grid.flatten())][-numBrightest]
        r_grid[r_grid > max_r] = max_r
        g_grid[g_grid > max_g] = max_g
        b_grid[b_grid > max_b] = max_b
        fuv_grid[fuv_grid > max_fuv] = max_fuv
        w4_grid[w4_grid > max_w4] = max_w4
        tot_max = np.amax([np.amax(r_grid),np.amax(g_grid),np.amax(b_grid)])
        stretch_image = 0.005 # increase to make dim pixels dimmer
        stretch_image2 = 0.001
        image_r = np.arcsinh((r_grid/np.amax(r_grid))/stretch_image) 
        image_g = np.arcsinh((g_grid/np.amax(g_grid))/stretch_image) 
        image_b = np.arcsinh((b_grid/np.amax(b_grid))/stretch_image) 
        image_fuv = np.arcsinh((fuv_grid/np.amax(fuv_grid))/stretch_image2) 
        image_w4 = np.arcsinh((w4_grid/np.amax(w4_grid))/stretch_image2) 
        fig = plt.figure()
        sizes = np.shape(r_grid)
        fig.set_size_inches(1. * sizes[0] / sizes[1], 1, forward = False)
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)
        red = (image_r*0.70+image_w4*0.3)/np.amax(image_r*0.70+image_w4*0.3)
        green = image_g/np.amax(image_g)*0.75
        blue = (image_b*0.70+image_fuv*0.3)/np.amax(image_b*0.70+image_fuv*0.3)
        image = np.transpose(np.asarray([red,green,blue]))
        ax.imshow(image, interpolation='none')
        plt.savefig(plotPath+'Images/composite/'+galaxy+'/axisRatio'+str(axisRatios[nameMask][i])+'.png', dpi=sizes[0])
        plt.close()

parser = argparse.ArgumentParser()
parser.add_argument("--ageSmooth") # if True, smooth ages based on number of star particles (makeTextFiles parameter)
parser.add_argument("--SF") # if True, star particles younger than 10 Myrs are assigned MAPPINGS-III SEDs (makeTextFiles parameter)
#parser.add_argument("--tauClear") # clearing time in Myrs for MAPPINGS-III f_PDR calculations (only matters if SF=True) (makeTextFiles parameter)
parser.add_argument("--clumps") # if True, add subgrid clumpiness to gas near MAPPINGS-III particles (makeTextFiles parameter)
parser.add_argument("--numCells") # number of cells along one dimension for subgrid clumping (only matters if clumps=True) (makeTextFiles parameter)
parser.add_argument("--numClumps") # number of clumps for subgrid clumping (only matters if clumps=True) (makeTextFiles parameter)
parser.add_argument("--numPhotons") # number of photon packages (SKIRT parameter)
parser.add_argument("--pixels") # number of pixels (square) for image (SKIRT parameter)
#parser.add_argument("--dustFraction") # dust to metal ratio (SKIRT parameter)
parser.add_argument("--maxTemp") # maximum temperature at which dust can form (SKIRT parameter)
parser.add_argument("--SSP") # simple stellar population model including IMF after underscore (SKIRT parameter)
args = parser.parse_args()

codePath = '/home/ntf229/sphRad/'
resultPath = '/scratch/ntf229/sphRad/' # store results here
fitsPath = '/scratch/ntf229/sphRad/resources/bestParamsFits/ageSmooth/'

massPath = resultPath+'resources/NIHAO/GlobalProps/stellarMasses/'
SFRPath = resultPath+'resources/NIHAO/GlobalProps/SFR/'
selectedPath = resultPath+'resources/selectedOrientations/'
if eval(args.ageSmooth):
    SFRPath += 'ageSmooth/'
else:
    SFRPath += 'noAgeSmooth/'

# Best parameters
tauClear = 2.5
dustFraction = 0.1

SKIRTPath, plotPath, noDustSKIRTPath, noDustPlotPath, particlePath = directoryStructure(tauClear, dustFraction)

# SKIRT broadband photometry names 
band_names = ['FUV', 'NUV', 'u', 'g', 'r', 'i', 'z', '2MASS_J', '2MASS_H', '2MASS_KS', 'W1', 
                'W2', 'W3', 'W4', 'PACS70', 'PACS100', 'PACS160', 'SPIRE250', 'SPIRE350', 'SPIRE500']

# Final sample (65) excluding g3.19e10 (conatins two galaxies)
galaxies = ['g1.88e10','g1.89e10','g1.90e10','g2.34e10','g2.63e10','g2.64e10','g2.80e10','g2.83e10','g2.94e10',
            'g3.44e10','g3.67e10','g3.93e10','g4.27e10','g4.48e10','g4.86e10','g4.94e10','g4.99e10',
            'g5.05e10','g6.12e10','g6.37e10','g6.77e10','g6.91e10','g7.12e10','g8.89e10','g9.59e10','g1.05e11',
            'g1.08e11','g1.37e11','g1.52e11','g1.57e11','g1.59e11','g1.64e11','g2.04e11','g2.19e11','g2.39e11',
            'g2.41e11','g2.42e11','g2.54e11','g3.06e11','g3.21e11','g3.23e11','g3.49e11','g3.55e11','g3.59e11',
            'g3.71e11','g4.90e11','g5.02e11','g5.31e11','g5.36e11','g5.38e11','g5.46e11','g5.55e11','g6.96e11',
            'g7.08e11','g7.44e11','g7.55e11','g7.66e11','g8.06e11','g8.13e11','g8.26e11','g8.28e11','g1.12e12',
            'g1.77e12','g1.92e12','g2.79e12']

summary = fitsio.read(fitsPath+'nihao-resolved-photometry.fits', ext='SUMMARY')

names = summary['name']
stellarMass = summary['stellar_mass']
SFR = summary['sfr']
dustMass = summary['dust_mass']
axisRatios = summary['axis_ratio']
bands = summary['bands'][0]
flux = summary['flux']
#flux_noDust = summary['flux_nodust']

# flux shape: (10, 20, 500, 500)

sSFR = SFR / stellarMass

sphMassMask = np.log10(stellarMass) > 9.5 

# skip low mass galaxies:
#for i in range(len(galaxies)):
#    nameMask = names == galaxies[i]
#    if sphMassMask[nameMask][0]:
#        print('including '+galaxies[i])
#        makeImages(galaxies[i])
#    else:
#        print('skipping '+galaxies[i])

for i in range(len(galaxies)):
    makeImages(galaxies[i])

print('done')
