import numpy as np
import os
import argparse
from astropy.io import fits
from os.path import expanduser

parser = argparse.ArgumentParser()
parser.add_argument("--ageSmooth") # if True, smooth ages based on number of star particles (makeTextFiles parameter)
parser.add_argument("--SF") # if True, star particles younger than 10 Myrs are assigned MAPPINGS-III SEDs (makeTextFiles parameter)
parser.add_argument("--tauClear") # clearing time in Myrs for MAPPINGS-III f_PDR calculations (only matters if SF=True) (makeTextFiles $
parser.add_argument("--clumps") # if True, add subgrid clumpiness to gas near MAPPINGS-III particles (makeTextFiles parameter)
parser.add_argument("--numCells") # number of cells along one dimension for subgrid clumping (only matters if clumps=True) (makeTextFil$
parser.add_argument("--numClumps") # number of clumps for subgrid clumping (only matters if clumps=True) (makeTextFiles parameter)
parser.add_argument("--dust") # include dust; True or False
parser.add_argument("--inc") # inclination angle (SKIRT parameter)
parser.add_argument("--az") # azimuth angle (SKIRT parameter)
parser.add_argument("--numPhotons") # number of photon packages (SKIRT parameter)
parser.add_argument("--pixels") # number of pixels (square) for image (SKIRT parameter)
parser.add_argument("--dustFraction") # dust to metal ratio (SKIRT parameter)
parser.add_argument("--maxTemp") # maximum temperature at which dust can form (SKIRT parameter)
parser.add_argument("--SSP") # simple stellar population model including IMF after underscore (SKIRT parameter)
args = parser.parse_args()

origDir = os.getcwd()
codePath=expanduser('~')+'/sphRad/'
resultPath = '/scratch/ntf229/sphRad/' # store results here

# Directory structure stores important parameters
photPath = resultPath+'resources/Photometry/'
SKIRTPath = resultPath+'resources/SKIRT/'
if eval(args.ageSmooth):
    photPath += 'ageSmooth/'
    SKIRTPath += 'ageSmooth/'
else:
    photPath += 'noAgeSmooth/'
    SKIRTPath += 'noAgeSmooth/'
if eval(args.SF):
    photPath += 'SF/tauClear'+args.tauClear+'/'
    SKIRTPath += 'SF/tauClear'+args.tauClear+'/'
else:
    photPath += 'noSF/'
    SKIRTPath += 'noSF/'
if eval(args.clumps):
    photPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
    SKIRTPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
else:
    photPath += 'noClumps/'
    SKIRTPath += 'noClumps/'
if eval(args.dust):
    photPath += 'dust/dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'
    SKIRTPath += 'dust/dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'
else:
    photPath += 'noDust/'
    SKIRTPath += 'noDust/'

photPath += 'numPhotons'+args.numPhotons+'/inc'+args.inc+'/az'+args.az+'/'+args.SSP+'/'
SKIRTPath += 'numPhotons'+args.numPhotons+'/inc'+args.inc+'/az'+args.az+'/'+args.SSP+'/'

band_names = ['FUV', 'NUV', 'u', 'g', 'r', 'i', 'z', '2MASS_J', '2MASS_H', '2MASS_KS', 'W1', 'W2', 'W3', 'W4', 'PACS70', 'PACS100', 'PACS160', 'SPIRE250', 'SPIRE350', 'SPIRE500']

# Final sample (66)
galaxies = ['g1.88e10','g1.89e10','g1.90e10','g2.34e10','g2.63e10','g2.64e10','g2.80e10','g2.83e10','g2.94e10',
            'g3.19e10','g3.44e10','g3.67e10','g3.93e10','g4.27e10','g4.48e10','g4.86e10','g4.94e10','g4.99e10',
            'g5.05e10','g6.12e10','g6.37e10','g6.77e10','g6.91e10','g7.12e10','g8.89e10','g9.59e10','g1.05e11',
            'g1.08e11','g1.37e11','g1.52e11','g1.57e11','g1.59e11','g1.64e11','g2.04e11','g2.19e11','g2.39e11',
            'g2.41e11','g2.42e11','g2.54e11','g3.06e11','g3.21e11','g3.23e11','g3.49e11','g3.55e11','g3.59e11',
            'g3.71e11','g4.90e11','g5.02e11','g5.31e11','g5.36e11','g5.38e11','g5.46e11','g5.55e11','g6.96e11',
            'g7.08e11','g7.44e11','g7.55e11','g7.66e11','g8.06e11','g8.13e11','g8.26e11','g8.28e11','g1.12e12',
            'g1.77e12','g1.92e12','g2.79e12']

# just g7.55e11 edge on 
#galaxies = ['g7.55e11']

for i in range(len(galaxies)):
	savePath = photPath+galaxies[i]+'/'
	os.system('mkdir -p '+savePath)
	file = fits.open(SKIRTPath+galaxies[i]+'/sph_broadband_total.fits')
	data = np.asarray(file[0].data) # (2000, 2000, 20) bands in MJy/st
	np.save(savePath+'phot.npy', data)
	sed = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_SED_sed.dat', unpack = True)
	wave = sed[0] # spatially integrated SED wavelengths in microns
	wave = wave * 1e4 # convert to Angstroms
	spec = sed[1] # spatially integrated SED fluxes in Janskys
	np.save(savePath+'spec.npy', spec)
	np.save(savePath+'wave.npy', wave)

print('done')
