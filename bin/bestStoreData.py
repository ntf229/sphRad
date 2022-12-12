import numpy as np
import os
from astropy.io import fits
import argparse
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from timeit import default_timer as timer
import datetime
import matplotlib as mpl
import matplotlib.patches as mpatches

def directoryStructure(tauClear, dustFraction):
	particlePath = resultPath+'resources/NIHAO/Particles/'
	SKIRTPath = resultPath+'resources/selectedOrientations_SKIRT/'
	noDustSKIRTPath = resultPath+'resources/selectedOrientations_SKIRT/'
	if eval(args.ageSmooth):
	    particlePath += 'ageSmooth/'
	    SKIRTPath += 'ageSmooth/'
	    noDustSKIRTPath += 'ageSmooth/'
	else:
	    particlePath += 'noAgeSmooth/'
	    SKIRTPath += 'noAgeSmooth/'
	    noDustSKIRTPath += 'noAgeSmooth/'
	if eval(args.SF):
	    particlePath += 'SF/tauClear'+np.format_float_positional(tauClear, trim='-')+'/'
	    SKIRTPath += 'SF/tauClear'+np.format_float_positional(tauClear, trim='-')+'/'
	    noDustSKIRTPath += 'SF/tauClear'+np.format_float_positional(tauClear, trim='-')+'/'
	else:
	    particlePath += 'noSF/'
	    SKIRTPath += 'noSF/'
	    noDustSKIRTPath += 'noSF/'
	SKIRTPath += 'dust/dustFraction'+str(dustFraction)+'/maxTemp'+args.maxTemp+'/'
	noDustSKIRTPath += 'noDust/'
	if eval(args.clumps):
	    particlePath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
	    SKIRTPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
	else:
	    particlePath += 'noClumps/'
	    SKIRTPath += 'noClumps/'
	SKIRTPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
	noDustSKIRTPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
	return SKIRTPath, noDustSKIRTPath, particlePath

def stellarMass(galaxy):
    if os.path.isfile(massPath+galaxy+'/stellarMass.npy'):
        SKIRT_stellarMass = float(np.load(massPath+galaxy+'/stellarMass.npy'))
    else:
        stars = np.load(particlePath+galaxy+'/stars.npy')
        youngStars = np.load(particlePath+galaxy+'/youngStars.npy')
        SKIRT_stellarMass = np.sum(stars[:,7]) + (np.sum(youngStars[:,7] * 1.e7)) # in Msun
        os.system('mkdir -p '+massPath+galaxy+'/')
        np.save(massPath+galaxy+'/stellarMass.npy', SKIRT_stellarMass)
    return SKIRT_stellarMass

def SFR(galaxy):
    if os.path.isfile(SFRPath+galaxy+'/SFR.npy'):
        SKIRT_SFR = float(np.load(SFRPath+galaxy+'/SFR.npy'))
    else:
        stars = np.load(particlePath+galaxy+'/stars.npy')
        youngStars = np.load(particlePath+galaxy+'/youngStars.npy')
        youngMask = stars[:,9] < 1.e8 # younger than 100 Myrs
        SKIRT_SFR = (np.sum(stars[youngMask,7]) / 1.e8) + np.sum(youngStars[:,7]) # in Msun per year
        os.system('mkdir -p '+SFRPath+galaxy+'/')
        np.save(SFRPath+galaxy+'/SFR.npy', SKIRT_SFR)
    return SKIRT_SFR

def dustMass(galaxy):
    if os.path.isfile(massPath+galaxy+'/maxTemp'+args.maxTemp+'/metalMass.npy'):
        metalMass = float(np.load(massPath+galaxy+'/maxTemp'+args.maxTemp+'/metalMass.npy'))
    else:
        gas = np.load(particlePath+galaxy+'/gas.npy')
        tempMask = gas[:,6] < float(args.maxTemp)
        ghostMask = np.asarray(gas[tempMask][:,4] > 0, dtype=bool) # mask out negative mass ghost particles
        metalMass = np.sum(gas[tempMask, 4][ghostMask] * gas[tempMask, 5][ghostMask]) # in Msun
        os.system('mkdir -p '+massPath+galaxy+'/maxTemp'+args.maxTemp+'/')
        np.save(massPath+galaxy+'/maxTemp'+args.maxTemp+'/metalMass.npy', metalMass)
    SKIRT_dustMass = metalMass * dustFraction # len(dustFractions)
    return SKIRT_dustMass

def getAvValues():
    AvValues = np.zeros((len(galaxies), numOrientations))
    for i in range(len(galaxies)):
        for j in range(numOrientations):
            instName = 'axisRatio'+str(np.round_(SKIRT_axisRatio[i,j], decimals = 4))
            sed = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_SED_'+instName+'_sed.dat', unpack = True)
            wave = sed[0] * 1e4 # spatially integrated SED wavelengths converted to Angstroms
            spec = sed[1] # spatially integrated SED fluxes in Janskys
            sed = np.loadtxt(noDustSKIRTPath+galaxies[i]+'/sph_SED_'+instName+'_sed.dat', unpack = True)
            noDustSpec = sed[1] 
            att_mask = (wave >= 912) & (wave <= 2e4)
            dustMags = 22.5 - 2.5*np.log10((spec/3631) * 1e9) # convert to Pogson magnitudes
            noDustMags = 22.5 - 2.5*np.log10((noDustSpec/3631) * 1e9) # convert to Pogson magnitudes
            attenuation = dustMags[att_mask] - noDustMags[att_mask]
            Av_index = np.abs(wave[att_mask] - 5500).argmin() # find wave index closest to 5500 angstroms (V)
            AvValues[i,j] = attenuation[Av_index]
    return AvValues

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
storePath = resultPath+'resources/bestData/'
massPath = resultPath+'resources/NIHAO/GlobalProps/stellarMasses/'
SFRPath = resultPath+'resources/NIHAO/GlobalProps/SFR/'
selectedPath = resultPath+'resources/selectedOrientations/'

if eval(args.ageSmooth):
    SFRPath += 'ageSmooth/'
else:
    SFRPath += 'noAgeSmooth/'

# Best parameters
tauClear = 3.
dustFraction = 0.12

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

# Initialize arrays for SKIRT spatially integrated photometry and axis ratios for color plots
numOrientations = 10
#SKIRT_flux = np.zeros((len(galaxies), numOrientations, len(band_names)))
#noDustSKIRT_flux = np.zeros((len(galaxies), numOrientations, len(band_names)))
SKIRT_axisRatio = np.zeros((len(galaxies), numOrientations))
SKIRT_stellarMass = np.zeros(len(galaxies))
SKIRT_SFR = np.zeros(len(galaxies))
SKIRT_dustMass = np.zeros(len(galaxies))
SKIRT_AvValues = np.zeros((len(galaxies), numOrientations))

SKIRTPath, noDustSKIRTPath, particlePath = directoryStructure(tauClear, dustFraction)

# Fill up SKIRT arrays 
for i in range(len(galaxies)):
    os.system('mkdir -p '+storePath+galaxies[i]+'/')
    print('starting '+galaxies[i])
    SKIRT_stellarMass[i] = stellarMass(galaxies[i])
    SKIRT_SFR[i] = SFR(galaxies[i])
    SKIRT_dustMass[i] = dustMass(galaxies[i])
    # import axis ratios from table
    selections = np.load(selectedPath+galaxies[i]+'/selectedAxisRatios.npy') # [inc, az, axisRatio]
    SKIRT_axisRatio[i,:] = selections[:,2]
    #np.save(storePath+galaxies[i]+'/axisRatios.npy', str(np.round_(SKIRT_axisRatio[i,:], decimals = 4)))
    np.save(storePath+galaxies[i]+'/axisRatios.npy', SKIRT_axisRatio[i,:])
    np.save(storePath+galaxies[i]+'/stellarMass.npy', SKIRT_stellarMass[i])
    np.save(storePath+galaxies[i]+'/SFR.npy', SKIRT_SFR[i])
    np.save(storePath+galaxies[i]+'/dustMass.npy', SKIRT_dustMass[i])
    #sphFaceOnMask[i,:] = SKIRT_axisRatio[i,:] == np.amax(SKIRT_axisRatio[i,:]) # most face-on orientation
    for j in range(numOrientations):
        axisRatioStr = str(np.round_(SKIRT_axisRatio[i,j], decimals = 4))
        instName = 'axisRatio'+axisRatioStr
        bb = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_broadband_'+instName+'_sed.dat', unpack = True)
        #SKIRT_flux[i,j,:] = bb[1] # spatially integrated broadband fluxes in Janskys
        noDustbb = np.loadtxt(noDustSKIRTPath+galaxies[i]+'/sph_broadband_'+instName+'_sed.dat', unpack = True)
        #noDustSKIRT_flux[i,j,:] = noDustbb[1] # spatially integrated broadband fluxes in Janskys
        sed = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_SED_'+instName+'_sed.dat', unpack = True)
        wave = sed[0] * 1e4 # spatially integrated SED wavelengths converted to Angstroms
        spec = sed[1] # spatially integrated SED fluxes in Janskys
        sed = np.loadtxt(noDustSKIRTPath+galaxies[i]+'/sph_SED_'+instName+'_sed.dat', unpack = True)
        noDustSpec = sed[1] 
        att_mask = (wave >= 912) & (wave <= 2e4)
        dustMags = 22.5 - 2.5*np.log10((spec/3631) * 1e9) # convert to Pogson magnitudes
        noDustMags = 22.5 - 2.5*np.log10((noDustSpec/3631) * 1e9) # convert to Pogson magnitudes
        attenuation = dustMags[att_mask] - noDustMags[att_mask]
        Av_index = np.abs(wave[att_mask] - 5500).argmin() # find wave index closest to 5500 angstroms (V)
        Av = attenuation[Av_index]
        np.save(storePath+galaxies[i]+'/BB_'+instName+'.npy', bb[1]) # spatially integrated broadband fluxes in Janskys
        np.save(storePath+galaxies[i]+'/noDustBB_'+instName+'.npy', noDustbb[1])
        np.save(storePath+galaxies[i]+'/Av_'+instName+'.npy', Av)
        np.save(storePath+galaxies[i]+'/specWavelengths_'+instName+'.npy', wave) # in Angstroms
        np.save(storePath+galaxies[i]+'/spec_'+instName+'.npy', spec) # in Janskys
        np.save(storePath+galaxies[i]+'/noDustSpec_'+instName+'.npy', noDustSpec)

#SKIRT_sSFR = SKIRT_SFR / SKIRT_stellarMass
#SKIRT_dustToStellar = SKIRT_dustMass / SKIRT_stellarMass

#SKIRT_AvValues = getAvValues() 	

print('done')
