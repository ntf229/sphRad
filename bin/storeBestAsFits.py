import numpy as np
import os
import argparse
from astropy.io import fits
from os.path import expanduser
import fitsio

def getStellarMass(galaxy):
    if os.path.isfile(massPath+galaxy+'/stellarMass.npy'):
        SKIRT_stellarMass = float(np.load(massPath+galaxy+'/stellarMass.npy'))
    else:
        stars = np.load(particlePath+galaxy+'/stars.npy')
        youngStars = np.load(particlePath+galaxy+'/youngStars.npy')
        SKIRT_stellarMass = np.sum(stars[:,7]) + (np.sum(youngStars[:,7] * 1.e7)) # in Msun
        os.system('mkdir -p '+massPath+galaxy+'/')
        np.save(massPath+galaxy+'/stellarMass.npy', SKIRT_stellarMass)
    return SKIRT_stellarMass

def getSFR(galaxy):
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

def getDustMass(galaxy):
    if os.path.isfile(massPath+galaxy+'/maxTemp'+args.maxTemp+'/metalMass.npy'):
        metalMass = float(np.load(massPath+galaxy+'/maxTemp'+args.maxTemp+'/metalMass.npy'))
    else:
        gas = np.load(particlePath+galaxy+'/gas.npy')
        tempMask = gas[:,6] < float(args.maxTemp)
        ghostMask = np.asarray(gas[tempMask][:,4] > 0, dtype=bool) # mask out negative mass ghost particles
        metalMass = np.sum(gas[tempMask, 4][ghostMask] * gas[tempMask, 5][ghostMask]) # in Msun
        os.system('mkdir -p '+massPath+galaxy+'/maxTemp'+args.maxTemp+'/')
        np.save(massPath+galaxy+'/maxTemp'+args.maxTemp+'/metalMass.npy', metalMass)
    SKIRT_dustMass = metalMass * dustFraction
    return SKIRT_dustMass

def getAv(galaxy, instName):
    att_mask = (wave >= 912) & (wave <= 2e4)
    dustMags = 22.5 - 2.5*np.log10((spec/3631) * 1e9) # convert to Pogson magnitudes
    noDustMags = 22.5 - 2.5*np.log10((noDustSpec/3631) * 1e9) # convert to Pogson magnitudes
    attenuation = dustMags[att_mask] - noDustMags[att_mask]
    Av_index = np.abs(wave[att_mask] - 5500).argmin() # find wave index closest to 5500 angstroms (V)
    Av = attenuation[Av_index]
    return Av

def reduceImageSize(image): # (x,y,band)
    newImage = np.zeros((int(len(image[:,0,0])/4), int(len(image[0,:,0])/4), len(image[0,0,:])))
    for i in range(len(newImage[0,0,:])):
        for j in range(len(newImage[:,0,0])):
            for k in range(len(newImage[0,:,0])):
                newImage[j,k,i] = np.mean(image[j*4:j*4+4,k*4:k*4+4,i])

parser = argparse.ArgumentParser()
parser.add_argument("--ageSmooth") # if True, smooth ages based on number of star particles (makeTextFiles parameter)
parser.add_argument("--SF") # if True, star particles younger than 10 Myrs are assigned MAPPINGS-III SEDs (makeTextFiles parameter)
parser.add_argument("--tauClear") # clearing time in Myrs for MAPPINGS-III f_PDR calculations (only matters if SF=True) (makeTextFiles $
parser.add_argument("--clumps") # if True, add subgrid clumpiness to gas near MAPPINGS-III particles (makeTextFiles parameter)
parser.add_argument("--numCells") # number of cells along one dimension for subgrid clumping (only matters if clumps=True) (makeTextFil$
parser.add_argument("--numClumps") # number of clumps for subgrid clumping (only matters if clumps=True) (makeTextFiles parameter)
parser.add_argument("--numPhotons") # number of photon packages (SKIRT parameter)
parser.add_argument("--pixels") # number of pixels (square) for image (SKIRT parameter)
parser.add_argument("--dustFraction") # dust to metal ratio (SKIRT parameter)
parser.add_argument("--maxTemp") # maximum temperature at which dust can form (SKIRT parameter)
parser.add_argument("--SSP") # simple stellar population model including IMF after underscore (SKIRT parameter)
args = parser.parse_args()

storeImages = False

numPixels = int(args.pixels)
reducedPixels = int(numPixels/4)

dustFraction = float(args.dustFraction)

origDir = os.getcwd()
codePath=expanduser('~')+'/sphRad/'
resultPath = '/scratch/ntf229/sphRad/' # store results here
selectedPath = resultPath+'resources/selectedOrientations/'
massPath = resultPath+'resources/NIHAO/GlobalProps/stellarMasses/'
SFRPath = resultPath+'resources/NIHAO/GlobalProps/SFR/'
savePath = resultPath+'resources/bestParamsFits/'
os.system('mkdir -p '+savePath)

# Directory structure stores important parameters
SKIRTPath = resultPath+'resources/bestParamsRT/'
particlePath = resultPath+'resources/NIHAO/Particles/'

if eval(args.ageSmooth):
    SKIRTPath += 'ageSmooth/'
    particlePath += 'ageSmooth/'
else:
    SKIRTPath += 'noAgeSmooth/'
    particlePath += 'noAgeSmooth/'
if eval(args.SF):
    SKIRTPath += 'SF/tauClear'+args.tauClear+'/'
    particlePath += 'SF/tauClear'+args.tauClear+'/'
else:
    SKIRTPath += 'noSF/'
    particlePath += 'noSF/'

noDustSKIRTPath = SKIRTPath + 'noDust/'
SKIRTPath += 'dust/dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'

if eval(args.clumps):
    SKIRTPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
    particlePath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
else:
    SKIRTPath += 'noClumps/'
    particlePath += 'noClumps/'

SKIRTPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
noDustSKIRTPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'

band_names = np.asarray(['FUV', 'NUV', 'u', 'g', 'r', 'i', 'z', '2MASS_J', 
              '2MASS_H', '2MASS_KS', 'W1', 'W2', 'W3', 'W4', 'PACS70', 
              'PACS100', 'PACS160', 'SPIRE250', 'SPIRE350', 'SPIRE500'])

# Final sample (65) excluding g3.19e10 (two galaxies)
names = ['g1.88e10','g1.89e10','g1.90e10','g2.34e10','g2.63e10','g2.64e10','g2.80e10','g2.83e10','g2.94e10',
            'g3.44e10','g3.67e10','g3.93e10','g4.27e10','g4.48e10','g4.86e10','g4.94e10','g4.99e10',
            'g5.05e10','g6.12e10','g6.37e10','g6.77e10','g6.91e10','g7.12e10','g8.89e10','g9.59e10','g1.05e11',
            'g1.08e11','g1.37e11','g1.52e11','g1.57e11','g1.59e11','g1.64e11','g2.04e11','g2.19e11','g2.39e11',
            'g2.41e11','g2.42e11','g2.54e11','g3.06e11','g3.21e11','g3.23e11','g3.49e11','g3.55e11','g3.59e11',
            'g3.71e11','g4.90e11','g5.02e11','g5.31e11','g5.36e11','g5.38e11','g5.46e11','g5.55e11','g6.96e11',
            'g7.08e11','g7.44e11','g7.55e11','g7.66e11','g8.06e11','g8.13e11','g8.26e11','g8.28e11','g1.12e12',
            'g1.77e12','g1.92e12','g2.79e12']

numOrientations = 10

galaxies_dtype = [('name', str, 20),
                  ('stellar_mass', np.float32),
                  ('sfr', np.float32),
                  ('dust_mass', np.float32),
                  ('axis_ratios', np.float32, numOrientations)]

summary_dtype = [('name', str, 20),
                 ('stellar_mass', np.float32),
                 ('sfr', np.float32),
                 ('dust_mass', np.float32),
                 ('axis_ratio', np.float32),
                 ('Av', np.float32),
                 ('bands', band_names.dtype, len(band_names)),
                 ('flux', np.float32, len(band_names)),
                 ('flux_nodust', np.float32, len(band_names))]

image_summary_dtype = [('name', str, 20),
                 ('stellar_mass', np.float32),
                 ('sfr', np.float32),
                 ('dust_mass', np.float32),
                 ('axis_ratio', np.float32),
                 ('Av', np.float32),
                 ('bands', band_names.dtype, len(band_names)),
                 ('flux', np.float32, reducedPixels*reducedPixels*len(band_names)),
                 ('flux_nodust', np.float32, reducedPixels*reducedPixels*len(band_names))]

galaxies = np.zeros(len(names), dtype=galaxies_dtype)
summary = np.zeros(len(names) * numOrientations, dtype=summary_dtype)
if storeImages:
    image_summary = np.zeros(len(names) * numOrientations, dtype=image_summary_dtype)

wave = None

indx = 0

for i in range(len(names)):
    stellarMass = getStellarMass(names[i])
    SFR = getSFR(names[i])
    dustMass = getDustMass(names[i])
    selections = np.load(selectedPath+names[i]+'/selectedIncAzAR.npy') # [inc, az, axisRatio]
    axisRatios = selections[:,2]
    galaxies['name'][i] = names[i]
    galaxies['stellar_mass'][i] = stellarMass
    galaxies['sfr'][i] = SFR
    galaxies['dust_mass'][i] = dustMass
    galaxies['axis_ratios'][i] = axisRatios
    for j in range(numOrientations): # loop through orientations
        instName = 'axisRatio' + str(np.round_(axisRatios[j], decimals=4)) # orientation naming system
        # Spatially integrated
        BB = np.loadtxt(SKIRTPath+names[i]+'/sph_broadband_'+instName+'_sed.dat', unpack = True)[1] # spatially integrated broadband fluxes in Janskys
        noDustBB = np.loadtxt(noDustSKIRTPath+names[i]+'/sph_broadband_'+instName+'_sed.dat', unpack = True)[1]  # spatially integrated broadband fluxes in Janskys
        sed = np.loadtxt(SKIRTPath+names[i]+'/sph_SED_'+instName+'_sed.dat', unpack = True)
        if (wave is None): # only need to do once
            wave = sed[0] * 1e4 # spatially integrated SED wavelengths converted to Angstroms
            spectrum = np.zeros((len(names) * numOrientations, len(wave)), dtype=np.float32)
            spectrum_nodust = np.zeros((len(names) * numOrientations, len(wave)), dtype=np.float32)
        spec = sed[1] # spatially integrated SED fluxes in Janskys
        noDustSed = np.loadtxt(noDustSKIRTPath+names[i]+'/sph_SED_'+instName+'_sed.dat', unpack = True)
        noDustSpec = noDustSed[1]
        Av = getAv(names[i], instName)
        # Store integrated data
        summary['name'][indx] = galaxies['name'][i]
        summary['stellar_mass'][indx] = galaxies['stellar_mass'][i]
        summary['sfr'][indx] = galaxies['sfr'][i]
        summary['dust_mass'][indx] = galaxies['dust_mass'][i]
        summary['axis_ratio'][indx] = galaxies['axis_ratios'][i, j]
        summary['Av'][indx] = Av
        summary['bands'][indx] = band_names
        summary['flux'][indx, :] = BB
        summary['flux_nodust'][indx, :] = noDustBB
        spectrum[indx, :] = spec
        spectrum_nodust[indx, :] = noDustSpec
        # Spatially resolved
        if storeImages:
            imageFile = fits.open(SKIRTPath+names[i]+'/sph_broadband_'+instName+'_total.fits')
            imageBB = np.asarray(imageFile[0].data) # (2000, 2000, 20) bands in MJy/st
            imageBB = reduceImageSize(imageBB) # (500, 500, 20)
            noDustImageFile = fits.open(noDustSKIRTPath+names[i]+'/sph_broadband_'+instName+'_total.fits')
            noDustImageBB = np.asarray(noDustImageFile[0].data) # (2000, 2000, 20) bands in MJy/st
            noDustImageBB = reduceImageSize(noDustImageBB) # (500, 500, 20)
            # Store resolved data
            image_summary['name'][indx] = galaxies['name'][i]
            image_summary['stellar_mass'][indx] = galaxies['stellar_mass'][i]
            image_summary['sfr'][indx] = galaxies['sfr'][i]
            image_summary['dust_mass'][indx] = galaxies['dust_mass'][i]
            image_summary['axis_ratio'][indx] = galaxies['axis_ratios'][i, j]
            image_summary['Av'][indx] = Av
            image_summary['bands'][indx] = band_names
            image_summary['flux'][indx, :] = imageBB
            image_summary['flux_nodust'][indx, :] = noDustImageBB
        # Increase indx
        indx += 1

fitsio.write(savePath+'nihao-integrated-seds.fits', galaxies, extname='GALAXIES', clobber=True)
fitsio.write(savePath+'nihao-integrated-seds.fits', summary, extname='SUMMARY', clobber=False)
fitsio.write(savePath+'nihao-integrated-seds.fits', wave, extname='WAVE', clobber=False)
fitsio.write(savePath+'nihao-integrated-seds.fits', spectrum, extname='SPEC', clobber=False)
fitsio.write(savePath+'nihao-integrated-seds.fits', spectrum_nodust, extname='SPECNODUST', clobber=False)

if storeImages:
    fitsio.write(savePath+'nihao-resolved-photometry.fits', galaxies, extname='GALAXIES', clobber=True)
    fitsio.write(savePath+'nihao-resolved-photometry.fits', image_summary, extname='SUMMARY', clobber=False)

print('done')
