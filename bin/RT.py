# run SKIRT RT from text files (need to run makeTextFiles.py with same arguments first)

import argparse
import os
import shutil
from os.path import expanduser
from timeit import default_timer as timer
import subprocess
import datetime
import numpy as np
import sys

print('The code has started')

parser = argparse.ArgumentParser()
parser.add_argument("--ageSmooth") # if True, smooth ages based on number of star particles (makeTextFiles parameter)
parser.add_argument("--SF") # if True, star particles younger than 10 Myrs are assigned MAPPINGS-III SEDs (makeTextFiles parameter)
parser.add_argument("--tauClear") # clearing time in Myrs for MAPPINGS-III f_PDR calculations (only matters if SF=True) (makeTextFiles parameter)
parser.add_argument("--clumps") # if True, add subgrid clumpiness to gas near MAPPINGS-III particles (makeTextFiles parameter)
parser.add_argument("--numCells") # number of cells along one dimension for subgrid clumping (only matters if clumps=True) (makeTextFiles parameter)
parser.add_argument("--numClumps") # number of clumps for subgrid clumping (only matters if clumps=True) (makeTextFiles parameter)
parser.add_argument("--inc") # inclination angle (SKIRT parameter)
parser.add_argument("--az") # azimuth angle (SKIRT parameter)
parser.add_argument("--numPhotons") # number of photon packages (SKIRT parameter)
parser.add_argument("--pixels") # number of pixels (square) for image (SKIRT parameter)
parser.add_argument("--dustFraction") # dust to metal ratio (SKIRT parameter)
parser.add_argument("--maxTemp") # maximum temperature at which dust can form (SKIRT parameter)
parser.add_argument("--SSP") # simple stellar population model including IMF after underscore (SKIRT parameter)
parser.add_argument("--galaxy") # name of galaxy 
args = parser.parse_args()

origDir = os.getcwd()
codePath=expanduser('~')+'/sphRad/'
resultPath = '/scratch/ntf229/sphRad/' # store results here

# Directory structure stores important parameters
textPath = resultPath+'resources/NIHAO/TextFiles/'
SKIRTPath = resultPath+'resources/SKIRT/'
noDustSKIRTPath = resultPath+'resources/SKIRT/'
if eval(args.ageSmooth):
    textPath += 'ageSmooth/'
    SKIRTPath += 'ageSmooth/'
    noDustSKIRTPath += 'ageSmooth/'
else:
    textPath += 'noAgeSmooth/'
    SKIRTPath += 'noAgeSmooth/'
    noDustSKIRTPath += 'noAgeSmooth/'
if eval(args.SF):
    textPath += 'SF/tauClear'+args.tauClear+'/'
    SKIRTPath += 'SF/tauClear'+args.tauClear+'/'
    noDustSKIRTPath += 'SF/tauClear'+args.tauClear+'/'
else:
	textPath += 'noSF/'
	SKIRTPath += 'noSF/'
	noDustSKIRTPath += 'noSF/'
if eval(args.clumps):
    textPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
    SKIRTPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
    noDustSKIRTPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
else:
	textPath += 'noClumps/'
	SKIRTPath += 'noClumps/'
	noDustSKIRTPath += 'noClumps/'

SKIRTPath += 'dust/dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'
noDustSKIRTPath += 'noDust/'

textPath += args.galaxy+'/'
SKIRTPath += 'numPhotons'+args.numPhotons+'/inc'+args.inc+'/az'+args.az+'/'+args.SSP+'/'+args.galaxy+'/'
noDustSKIRTPath += 'numPhotons'+args.numPhotons+'/inc'+args.inc+'/az'+args.az+'/'+args.SSP+'/'+args.galaxy+'/'

start = timer()

# Including dust
if os.path.isfile(SKIRTPath+'sph_broadband_total.fits'):
    print('SKIRT dust SED already exists')
else:
    print('Generating dust SKIRT SED')
    os.system('mkdir -p '+SKIRTPath)
    # copy dust and radiation text files to SKIRT directory
    os.system('cp '+textPath+'gas.txt '+SKIRTPath+'gas.txt')
    os.system('cp '+textPath+'stars.txt '+SKIRTPath+'stars.txt')
    if eval(args.SF):
        os.system('cp '+textPath+'youngStars.txt '+SKIRTPath+'youngStars.txt')
    else:
        os.system('touch '+SKIRTPath+'youngStars.txt') # create empty text file
    # move ski file to SKIRT directory
    os.system('cp '+codePath+'resources/sph_template.ski '+SKIRTPath+'sph.ski')
    # calculate size of galaxy image from text files
    gas = np.loadtxt(textPath+'gas.txt')
    stars = np.loadtxt(textPath+'stars.txt')
    
    xLengthGas = (np.amax(gas[:,0]) - np.amin(gas[:,0]))
    yLengthGas = (np.amax(gas[:,1]) - np.amin(gas[:,1]))
    zLengthGas = (np.amax(gas[:,2]) - np.amin(gas[:,2]))
    xLengthStars = (np.amax(stars[:,0]) - np.amin(stars[:,0]))
    yLengthStars = (np.amax(stars[:,1]) - np.amin(stars[:,1]))
    zLengthStars = (np.amax(stars[:,2]) - np.amin(stars[:,2]))
    maxLength = np.amax([xLengthGas, yLengthGas, zLengthGas, xLengthStars, yLengthStars, zLengthStars])
    
    # change values in newly created .ski file to argparse values
    os.system('python '+codePath+'python/modify_ski.py --filePath='+SKIRTPath+'sph.ski --inc='+args.inc+' --az='+args.az+' --numPhotons='+args.numPhotons+' --pixels='+args.pixels+' --size='+str(maxLength)+' --dustFraction='+args.dustFraction+' --maxTemp='+args.maxTemp+' --SSP='+args.SSP)
    
    # go to SKIRT directory and run, then cd back
    os.chdir(SKIRTPath)
    os.system('skirt sph.ski')
    os.system('python -m pts.do plot_density_cuts .')
    
    # delete radiation and dust text files
    os.system('rm stars.txt')
    os.system('rm gas.txt')
    os.system('rm youngStars.txt')

# No dust
if os.path.isfile(noDustSKIRTPath+'sph_broadband_total.fits'):
    print('SKIRT no dust SED already exists')
else:
    print('Generating no dust SKIRT SED')
    os.system('mkdir -p '+noDustSKIRTPath)
    os.system('touch '+noDustSKIRTPath+'gas.txt') # create empty text file
    # copy stars text files to SKIRT directory
    os.system('cp '+textPath+'stars.txt '+noDustSKIRTPath+'stars.txt')
    if eval(args.SF):
        os.system('cp '+textPath+'youngStars_f_PDR0.txt '+noDustSKIRTPath+'youngStars.txt')
    else:
        os.system('touch '+noDustSKIRTPath+'youngStars.txt') # create empty text file
    # move ski file to SKIRT directory
    os.system('cp '+codePath+'resources/sph_template.ski '+noDustSKIRTPath+'sph.ski')
    
    # calculate size of galaxy image from text files
    stars = np.loadtxt(noDustSKIRTPath+'stars.txt')
    
    xLengthStars = (np.amax(stars[:,0]) - np.amin(stars[:,0]))
    yLengthStars = (np.amax(stars[:,1]) - np.amin(stars[:,1]))
    zLengthStars = (np.amax(stars[:,2]) - np.amin(stars[:,2]))
    maxLength = np.amax([xLengthStars, yLengthStars, zLengthStars])
    
    # change values in newly created .ski file to argparse values
    os.system('python '+codePath+'python/modify_ski.py --filePath='+noDustSKIRTPath+'sph.ski --inc='+args.inc+' --az='+args.az+' --numPhotons='+args.numPhotons+' --pixels='+args.pixels+' --size='+str(maxLength)+' --dustFraction='+args.dustFraction+' --maxTemp='+args.maxTemp+' --SSP='+args.SSP)
    
    # go to SKIRT directory and run, then cd back
    os.chdir(noDustSKIRTPath)
    os.system('skirt sph.ski')
    os.system('python -m pts.do plot_density_cuts .')
    
    # delete radiation and dust text files
    os.system('rm stars.txt')
    os.system('rm gas.txt')
    os.system('rm youngStars.txt')

end = timer()
time_SKIRT = end - start
time_SKIRT = str(datetime.timedelta(seconds=time_SKIRT))
print('Time to get SKIRT SED:', time_SKIRT)  

