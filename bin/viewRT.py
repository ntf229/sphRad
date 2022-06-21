# run SKIRT RT from text files (need to run makeTextFiles.py with same arguments first)
# uniformly samples viewing angles (num = args.orientations)

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
parser.add_argument("--dust") # include dust; True or False
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
SKIRTPath = resultPath+'resources/SKIRT_views/'
if eval(args.ageSmooth):
    textPath += 'ageSmooth/'
    SKIRTPath += 'ageSmooth/'
else:
    textPath += 'noAgeSmooth/'
    SKIRTPath += 'noAgeSmooth/'
if eval(args.SF):
    textPath += 'SF/tauClear'+args.tauClear+'/'
    SKIRTPath += 'SF/tauClear'+args.tauClear+'/'
else:
	textPath += 'noSF/'
	SKIRTPath += 'noSF/'
if eval(args.clumps):
    textPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
    SKIRTPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
else:
	textPath += 'noClumps/'
	SKIRTPath += 'noClumps/'
if eval(args.dust):
    SKIRTPath += 'dust/dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'
else:
    SKIRTPath += 'noDust/'

textPath += args.galaxy+'/'
SKIRTPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'+args.galaxy+'/'

start = timer()

if os.path.isfile(SKIRTPath+'sph_broadband_total.fits'):
  print('SKIRT generated SED already exists')
else:
  print('Generating SKIRT SED')
  os.system('mkdir -p '+SKIRTPath)
  # copy dust and radiation text files to SKIRT directory
  os.system('cp '+textPath+'gas.txt '+SKIRTPath+'gas.txt')
  os.system('cp '+textPath+'stars.txt '+SKIRTPath+'stars.txt')
  if eval(args.SF):
    os.system('cp '+textPath+'youngStars.txt '+SKIRTPath+'youngStars.txt')
  else:
    os.system('touch '+SKIRTPath+'youngStars.txt')
    print('Created empty youngStars.txt file')

  # move ski file to SKIRT directory
  os.system('cp '+codePath+'resources/1e2views.ski '+SKIRTPath+'sph.ski')

  # calculate size of galaxy image from text files
  gas = np.loadtxt(SKIRTPath+'gas.txt')
  stars = np.loadtxt(SKIRTPath+'stars.txt')

  xLengthGas = (np.amax(gas[:,0]) - np.amin(gas[:,0]))
  yLengthGas = (np.amax(gas[:,1]) - np.amin(gas[:,1]))
  zLengthGas = (np.amax(gas[:,2]) - np.amin(gas[:,2]))

  xLengthStars = (np.amax(stars[:,0]) - np.amin(stars[:,0]))
  yLengthStars = (np.amax(stars[:,1]) - np.amin(stars[:,1]))
  zLengthStars = (np.amax(stars[:,2]) - np.amin(stars[:,2]))

  maxLength = np.amax([xLengthGas, yLengthGas, zLengthGas, xLengthStars, yLengthStars, zLengthStars])

  print('maxLength:', maxLength)

  # change values in newly created .ski file to argparse values
  os.system('python '+codePath+'python/modify_ski_view.py --filePath='+SKIRTPath+'sph.ski --numPhotons='+args.numPhotons+' --pixels='+args.pixels+' --size='+str(maxLength)+' --dustFraction='+args.dustFraction+' --maxTemp='+args.maxTemp+' --SSP='+args.SSP)

  if eval(args.dust):
    print('Including dust')
  else:
    os.system('rm '+SKIRTPath+'gas.txt')
    os.system('touch '+SKIRTPath+'gas.txt')
    print('Created empty dust.txt file')

  # go to SKIRT directory and run, then cd back
  os.chdir(SKIRTPath)
  os.system('skirt sph.ski')
  os.system('python -m pts.do plot_density_cuts .')
  #os.system('python -m pts.do plot_seds .')
  #os.system('python -m pts.do make_images .')

  # delete radiation and dust text files
  os.system('rm stars.txt')
  os.system('rm gas.txt')
  os.system('rm youngStars.txt')

end = timer()
time_SKIRT = end - start
time_SKIRT = str(datetime.timedelta(seconds=time_SKIRT))
print('Time to get SKIRT SED:', time_SKIRT)  

