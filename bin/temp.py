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
parser.add_argument("--dust") # include dust; True or False
parser.add_argument("--inc") # inclination angle (SKIRT parameter)
parser.add_argument("--numPhotons") # number of photon packages (SKIRT parameter)
parser.add_argument("--pixels") # number of pixels (square) for image (SKIRT parameter)
parser.add_argument("--dustFraction") # dust to metal ratio (SKIRT parameter)
parser.add_argument("--maxTemp") # maximum temperature at which dust can form (SKIRT parameter)
parser.add_argument("--galaxy") # name of galaxy 
args = parser.parse_args()

origDir = os.getcwd()
codePath=expanduser('~')+'/sphRad/'
resultPath = '/scratch/ntf229/sphRad/' # store results here

# Directory structure stores important parameters
textPath = resultPath+'resources/NIHAO/TextFiles/'
SKIRTPath = resultPath+'resources/SKIRT/'
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
SKIRTPath += 'numPhotons'+args.numPhotons+'/inc'+args.inc+'/'+args.galaxy+'/'

os.chdir(SKIRTPath)
os.system('python -m pts.do plot_density_cuts .')

print('done')

