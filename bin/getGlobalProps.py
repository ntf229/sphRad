# store global properties of NIHAO galaxies as numpy array 

import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--dustFraction") # dust to metal ratio (SKIRT parameter)
parser.add_argument("--maxTemp") # maximum temperature at which dust can form (SKIRT parameter)
args = parser.parse_args()

codePath = '/home/ntf229/sphRad/'
resultPath = '/scratch/ntf229/sphRad/' # store results here

# Directory structure stores important parameters
textPath = resultPath+'resources/NIHAO/TextFiles/noAgeSmooth/noSF/noClumps/'
plotPath = resultPath+'resources/NIHAO/GlobalProps/'

plotPath += 'dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'

os.system('mkdir -p '+plotPath)

# Final sample (66)
galaxies = ['g1.88e10','g1.89e10','g1.90e10','g2.34e10','g2.63e10','g2.64e10','g2.80e10','g2.83e10','g2.94e10',
            'g3.19e10','g3.44e10','g3.67e10','g3.93e10','g4.27e10','g4.48e10','g4.86e10','g4.94e10','g4.99e10',
            'g5.05e10','g6.12e10','g6.37e10','g6.77e10','g6.91e10','g7.12e10','g8.89e10','g9.59e10','g1.05e11',
            'g1.08e11','g1.37e11','g1.52e11','g1.57e11','g1.59e11','g1.64e11','g2.04e11','g2.19e11','g2.39e11',
            'g2.41e11','g2.42e11','g2.54e11','g3.06e11','g3.21e11','g3.23e11','g3.49e11','g3.55e11','g3.59e11',
            'g3.71e11','g4.90e11','g5.02e11','g5.31e11','g5.36e11','g5.38e11','g5.46e11','g5.55e11','g6.96e11',
            'g7.08e11','g7.44e11','g7.55e11','g7.66e11','g8.06e11','g8.13e11','g8.26e11','g8.28e11','g1.12e12',
            'g1.77e12','g1.92e12','g2.79e12']

params = ['Name', 'Stellar Mass', 'Gas Mass', 'Dust Mass', 'SFR']

# initialize numpy array to store global parameters
GlobalProps = np.empty((len(galaxies), len(params)), dtype=object)

for i in range(len(galaxies)):
    stars = np.loadtxt(textPath+galaxies[i]+'/stars.txt')
    gas = np.loadtxt(textPath+galaxies[i]+'/gas.txt')
    starMasses = stars[:,7] # units of M_sun
    stellarMass = np.sum(starMasses)
    starAges = stars[:,9] # units of years
    gasMasses = gas[:,4] # units M_sun
    gasMass = np.sum(gasMasses)
    gasTemp = gas[:,6]
    gasMetals = gas[:,5]
    dustMask = gasTemp < float(args.maxTemp)
    dustMass = np.sum(gasMasses[dustMask] * gasMetals[dustMask]) * float(args.dustFraction)
    ageMask = starAges < 1.e8 # less than 100 Myrs old
    SFR = np.sum(starMasses[ageMask]) / 1.e8 # M_sun / year averaged over 100 Myrs
    GlobalProps[i,0] = galaxies[i]
    GlobalProps[i,1] = str(stellarMass)
    GlobalProps[i,2] = str(gasMass)
    GlobalProps[i,3] = str(dustMass)
    GlobalProps[i,4] = str(SFR)

np.save(plotPath+'GlobalProps.npy', GlobalProps)

print('done')
