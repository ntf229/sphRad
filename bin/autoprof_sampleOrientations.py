import argparse
import os
import shutil
from os.path import expanduser
from timeit import default_timer as timer
import subprocess
import datetime
import numpy as np
import sys
import xml.etree.ElementTree as ET
from copy import deepcopy

parser = argparse.ArgumentParser()
parser.add_argument("--num") # number of orientations 
parser.add_argument("--galaxy") # name of galaxy 
args = parser.parse_args()

num = int(args.num)-1 # first one is inc=0 az=0, already included in original ski file

origDir = os.getcwd()
codePath=expanduser('~')+'/sphRad/'
resultPath = '/scratch/ntf229/sphRad/' # store results here

# Directory structure stores important parameters
textPath = resultPath+'resources/NIHAO/TextFiles/noAgeSmooth/noSF/noClumps/'
SKIRTPath = resultPath+'resources/sampleOrientations_SKIRT/'

textPath += args.galaxy+'/'
SKIRTPath += args.galaxy+'/'

start = timer()

# sample orientations uniformly on the sphere 
a = np.random.uniform(low=-1.0, high=1.0, size=num)
inc = np.arccos(a) * 180 / np.pi 
az = np.random.uniform(low=0.0, high=360.0, size=num)

# round inc and az to 2 decimals
inc = np.round_(inc, decimals = 2)
az = np.round_(az, decimals = 2)

# calculate size of galaxy image from text files
stars = np.loadtxt(textPath+'stars.txt') 
xLengthStars = (np.amax(stars[:,0]) - np.amin(stars[:,0]))
yLengthStars = (np.amax(stars[:,1]) - np.amin(stars[:,1]))
zLengthStars = (np.amax(stars[:,2]) - np.amin(stars[:,2]))
maxLength = np.amax([xLengthStars, yLengthStars, zLengthStars])

os.system('mkdir -p '+SKIRTPath)
# copy stars text files to SKIRT directory
os.system('cp '+textPath+'stars.txt '+SKIRTPath+'stars.txt')
# move ski file to SKIRT directory
os.system('cp '+codePath+'resources/sampleOrientations_sph_template.ski '+SKIRTPath+'sph.ski')
# create instruments with sampled inc and az 
tree = ET.parse(SKIRTPath+'sph.ski')
root = tree.getroot()
for child in root.iter('FullInstrument'): 
    fullBB = deepcopy(child)
for child in root.iter('instruments'):
    for i in range(num):
        fullBB.set('inclination', str(inc[i])+' deg')
        fullBB.set('azimuth', str(az[i])+' deg')
        fullBB.set('instrumentName', 'inc'+str(inc[i])+'_az'+str(az[i]))
        child.insert(0,deepcopy(fullBB))
tree.write(SKIRTPath+'sph.ski', encoding='UTF-8', xml_declaration=True)
# change parameter  values in newly created .ski file 
os.system('python '+codePath+'python/sampleOrientations_modify_ski.py --filePath='+
          SKIRTPath+'sph.ski --size='+str(maxLength))
# go to SKIRT directory and run, then cd back
os.chdir(SKIRTPath)
os.system('skirt sph.ski')        
# delete radiation text files
os.system('rm stars.txt')

# Axis Ratios
print('starting autoprof')
for i in range(num+1): # including inc=0 az=0 
    if i == 0:
        name = 'inc0_az0'
        current_inc = 0
        current_az = 0
        if os.path.isfile(SKIRTPath+'axisRatios.npy'):
            continue # don't need to re-run inc=0 az=0 orientation
    else:
        name = 'inc'+str(inc[i-1])+'_az'+str(az[i-1])
        current_inc = inc[i-1]
        current_az = az[i-1]
    os.system('mkdir -p '+SKIRTPath+'autoprof/'+name+'/')
    os.system('python '+codePath+'python/orientations_mkprof.py --galaxy='+args.galaxy+
              ' --path='+SKIRTPath+' --save='+SKIRTPath+'autoprof/'+name+'/ --name='+name)
    axisRatioFile = np.loadtxt(SKIRTPath+'autoprof/'+name+'/'+args.galaxy+'-r-ba.csv', dtype=str)
    axisRatio = float(axisRatioFile[1,1].split(',')[0])
    os.system('rm '+SKIRTPath+'autoprof/'+name+'/galaxy-r.fits')
    os.system('rm '+SKIRTPath+'autoprof/'+name+'/photometry_ellipse_galaxy-r.jpg')
    #os.system('rm '+SKIRTPath+'autoprof/'+name+'/initialize_ellipse_galaxy-r.jpg')
    os.system('rm '+SKIRTPath+'sph_'+name+'_sed.dat')
    os.system('rm '+SKIRTPath+'sph_'+name+'_total.fits')
    if os.path.isfile(SKIRTPath+'axisRatios.npy'):
        axisRatios = np.load(SKIRTPath+'axisRatios.npy')
        axisRatios = np.vstack([axisRatios, [current_inc, current_az, axisRatio]])
        np.save(SKIRTPath+'axisRatios.npy', axisRatios)
    else:
        axisRatios = np.asarray([current_inc, current_az, axisRatio])
        np.save(SKIRTPath+'axisRatios.npy', axisRatios)

end = timer()
time_SKIRT = end - start
time_SKIRT = str(datetime.timedelta(seconds=time_SKIRT))
print('Time to finish:', time_SKIRT)  

