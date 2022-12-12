import argparse
import os
from os.path import expanduser
from timeit import default_timer as timer
import numpy as np
import sys
import os
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from matplotlib.patches import Ellipse
import datetime

def getData(galaxy, fileName):
	inc = fileName.split('inc')[1].split('_az')[0]
	az = fileName.split('az')[1].split('_total')[0]
	images_file = fits.open(SKIRTPath+galaxies[g]+'/'+fileName)
	images = np.asarray(images_file[0].data) # cube of broadbands in MJy/sr
	r = images[0, :, :] # r band
	return inc, az, r

def plotEllipse(galaxy, inc, az, r, center, phi, ba):
	os.system('mkdir -p '+savePath+galaxy+'/')
	fig = plt.figure()
	sizes = np.shape(r)
	fig.set_size_inches(1. * sizes[0] / sizes[1], 1, forward = False)
	ax = plt.Axes(fig, [0., 0., 1., 1.])
	#ax.set_axis_off()
	fig.add_axes(ax)
	ax.imshow(r,interpolation='none', origin='lower')
	ax.add_patch(Ellipse((center[1], center[0]), width=100, height=100*ba, angle=90 - 1*phi,
		edgecolor='black',
		facecolor='none',
		linewidth=0.2))
	axisRatioStr = format(ba, '.4f')
	plt.savefig(savePath+galaxy+'/axisRatio'+axisRatioStr+'_inc'+inc+'_az'+az+'_ellipse.png', dpi=sizes[0])
	plt.close()

parser = argparse.ArgumentParser()
parser.add_argument("--num") # number of orientations 
args = parser.parse_args()

num = int(args.num)-1 # first one is inc=0 az=0, already included in original ski file

origDir = os.getcwd()
codePath=expanduser('~')+'/sphRad/'
resultPath = '/scratch/ntf229/sphRad/' # store results here

# Directory structure stores important parameters
SKIRTPath = resultPath+'resources/sampleOrientations_SKIRT/'
savePath = resultPath+'resources/sampledAxisRatios/'

start = timer()

# Final sample (65) excluding g3.19e10 (conatins two galaxies)
galaxies = ['g1.88e10','g1.89e10','g1.90e10','g2.34e10','g2.63e10','g2.64e10','g2.80e10','g2.83e10','g2.94e10',
			'g3.44e10','g3.67e10','g3.93e10','g4.27e10','g4.48e10','g4.86e10','g4.94e10','g4.99e10',
			'g5.05e10','g6.12e10','g6.37e10','g6.77e10','g6.91e10','g7.12e10','g8.89e10','g9.59e10','g1.05e11',
			'g1.08e11','g1.37e11','g1.52e11','g1.57e11','g1.59e11','g1.64e11','g2.04e11','g2.19e11','g2.39e11',
			'g2.41e11','g2.42e11','g2.54e11','g3.06e11','g3.21e11','g3.23e11','g3.49e11','g3.55e11','g3.59e11',
			'g3.71e11','g4.90e11','g5.02e11','g5.31e11','g5.36e11','g5.38e11','g5.46e11','g5.55e11','g6.96e11',
			'g7.08e11','g7.44e11','g7.55e11','g7.66e11','g8.06e11','g8.13e11','g8.26e11','g8.28e11','g1.12e12',
			'g1.77e12','g1.92e12','g2.79e12']

#galaxies = ['g1.37e11','g5.02e11','g7.66e11','g8.26e11','g2.79e12']

numCenters = 50 # get average position of numCenters brightest pixels
numOrientations = int(args.num)

stretch = 0.01 # decreasing stretch makes dim pixels brighter
stretch_image = 0.001

minCut = 0.01

incAzAR = np.zeros((len(galaxies), numOrientations, 3)) # inc, az, axis ratio

length = 250 # number of image pixels
positions = np.linspace(0, length-1, num=length)
xGrid = np.zeros((length, length), dtype=int)
yGrid = np.zeros((length, length), dtype=int)
for j in range(length):
	for k in range(length):
		xGrid[j,k] = int(positions[j])
		yGrid[j,k] = int(positions[k])

for g in range(len(galaxies)):
	fileNames = np.asarray(os.listdir(SKIRTPath+galaxies[g]+'/'))
	fitsMask = np.char.find(fileNames, '.fits') != -1
	fileNames = fileNames[fitsMask] # only includes fits files
	inc = []
	az = []
	for i in range(numOrientations):
		inc, az, r = getData(galaxies[g], fileNames[i])
		temp_r = r.copy()
		centers = np.zeros((numCenters, 2), dtype=int)
		for j in range(numCenters):
			centers[j] = np.unravel_index(np.argmax(temp_r), temp_r.shape)
			maxCut = temp_r[centers[j,0], centers[j,1]] # last loop will set maximum allowed value
			temp_r[centers[j,0], centers[j,1]] = 0 # so we don't get the same index over and over
		center = [int(np.mean(centers[:,0])), int(np.mean(centers[:,1]))]
		norm_r = r/np.amax(r)
		maxCut /= np.amax(r)
		norm_r[norm_r > maxCut] = maxCut
		norm_r = norm_r/np.amax(norm_r) # re-normalize
		norm_r[norm_r<minCut] = 0
		scaled_r = np.arcsinh(norm_r/stretch)
		image_r = np.arcsinh((r/np.amax(r))/stretch_image) # no cuts
		radius = np.zeros((length, length))
		# shift grid to be centered at center
		xGridShift = xGrid - center[0]
		yGridShift = yGrid - center[1]
		radius = np.sqrt(xGridShift**2 + yGridShift**2)
		centerMask = radius > 7
		Mxx = np.sum(scaled_r[centerMask] * xGridShift[centerMask]**2 / radius[centerMask]**2) / np.sum(scaled_r[centerMask])
		Mxy = np.sum(scaled_r[centerMask] * xGridShift[centerMask] * yGridShift[centerMask] / radius[centerMask]**2) / np.sum(scaled_r[centerMask])
		Myy = np.sum(scaled_r[centerMask] * yGridShift[centerMask]**2 / radius[centerMask]**2) / np.sum(scaled_r[centerMask])
		q = 2 * Mxx - 1
		u = 2 * Mxy
		phi = 0.5 * np.arctan2(u, q) * 180. / np.pi
		yy = q**2 + u**2
		ba = (1. + yy - 2. * np.sqrt(yy)) / (1 - yy)
		plotEllipse(galaxies[g], inc, az, image_r, center, phi, ba)
		incAzAR[g, i, 0] = float(inc)
		incAzAR[g, i, 1] = float(az)
		incAzAR[g, i, 2] = float(ba)

np.save(savePath+'sampledIncAzAR.npy', incAzAR)

end = timer()
time_SKIRT = end - start
time_SKIRT = str(datetime.timedelta(seconds=time_SKIRT))
print('Time to finish:', time_SKIRT)  

