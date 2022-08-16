# selects and stores orientations for each galaxy that span the range of axis ratios

import numpy as np
import os
import matplotlib.pyplot as plt

codePath = '/home/ntf229/sphRad/'
resultPath = '/scratch/ntf229/sphRad/' # store results here

# Directory structure stores important parameters
dataPath = resultPath+'resources/sampleOrientations_SKIRT/'
selectedPath = resultPath+'resources/selectedOrientations/'
plotPath = resultPath+'resources/Plots/AxisRatioDist/'
massPath = resultPath+'resources/NIHAO/GlobalProps/stellarMasses/'

os.system('mkdir -p '+plotPath)

# Final sample (65) excluding g3.19e10 (conatins two galaxies)
galaxies = ['g1.88e10','g1.89e10','g1.90e10','g2.34e10','g2.63e10','g2.64e10','g2.80e10','g2.83e10','g2.94e10',
            'g3.44e10','g3.67e10','g3.93e10','g4.27e10','g4.48e10','g4.86e10','g4.94e10','g4.99e10',
            'g5.05e10','g6.12e10','g6.37e10','g6.77e10','g6.91e10','g7.12e10','g8.89e10','g9.59e10','g1.05e11',
            'g1.08e11','g1.37e11','g1.52e11','g1.57e11','g1.59e11','g1.64e11','g2.04e11','g2.19e11','g2.39e11',
            'g2.41e11','g2.42e11','g2.54e11','g3.06e11','g3.21e11','g3.23e11','g3.49e11','g3.55e11','g3.59e11',
            'g3.71e11','g4.90e11','g5.02e11','g5.31e11','g5.36e11','g5.38e11','g5.46e11','g5.55e11','g6.96e11',
            'g7.08e11','g7.44e11','g7.55e11','g7.66e11','g8.06e11','g8.13e11','g8.26e11','g8.28e11','g1.12e12',
            'g1.77e12','g1.92e12','g2.79e12']

selectedInc = np.zeros((len(galaxies), 10)) 
selectedAz = np.zeros((len(galaxies), 10)) 
selectedAxisRatio = np.zeros((len(galaxies), 10)) 

plt.figure(figsize=(10,8))
for i in range(len(galaxies)):
    axisRatios = np.load(dataPath+galaxies[i]+'/axisRatios.npy')[:,2]
    inc = np.load(dataPath+galaxies[i]+'/axisRatios.npy')[:,0]
    az = np.load(dataPath+galaxies[i]+'/axisRatios.npy')[:,1]
    stellarMass = float(np.load(massPath+galaxies[i]+'/stellarMass.npy'))
    plt.scatter(np.zeros(len(axisRatios))+stellarMass, axisRatios, marker='o', s=5, c='k')
    axisRatioBins = np.linspace(np.amin(axisRatios[axisRatios > 0.11]), 1, num=10)
    for j in range(len(axisRatioBins)):
        index = np.abs(axisRatios - axisRatioBins[j]).argmin()
        selectedInc[i,j] = inc[index]
        selectedAz[i,j] = az[index]
        selectedAxisRatio[i,j] = np.round_(axisRatios[index], decimals = 4)
        # remove selections to avoid duplicates
        inc = np.delete(inc, index)
        az = np.delete(az, index) 
        axisRatios = np.delete(axisRatios, index) 
    # save selections as numpy array 
    os.system('mkdir -p '+selectedPath+galaxies[i]+'/')
    np.save(selectedPath+galaxies[i]+'/selectedAxisRatios.npy', 
            np.c_[selectedInc[i,:], selectedAz[i,:], selectedAxisRatio[i,:]])

plt.xlabel('Stellar Mass ('+r'$M_{\odot}$)', fontsize=16)
plt.ylabel('Axis Ratio', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xscale('log')
plt.savefig(plotPath+'axisRatioDist.png',dpi=300)
plt.close()

# copy selected autoprof images for visual inspection
#for i in range(len(galaxies)):
#    os.system('mkdir -p '+plotPath+'Selected/'+galaxies[i])
#    for j in range(len(axisRatioBins)):
#        name = 'inc'+str(selectedInc[i,j])+'_az'+str(selectedAz[i,j])
#        if selectedInc[i,j] == 0:
#            if selectedAz[i,j] == 0:
#                name = 'inc0_az0'
#        os.system('cp '+dataPath+galaxies[i]+'/autoprof/'+name+'/fit_ellipse_galaxy-r.jpg '+
#        plotPath+'Selected/'+galaxies[i]+'/axisRatio'+str(selectedAxisRatio[i,j])+'.jpg')
    
    
