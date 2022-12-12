# selects and stores orientations for each galaxy that span the range of axis ratios

import numpy as np
import os

codePath = '/home/ntf229/sphRad/'
resultPath = '/scratch/ntf229/sphRad/' # store results here

# Directory structure stores important parameters
sampledPath = resultPath+'resources/sampledAxisRatios/'
selectedPath = resultPath+'resources/selectedOrientations/'

# Final sample (65) excluding g3.19e10 (conatins two galaxies)
galaxies = ['g1.88e10','g1.89e10','g1.90e10','g2.34e10','g2.63e10','g2.64e10','g2.80e10','g2.83e10','g2.94e10',
            'g3.44e10','g3.67e10','g3.93e10','g4.27e10','g4.48e10','g4.86e10','g4.94e10','g4.99e10',
            'g5.05e10','g6.12e10','g6.37e10','g6.77e10','g6.91e10','g7.12e10','g8.89e10','g9.59e10','g1.05e11',
            'g1.08e11','g1.37e11','g1.52e11','g1.57e11','g1.59e11','g1.64e11','g2.04e11','g2.19e11','g2.39e11',
            'g2.41e11','g2.42e11','g2.54e11','g3.06e11','g3.21e11','g3.23e11','g3.49e11','g3.55e11','g3.59e11',
            'g3.71e11','g4.90e11','g5.02e11','g5.31e11','g5.36e11','g5.38e11','g5.46e11','g5.55e11','g6.96e11',
            'g7.08e11','g7.44e11','g7.55e11','g7.66e11','g8.06e11','g8.13e11','g8.26e11','g8.28e11','g1.12e12',
            'g1.77e12','g1.92e12','g2.79e12']

sampledIncAzAR = np.load(sampledPath+'sampledIncAzAR.npy')
#selectedIncAzAr = np.zeros(sampledIncAzAR.shape)

for g in range(len(galaxies)):
    selectedIncAzAr = np.zeros((10, 3))
    os.system('mkdir -p '+selectedPath+galaxies[g]+'/') 
    incs = sampledIncAzAR[g,:,0]
    azs = sampledIncAzAR[g,:,1]
    axisRatios = sampledIncAzAR[g,:,2]
    axisRatioBins = np.linspace(np.amin(axisRatios), 1, num=10)
    for j in range(len(axisRatioBins)):
        # find axis ratios closest to bins
        index = np.abs(axisRatios - axisRatioBins[j]).argmin()
        #selectedIncAzAr[g,j,0] = incs[index]
        #selectedIncAzAr[g,j,1] = azs[index]
        #selectedIncAzAr[g,j,2] = axisRatios[index]
        selectedIncAzAr[j,0] = incs[index]
        selectedIncAzAr[j,1] = azs[index]
        selectedIncAzAr[j,2] = axisRatios[index]
        # remove selections to avoid duplicates
        incs = np.delete(incs, index)
        azs = np.delete(azs, index) 
        axisRatios = np.delete(axisRatios, index) 
        # copy ellipse images
        incStr = str(selectedIncAzAr[j,0])
        azStr = str(selectedIncAzAr[j,1])
        if incStr == '0.0':
            incStr = '0'
        if azStr == '0.0':
            azStr = '0'
        axisRatioStr = format(selectedIncAzAr[j,2], '.4f')
        os.system('cp '+sampledPath+galaxies[g]+'/axisRatio'+axisRatioStr+'_inc'+incStr+'_az'+
                    azStr+'_ellipse.png '+selectedPath+galaxies[g]+'/')
    # save selections as numpy array
    np.save(selectedPath+galaxies[g]+'/selectedIncAzAR.npy', selectedIncAzAr)


