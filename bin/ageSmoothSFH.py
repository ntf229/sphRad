import numpy as np
import matplotlib.pyplot as plt
import os

codePath = '/home/ntf229/sphRad/'
resultPath = '/scratch/ntf229/sphRad/' # store results here

# Directory structure stores important parameters
noAgeSmoothPath = resultPath+'resources/NIHAO/Particles/noAgeSmooth/noSF/noClumps/'
ageSmoothPath = resultPath+'resources/NIHAO/Particles/ageSmooth/noSF/noClumps/'
plotPath = resultPath+'resources/NIHAO/ageSmoothSFHPlots/'

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

for i in range(len(galaxies)):
    noAgeSmoothStars = np.load(noAgeSmoothPath+galaxies[i]+'/stars.npy')
    ageSmoothStars = np.load(ageSmoothPath+galaxies[i]+'/stars.npy')
    noAgeSmoothMasses = noAgeSmoothStars[:,7] # in M_sun
    ageSmoothMasses = ageSmoothStars[:,7] # in M_sun
    noAgeSmoothTotMass = np.sum(noAgeSmoothMasses) # in M_sun
    ageSmoothTotMass = np.sum(ageSmoothMasses) # in M_sun
    noAgeSmoothAges = noAgeSmoothStars[:,9]/1e9 # in Gyrs
    ageSmoothAges = ageSmoothStars[:,9]/1e9 # in Gyrs
    # SFH grids
    binGrid1 = np.linspace(0,14,num=100)
    yearsPerBin1 = 1.4e10 / 100
    binGrid2 = np.linspace(0,14,num=1000)
    yearsPerBin2 = 1.4e10 / 1000
    # Make histograms    
    plt.figure(figsize=(10,8))
    #counts, bins = np.histogram(noAgeSmoothAges, bins=binGrid1, weights=noAgeSmoothMasses, density=True)
    #plt.hist(bins[:-1], bins, weights=counts*noAgeSmoothTotMass, alpha=1, label='Before Smoothing', color='red')
    #newCounts, newBins = np.histogram(ageSmoothAges, bins=binGrid2, weights=ageSmoothMasses, density=True)
    #plt.hist(newBins[:-1], newBins, weights=newCounts*ageSmoothTotMass, alpha=0.5, label='After Smoothing', color='blue')
    counts, bins = np.histogram(noAgeSmoothAges, bins=binGrid1, weights=noAgeSmoothMasses/yearsPerBin1, density=False)
    plt.hist(bins[:-1], bins, weights=counts, alpha=1, label='Before Smoothing', color='red')
    newCounts, newBins = np.histogram(ageSmoothAges, bins=binGrid2, weights=ageSmoothMasses/yearsPerBin2, density=False)
    plt.hist(newBins[:-1], newBins, weights=newCounts, alpha=0.5, label='After Smoothing', color='blue')
    plt.xlabel('Age [Gyrs]', fontsize=28)
    plt.ylabel('SFH ['+r'$M_{\odot} \, / \, year$]', fontsize=28)
    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)
    plt.yscale('log')
    plt.legend(fontsize=20)
    plt.savefig(plotPath+galaxies[i]+'_SFHs.png', dpi=300, bbox_inches='tight', pad_inches=0.25)
    plt.close()


print('done')
