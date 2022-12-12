import numpy as np
import matplotlib.pyplot as plt
import os

codePath = '/home/ntf229/sphRad/'
resultPath = '/scratch/ntf229/sphRad/' # store results here

# Directory structure stores important parameters
noAgeSmoothPath = resultPath+'resources/NIHAO/Particles/noAgeSmooth/noSF/noClumps/'
ageSmoothPath = resultPath+'resources/NIHAO/Particles/ageSmooth/noSF/noClumps/'
plotPath = resultPath+'resources/NIHAO/ageSmoothSFRPlots/'

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

noAgeSmoothStellarMass = np.zeros(len(galaxies))
ageSmoothStellarMass = np.zeros(len(galaxies))
noAgeSmoothSFR100 = np.zeros(len(galaxies))
ageSmoothSFR100 = np.zeros(len(galaxies))
noAgeSmoothSFR10 = np.zeros(len(galaxies))
ageSmoothSFR10 = np.zeros(len(galaxies))

if os.path.isfile(plotPath+'ageSmoothSFR10.npy'):
    print('loading masses and SFRs')
    noAgeSmoothStellarMass = np.load(plotPath+'noAgeSmoothStellarMass.npy')
    ageSmoothStellarMass = np.load(plotPath+'ageSmoothStellarMass.npy')
    noAgeSmoothSFR100 = np.load(plotPath+'noAgeSmoothSFR100.npy')
    ageSmoothSFR100 = np.load(plotPath+'ageSmoothSFR100.npy')
    noAgeSmoothSFR10 = np.load(plotPath+'noAgeSmoothSFR10.npy')
    ageSmoothSFR10 = np.load(plotPath+'ageSmoothSFR10.npy')
else:
    print('loading particle data')
    for i in range(len(galaxies)):
        noAgeSmoothStars = np.load(noAgeSmoothPath+galaxies[i]+'/stars.npy')
        ageSmoothStars = np.load(ageSmoothPath+galaxies[i]+'/stars.npy')
        noAgeSmoothStellarMasses = noAgeSmoothStars[:,7] # in M_sun
        ageSmoothStellarMasses = ageSmoothStars[:,7] # in M_sun
        noAgeSmoothStellarMass[i] = np.sum(noAgeSmoothStellarMasses) # in M_sun
        ageSmoothStellarMass[i] = np.sum(ageSmoothStellarMasses) # in M_sun
        noAgeSmoothAges = noAgeSmoothStars[:,9] # in years
        ageSmoothAges = ageSmoothStars[:,9] # in years
        noAgeSmoothSFR100[i] = np.sum(noAgeSmoothStellarMasses[noAgeSmoothAges < 1.e8]) / 1.e8 # M_sun / year (averaged over 100 Myrs)
        ageSmoothSFR100[i] = np.sum(ageSmoothStellarMasses[ageSmoothAges < 1.e8]) / 1.e8 # M_sun / year (averaged over 100 Myrs)
        noAgeSmoothSFR10[i] = np.sum(noAgeSmoothStellarMasses[noAgeSmoothAges < 1.e7]) / 1.e7 # M_sun / year (averaged over 10 Myrs)
        ageSmoothSFR10[i] = np.sum(ageSmoothStellarMasses[ageSmoothAges < 1.e7]) / 1.e7 # M_sun / year (averaged over 10 Myrs)
    np.save(plotPath+'noAgeSmoothStellarMass.npy', noAgeSmoothStellarMass)
    np.save(plotPath+'ageSmoothStellarMass.npy', ageSmoothStellarMass)
    np.save(plotPath+'noAgeSmoothSFR100.npy', noAgeSmoothSFR100)
    np.save(plotPath+'ageSmoothSFR100.npy', ageSmoothSFR100)
    np.save(plotPath+'noAgeSmoothSFR10.npy', noAgeSmoothSFR10)
    np.save(plotPath+'ageSmoothSFR10.npy', ageSmoothSFR10)

plt.figure(figsize=(10,8))
plt.scatter(np.log10(noAgeSmoothStellarMass), noAgeSmoothSFR100/noAgeSmoothStellarMass, color='red', label='Original')
plt.scatter(np.log10(ageSmoothStellarMass), ageSmoothSFR100/ageSmoothStellarMass, color='blue', label='Smoothed Ages')
plt.legend(fontsize=16)
plt.xlabel('log(Stellar Mass / '+r'$M_{\odot}$)', fontsize=16)
plt.ylabel('sSFR / '+r'$yr^{-1}$ Averaged Over 100 Myrs', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig(plotPath+'sSFR100_stellarMass.png',dpi=300)
plt.close()

plt.figure(figsize=(10,8))
plt.scatter(np.log10(noAgeSmoothStellarMass), noAgeSmoothSFR10/noAgeSmoothStellarMass, color='red', label='Original')
plt.scatter(np.log10(ageSmoothStellarMass), ageSmoothSFR10/ageSmoothStellarMass, color='blue', label='Smoothed Ages')
plt.legend(fontsize=16)
plt.xlabel('log(Stellar Mass / '+r'$M_{\odot}$)', fontsize=16)
plt.ylabel('sSFR / '+r'$yr^{-1}$ Averaged Over 10 Myrs', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig(plotPath+'sSFR10_stellarMass.png',dpi=300)
plt.close()

print('done')
