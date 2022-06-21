import numpy as np
import matplotlib.pyplot as plt

# All SecondBatch galaxies (69)
#galaxies = ['g1.12e12','g1.92e12','g2.39e11','g2.79e12','g3.23e11','g3.49e11','g3.59e11','g3.61e11','g5.31e11',
#            'g5.36e11','g5.38e11','g5.55e11','g7.08e11','g7.44e11','g8.26e11','g8.28e11','g1.05e11','g1.23e10',
#            'g1.50e10','g1.88e10','g1.95e10','g2.37e10','g2.57e11','g2.83e10','g3.21e11','g4.36e09','g4.99e09',
#            'g5.84e09','g7.12e10','g9.59e10','g1.08e11','g1.37e11','g1.52e11','g1.89e10','g2.04e11','g2.39e10',
#            'g2.63e10','g2.94e10','g3.44e10','g4.48e10','g5.22e09','g6.12e10','g7.34e09','g1.09e10','g1.44e10',
#            'g1.57e10','g1.90e10','g2.09e10','g2.41e11','g2.64e10','g3.06e11','g3.67e10','g4.86e10','g5.41e09',
#            'g6.96e10','g8.89e10','g1.18e10','g1.47e10','g1.64e11','g1.92e10','g2.34e10','g2.54e11','g2.80e10',
#            'g3.19e10','g3.93e10','g4.94e10','g5.59e09','g7.05e09','g9.26e09']

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
	textFilePath = '/scratch/ntf229/RT_fit/resources/NIHAO/TextFiles/'+galaxies[i]+'/'
	stars = np.loadtxt(textFilePath+'stars.txt')
	starMasses = stars[:,7] # units of M_sun
	trueMass = np.sum(starMasses)
	starAges = stars[:,9] # units of years

	plt.figure(figsize=(10,8))

	counts, bins = np.histogram(starAges,bins=30,weights=starMasses,density=True)
	plt.hist(bins[:-1], bins, weights=counts*trueMass)

	plt.xlabel('Age (years)', fontsize=16)
	plt.ylabel('SFH ('+r'$M_{\odot}$'+'/year)', fontsize=16)
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.savefig('/scratch/ntf229/RT_fit/resources/NIHAO/SFHs/'+galaxies[i]+'_SFH.png',dpi=300)
	plt.close()
