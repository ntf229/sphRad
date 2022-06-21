import numpy as np
import matplotlib.pyplot as plt
import os

partialDust = True

# Final sample (66)
galaxies = ['g1.88e10','g1.89e10','g1.90e10','g2.34e10','g2.63e10','g2.64e10','g2.80e10','g2.83e10','g2.94e10',
            'g3.19e10','g3.44e10','g3.67e10','g3.93e10','g4.27e10','g4.48e10','g4.86e10','g4.94e10','g4.99e10',
            'g5.05e10','g6.12e10','g6.37e10','g6.77e10','g6.91e10','g7.12e10','g8.89e10','g9.59e10','g1.05e11',
            'g1.08e11','g1.37e11','g1.52e11','g1.57e11','g1.59e11','g1.64e11','g2.04e11','g2.19e11','g2.39e11',
            'g2.41e11','g2.42e11','g2.54e11','g3.06e11','g3.21e11','g3.23e11','g3.49e11','g3.55e11','g3.59e11',
            'g3.71e11','g4.90e11','g5.02e11','g5.31e11','g5.36e11','g5.38e11','g5.46e11','g5.55e11','g6.96e11',
            'g7.08e11','g7.44e11','g7.55e11','g7.66e11','g8.06e11','g8.13e11','g8.26e11','g8.28e11','g1.12e12',
            'g1.77e12','g1.92e12','g2.79e12']

if partialDust:
    path = '/scratch/ntf229/RT_fit/resources/SKIRT_partialDust/'
else:
    path = '/scratch/ntf229/RT_fit/resources/SKIRT/'
path2 = '/maxLevel13/wavelengths601/numPhotons1e9/inc0/dust/dustFraction0.1/maxTemp16000/'

os.system('mkdir -p '+path+'SEDs')

for i in range(len(galaxies)):

	plt.figure(figsize=(10,8))

	wave_correct = np.load(path+galaxies[i]+path2+'wave.npy') # units of microns
	wave_correct = wave_correct * 1e4 # convert to Angstroms
	spec_correct = np.load(path+galaxies[i]+path2+'spec.npy') # units of Jy
	spec_correct = spec_correct / 3631 # convert to maggies
	
	mask_correct = (wave_correct >= 1e3) & (wave_correct <= 1e7)
	
	plt.plot(wave_correct[mask_correct], spec_correct[mask_correct], label="Actual Spectrum",color='blue',alpha=1)
	plt.xlabel('Wavelength ('+r'$\AA$'+')', fontsize=16)
	plt.ylabel('Flux (Maggies)', fontsize=16)
	
	plt.xscale('log')
	plt.yscale('log')
	
	#plt.legend(fontsize=16)
	
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.savefig(path+'SEDs/'+galaxies[i]+'_SED.png', dpi=300)
	plt.close()

print('done')
