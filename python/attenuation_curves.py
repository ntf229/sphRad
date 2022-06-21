# plots attenuation curves of NIHAO galaxy spectra from SKIRT

import numpy as np
import matplotlib.pyplot as plt

plotPath = '/scratch/ntf229/RT_fit/plots/dustFraction0.2/maxTemp8000/'

# 16 galaxies selected by hand from the NIHAO Suite                                                                                                                             
#galaxies = ['g1.12e12','g1.92e12','g2.39e11','g2.79e12','g3.23e11','g3.49e11','g3.59e11','g3.61e11','g5.31e11','g5.36e11','g5.38e11','g5.55e11','g7.08e11','g7.44e11','g8.26e11','g8.28e11']
galaxies = ['g2.39e11','g3.23e11','g3.49e11','g3.59e11','g3.61e11','g5.31e11','g5.36e11','g5.38e11','g5.55e11','g7.08e11','g7.44e11','g8.26e11','g8.28e11','g1.12e12','g1.92e12','g2.79e12']


# SKIRT wavelength grid is the same for all galaxies                                                                                                                            
wave_SKIRT = np.load('/scratch/ntf229/RT_fit/resources/SKIRT/'+galaxies[0]+'/maxLevel13/wavelengths250/numPhotons3e8/inc0/dust/dustFraction0.2/maxTemp8000/wave.npy') * 1e4 # convert from microns to A     
att_mask_SKIRT = (wave_SKIRT >= 912) & (wave_SKIRT <= 5e4)

plt.figure(figsize=(10,8))

for i in range(len(galaxies)):

    # Calculate SKIRT attenuation curve (requires dust and nodust spectra)
    dustSpec_SKIRT = np.load('/scratch/ntf229/RT_fit/resources/SKIRT/'+galaxies[i]+'/maxLevel13/wavelengths250/numPhotons3e8/inc0/dust/dustFraction0.2/maxTemp8000/spec.npy') / 3631 # convert from Jy to maggies
    nodustSpec_SKIRT = np.load('/scratch/ntf229/RT_fit/resources/SKIRT/'+galaxies[i]+'/maxLevel13/wavelengths250/numPhotons3e8/inc0/nodust/spec.npy') / 3631
    dustSpec_SKIRT = 22.5 - 2.5*np.log10(dustSpec_SKIRT * 1e9) # convert to Pogson magnitudes
    nodustSpec_SKIRT = 22.5 - 2.5*np.log10(nodustSpec_SKIRT * 1e9) # convert to Pogson magnitudes
    att_curve_SKIRT = dustSpec_SKIRT[att_mask_SKIRT] - nodustSpec_SKIRT[att_mask_SKIRT] # x-axis for this is wave_SKIRT[att_mask_SKIRT]

    plt.plot(wave_SKIRT[att_mask_SKIRT], att_curve_SKIRT,label=galaxies[i])

plt.xlabel('Wavelength ('+r'$\AA$'+')',fontsize=16)
plt.ylabel(r'$A_{\lambda}$',fontsize=16)
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize=12)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig(plotPath+'SKIRT_attenuation_curves.png',dpi=600)
