# Make color-color plots for SKIRT generated spectra 

import numpy as np
import matplotlib.pyplot as plt
from sedpy import observate
import os

dust = True # True means use SKIRT with dust spectra 
partialDust = False

file = np.loadtxt('/home/ntf229/DustPedia/data/DustPedia_Aperture_Photometry_2.2.csv',dtype=str)
if dust:
    plot_path = '/scratch/ntf229/RT_fit/plots/color_plots/'
elif partialDust:
    plot_path = '/scratch/ntf229/RT_fit/plots/color_plots_partialDust/'
else:
    plot_path = '/scratch/ntf229/RT_fit/plots/color_plots_nodust/'

os.system('mkdir -p '+plot_path)

# file[0] are column names
# columns are separated by commas
# empty columns are still separated by a comma
# All the photometric measurements are in Jy

dp_num = len(file)
print('Number of galaxies in sample:', dp_num)

#bands = ["galex_FUV", "galex_NUV", "sdss_z0", "wise_w1", "wise_w3", "wise_w4", "herschel_pacs_100", "herschel_spire_250", "herschel_spire_500"]

#bands = ["galex_FUV", "galex_NUV", "sdss_u0", "sdss_r0", "sdss_z0", "wise_w1", "wise_w3", "wise_w4", "herschel_pacs_100", "herschel_spire_250", "herschel_spire_500"]
#dp_flux_index = np.asarray([7,10,13,19,25,37,43,46,73,79,85]) # index of DustPedia flux corresponding to the above bands

bands = ["galex_FUV", "galex_NUV", "sdss_u0", "sdss_g0", "sdss_r0", "sdss_i0", 
        "sdss_z0", "wise_w1", "wise_w2", "wise_w3", "wise_w4", "spitzer_irac_ch1", 
        "spitzer_irac_ch2", "spitzer_irac_ch3", "spitzer_irac_ch4", "spitzer_mips_24", 
        "spitzer_mips_70", "spitzer_mips_160", "herschel_pacs_100", "herschel_spire_250", 
        "herschel_spire_500"]
dp_flux_index = np.asarray([7,10,13,16,19,22,25,37,40,43,46,49,52,55,58,61,64,67,73,79,85]) # index of DustPedia flux corresponding to the above bands

dp_err_index = dp_flux_index + 1 # index of DustPedia flux errors corresponding to the above bands

# ratio_indicies = [[1,2,3,4],[5,6,7,8]] will give 2 plots, first plot with x-axis 1/2 and y-axis 3/4, where the numbers indicate the above bands
#ratio_indicies = [[4,8,6,8], [1,0,9,0]] # two plots: (1, 10, 100 microns) and (NUV/FUV vs. SPIRE_250/FUV)
#plot_names = ["1_10_100_microns", "UV_dust_emission"]
#ratio_indicies = [[4,8,6,8], [1,0,9,0], [0,4,1,4], [2,4,3,4], [5,4,6,4], [7,4,8,4], [9,4,10,4], [0,4,8,4], [1,4,8,4]]
#plot_names = ["1_10_100_microns", "UV_dust_emission","FUV_NUV","r_u","w1_w3","w4_pacs100","spire250_500","FUV_FIR","NUV_FIR"]
ratio_indicies = [[0,6,7,6], [0,6,8,6], [0,6,9,6], [0,6,10,6], [0,6,11,6], [0,6,12,6], 
                  [0,6,13,6], [0,6,14,6], [0,6,15,6], [0,6,16,6], [0,6,17,6], [0,6,18,6], 
                  [0,6,19,6], [0,6,20,6], [8,6,18,6]]
plot_names = ["w1","w2","w3","w4","spitzer3.6","spitzer4.5","spitzer5.8","spitzer8",
              "spitzer24","spitzer70","spitzer160","pacs100","spire250","spire500","w2_pacs100"]

dp_flux = np.zeros((dp_num, len(bands))) # first index specifies galaxy, second index specifies band 
dp_err = np.zeros((dp_num, len(bands)))

dp_bool = np.zeros((dp_num, len(bands))) # 1 indicates a negative flux (to be colored red) 

for i in range(len(file)-1):
    g = i+1 # index of current galaxy
    params = file[g].split(',')

    for j in range(len(bands)):
        if params[int(dp_flux_index[j])]:
            dp_flux[i,j] = float(params[int(dp_flux_index[j])])
            dp_err[i,j] = float(params[int(dp_err_index[j])])
            if dp_flux[i,j] <= 0:
                dp_flux[i,j] = 2 * dp_err[i,j] # if flux is negative, set to 2*sigma 
                dp_bool[i,j] = 1 # signals a negative flux for this galaxy/band
            elif dp_flux[i,j] - dp_err[i,j] <= 0: # if flux is consistent with 0 
                dp_bool[i,j] = 1
        else:
            dp_flux[i,j] = float("NaN") # set flux and error to NaN if photometry not available 
            dp_err[i,j] = float("NaN")

# Final sample (66)
galaxies = ['g1.88e10','g1.89e10','g1.90e10','g2.34e10','g2.63e10','g2.64e10','g2.80e10','g2.83e10','g2.94e10',
            'g3.19e10','g3.44e10','g3.67e10','g3.93e10','g4.27e10','g4.48e10','g4.86e10','g4.94e10','g4.99e10',
            'g5.05e10','g6.12e10','g6.37e10','g6.77e10','g6.91e10','g7.12e10','g8.89e10','g9.59e10','g1.05e11',
            'g1.08e11','g1.37e11','g1.52e11','g1.57e11','g1.59e11','g1.64e11','g2.04e11','g2.19e11','g2.39e11',
            'g2.41e11','g2.42e11','g2.54e11','g3.06e11','g3.21e11','g3.23e11','g3.49e11','g3.55e11','g3.59e11',
            'g3.71e11','g4.90e11','g5.02e11','g5.31e11','g5.36e11','g5.38e11','g5.46e11','g5.55e11','g6.96e11',
            'g7.08e11','g7.44e11','g7.55e11','g7.66e11','g8.06e11','g8.13e11','g8.26e11','g8.28e11','g1.12e12',
            'g1.77e12','g1.92e12','g2.79e12']

dustFraction = ['0.2'] 
#dustFraction = ['0.2']
#dustFraction = ['0.05','0.1','0.2','0.3','0.4']
#dustFraction = ['0.2','0.4','0.6','0.8','1']

SKIRT_flux = np.zeros((len(galaxies), len(dustFraction), len(bands)))

if dust:
    SKIRT_path1 = '/scratch/ntf229/RT_fit/resources/SKIRT/'
    SKIRT_path2 = '/numPhotons1e9/inc0/dust/dustFraction'
    #SKIRT_path2 = '/maxLevel13/numPhotons1e9/inc0/dust/dustFraction'
    SKIRT_path3 = '/maxTemp16000/'
elif partialDust:
    SKIRT_path1 = '/scratch/ntf229/RT_fit/resources/SKIRT_partialDust/'
    SKIRT_path2 = '/maxLevel13/numPhotons1e9/inc0/dust/dustFraction'
    SKIRT_path3 = '/maxTemp16000/'
else:
    SKIRT_path1 = '/scratch/ntf229/RT_fit/resources/SKIRT/'
    SKIRT_path2 = '/maxLevel13/numPhotons1e9/inc0/nodust/'

avg_gas_metals = np.zeros(len(galaxies))
dust_mass_fraction = np.zeros(len(galaxies))
totStarMass = np.zeros(len(galaxies))
numGasPart = np.zeros(len(galaxies))

for i in range(len(galaxies)):
    for j in range(len(dustFraction)):
        if dust or partialDust:
            wave = np.load(SKIRT_path1+galaxies[i]+SKIRT_path2+dustFraction[j]+SKIRT_path3+'wave.npy') * 1e4 # convert to Angstroms
            spec = np.load(SKIRT_path1+galaxies[i]+SKIRT_path2+dustFraction[j]+SKIRT_path3+'spec.npy') # in Jy
        else:
            wave = np.load(SKIRT_path1+galaxies[i]+SKIRT_path2+'wave.npy') * 1e4 # convert to Angstroms
            spec = np.load(SKIRT_path1+galaxies[i]+SKIRT_path2+'spec.npy') # in Jy
        filterlist = observate.load_filters(bands)
        f_lambda_cgs = (1/33333) * (1/(wave**2)) * spec 
        mags = observate.getSED(wave, f_lambda_cgs, filterlist=filterlist) # in AB Magnitudes
        jy = 1e26 * 10**(-1*(mags+48.6)/2.5) # convert from AB mags to Janskys (same order as bands)
        for k in range(len(bands)):
            SKIRT_flux[i,j,k] = jy[k]
        # load gas text files
        textFilePath = '/scratch/ntf229/RT_fit/resources/NIHAO/TextFiles/'+galaxies[i]+'/'
        gas = np.loadtxt(textFilePath+'gas.txt')
        stars = np.loadtxt(textFilePath+'stars.txt')
        #youngStars = np.loadtxt(textFilePath+'youngStars.txt')
        starMasses = stars[:,7] # units of M_sun
        #youngStarMasses = youngStars[:,7] * 1e7
        totStarMass[i] = np.sum(starMasses)
        #gasDensity = np.loadtxt(textFilePath+'gas_density.txt')
        gasMasses = gas[:,4] # units of M_sun
        #trueMass = np.sum(gasMasses)
        gasTemp = gas[:,6] # units of K
        gasMetals = gas[:,5]
        numGasPart[i] = len(gasMasses)
        #avg_gas_metals[i] = np.average(gasMetals, weights=gasMasses)
        #dust_mass_fraction[i] = 0.2 * np.sum(gasMetals*gasMasses) / totStarMass[i] # dust mass divided by total stellar mass (no temp cuts)
        for k in range(len(gasMasses)):
            if gasTemp[k] <= 16000: # including temperature cut 
                dust_mass_fraction[i] += gasMetals[k]*gasMasses[k]
        dust_mass_fraction[i] = float(dustFraction[j]) * dust_mass_fraction[i] / totStarMass[i] # dust mass divided by total stellar mass (with temp cuts)
        #dust_mass_fraction[i] =  float(dustFraction[j]) * dust_mass_fraction[i] # dust mass (with temp cuts)
        


# Make plots
for i in range(len(ratio_indicies)): # number of plots to make 
    plt.figure(figsize=(10,8))
    
    # Propogation of errors
    xerr = np.absolute(dp_flux[:,ratio_indicies[i][0]] / dp_flux[:,ratio_indicies[i][1]]) * np.sqrt( 
    			(dp_err[:,ratio_indicies[i][0]]/dp_flux[:,ratio_indicies[i][0]])**2 + (dp_err[:,ratio_indicies[i][1]]/dp_flux[:,ratio_indicies[i][1]])**2 )
    yerr = np.absolute(dp_flux[:,ratio_indicies[i][2]] / dp_flux[:,ratio_indicies[i][3]]) * np.sqrt( 
    			(dp_err[:,ratio_indicies[i][2]]/dp_flux[:,ratio_indicies[i][2]])**2 + (dp_err[:,ratio_indicies[i][3]]/dp_flux[:,ratio_indicies[i][3]])**2 )
    
    colors = np.empty(dp_num,dtype=str)
    for d in range(dp_num):
        if dp_bool[d,ratio_indicies[i][0]] or dp_bool[d,ratio_indicies[i][1]] or dp_bool[d,ratio_indicies[i][2]] or dp_bool[d,ratio_indicies[i][3]]:
            colors[d] = 'red'
        else:
            colors[d] = 'k'
        plt.errorbar(dp_flux[d,ratio_indicies[i][0]] / dp_flux[d,ratio_indicies[i][1]], dp_flux[d,ratio_indicies[i][2]] / dp_flux[d,ratio_indicies[i][3]], 
                    xerr=xerr[d], yerr=yerr[d], elinewidth=0.2, marker='o',
                    markersize=2, linewidth=0, color=colors[d], zorder=0, alpha=0.3)
    
    #plt.errorbar(dp_flux[:,ratio_indicies[i][0]] / dp_flux[:,ratio_indicies[i][1]], dp_flux[:,ratio_indicies[i][2]] / dp_flux[:,ratio_indicies[i][3]], 
    #			 xerr=xerr, yerr=yerr, elinewidth=0.2, marker='o', 
    #			 markersize=2, linewidth=0, color=colors, label='DustPedia',zorder=0,alpha=0.3)
    for j in range(len(dustFraction)):
        #plt.errorbar(SKIRT_flux[:,j,ratio_indicies[i][0]] / SKIRT_flux[:,j,ratio_indicies[i][1]], 
        #            SKIRT_flux[:,j,ratio_indicies[i][2]] / SKIRT_flux[:,j,ratio_indicies[i][3]], marker='o', markersize=5, 
        #            linewidth=0, label='SKIRT dust-to-metal ratio '+dustFraction[j],zorder=10,alpha=0.7)
        plt.scatter(SKIRT_flux[:,j,ratio_indicies[i][0]] / SKIRT_flux[:,j,ratio_indicies[i][1]], 
                            SKIRT_flux[:,j,ratio_indicies[i][2]] / SKIRT_flux[:,j,ratio_indicies[i][3]], marker='o', s=8, 
                            label='SKIRT dust-to-metal ratio '+dustFraction[j],zorder=10,alpha=0.7,cmap='rainbow',c=np.log10(dust_mass_fraction))
        for k in range(len(galaxies)):
            plt.annotate(galaxies[k],(SKIRT_flux[k,j,ratio_indicies[i][0]] / SKIRT_flux[k,j,ratio_indicies[i][1]], 
                        SKIRT_flux[k,j,ratio_indicies[i][2]] / SKIRT_flux[k,j,ratio_indicies[i][3]]),fontsize=5)
    plt.xlabel(bands[ratio_indicies[i][0]]+' / '+bands[ratio_indicies[i][1]], fontsize=16)
    plt.ylabel(bands[ratio_indicies[i][2]]+' / '+bands[ratio_indicies[i][3]], fontsize=16)
    plt.xscale('log')
    plt.yscale('log')
    #plt.legend(fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.colorbar(label='Gas Metallicity (Mass Fraction)')
    plt.colorbar(label='Log(Dust Mass / Stellar Mass)')
    #plt.colorbar(label='Number of Gas Particles')
    #plt.xlim((0.08, 50))
    #plt.ylim((8e-4, 1))
    plt.savefig(plot_path+plot_names[i]+'.png',dpi=600)
    plt.close()


print('done')
