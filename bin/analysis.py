import numpy as np
import os
from astropy.io import fits
import argparse
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from timeit import default_timer as timer
import datetime

parser = argparse.ArgumentParser()
parser.add_argument("--ageSmooth") # if True, smooth ages based on number of star particles (makeTextFiles parameter)
parser.add_argument("--SF") # if True, star particles younger than 10 Myrs are assigned MAPPINGS-III SEDs (makeTextFiles parameter)
parser.add_argument("--tauClear") # clearing time in Myrs for MAPPINGS-III f_PDR calculations (only matters if SF=True) (makeTextFiles parameter)
parser.add_argument("--clumps") # if True, add subgrid clumpiness to gas near MAPPINGS-III particles (makeTextFiles parameter)
parser.add_argument("--numCells") # number of cells along one dimension for subgrid clumping (only matters if clumps=True) (makeTextFiles parameter)
parser.add_argument("--numClumps") # number of clumps for subgrid clumping (only matters if clumps=True) (makeTextFiles parameter)
parser.add_argument("--inc") # inclination angle (SKIRT parameter)
parser.add_argument("--az") # azimuth angle (SKIRT parameter)
parser.add_argument("--numPhotons") # number of photon packages (SKIRT parameter)
parser.add_argument("--pixels") # number of pixels (square) for image (SKIRT parameter)
parser.add_argument("--dustFraction") # dust to metal ratio (SKIRT parameter)
parser.add_argument("--maxTemp") # maximum temperature at which dust can form (SKIRT parameter)
parser.add_argument("--SSP") # simple stellar population model including IMF after underscore (SKIRT parameter)
args = parser.parse_args()

codePath = '/home/ntf229/sphRad/'
resultPath = '/scratch/ntf229/sphRad/' # store results here

# Directory structure stores important parameters
textPath = resultPath+'resources/NIHAO/TextFiles/'
SKIRTPath = resultPath+'resources/SKIRT/'
noDustSKIRTPath = resultPath+'resources/SKIRT/'
plotPath = resultPath+'resources/Plots/'
massPath = resultPath+'resources/NIHAO/GlobalProps/stellarMasses/'
if eval(args.ageSmooth):
    textPath += 'ageSmooth/'
    SKIRTPath += 'ageSmooth/'
    noDustSKIRTPath += 'ageSmooth/'
    plotPath += 'ageSmooth/'
else:
    textPath += 'noAgeSmooth/'
    SKIRTPath += 'noAgeSmooth/'
    noDustSKIRTPath += 'noAgeSmooth/'
    plotPath += 'noAgeSmooth/'
if eval(args.SF):
    textPath += 'SF/tauClear'+args.tauClear+'/'
    SKIRTPath += 'SF/tauClear'+args.tauClear+'/'
    noDustSKIRTPath += 'SF/tauClear'+args.tauClear+'/'
    plotPath += 'SF/tauClear'+args.tauClear+'/'
else:
    textPath += 'noSF/'
    SKIRTPath += 'noSF/'
    noDustSKIRTPath += 'noSF/'
    plotPath += 'noSF/'
if eval(args.clumps):
    textPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
    SKIRTPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
    noDustSKIRTPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
    plotPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
else:
    textPath += 'noClumps/'
    SKIRTPath += 'noClumps/'
    noDustSKIRTPath += 'noClumps/'
    plotPath += 'noClumps/'

noDustPlotPath = plotPath

SKIRTPath += 'dust/dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'
plotPath += 'dust/dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'
noDustSKIRTPath += 'noDust/'
noDustPlotPath += 'noDust/'

SKIRTPath += 'numPhotons'+args.numPhotons+'/inc'+args.inc+'/az'+args.az+'/'+args.SSP+'/'
plotPath += 'numPhotons'+args.numPhotons+'/inc'+args.inc+'/az'+args.az+'/'+args.SSP+'/'
noDustSKIRTPath += 'numPhotons'+args.numPhotons+'/inc'+args.inc+'/az'+args.az+'/'+args.SSP+'/'
noDustPlotPath += 'numPhotons'+args.numPhotons+'/inc'+args.inc+'/az'+args.az+'/'+args.SSP+'/'

# SKIRT broadband photometry names 
band_names = ['FUV', 'NUV', 'u', 'g', 'r', 'i', 'z', '2MASS_J', '2MASS_H', '2MASS_KS', 'W1', 
              'W2', 'W3', 'W4', 'PACS70', 'PACS100', 'PACS160', 'SPIRE250', 'SPIRE350', 'SPIRE500']

# DustPedia data
dp_file = np.loadtxt(codePath+'resources/DustPedia_Aperture_Photometry_2.2.csv',dtype=str)
dp_file2 = np.loadtxt(codePath+'resources/DustPedia_HyperLEDA_Herschel.csv',dtype=str) 
dp_file3 = np.loadtxt(codePath+'resources/dustpedia_cigale_results_final_version.csv',dtype=str)
dp_num = len(dp_file) - 1 # first line is header
dp_flux_index = np.asarray([7,10,13,16,19,22,25,28,31,34,37,40,43,46,70,73,76,79,82,85]) # DP indicies corresponding to band_names
dp_err_index = dp_flux_index + 1 # index of DustPedia flux errors corresponding to band_names

# color-color plot info
#ratio_indicies = [[0,6,19,6], [0,6,15,6], [0,6,13,6], [0,6,11,6]]
#plot_names = ["spire500","pacs100","wise4","wise2"]
#xlims = [[2e-4, 0.6], [2e-4,1], [2e-4,1], [2e-4,1]]
#ylims = [[1e-1, 50], [1e-2,5e2], [6e-3,7e1], [5e-2,5]]
ratio_indicies = [[0,6,7,6], [0,6,8,6], [0,6,9,6], [0,6,10,6], [0,6,11,6], [0,6,12,6], [0,6,13,6], 
                  [0,6,14,6], [0,6,15,6], [0,6,16,6], [0,6,17,6], [0,6,18,6], [0,6,19,6]]
plot_names = ["J","H","K","W1","W2","W3","W4","pacs70","pacs100","pacs160","spire250","spire350","spire500"]

# axis ratio color plot info
#axisRatio_color_indices = [[0,6], [1,6], [2,6], [3,6], [4,6], [5,6]] # relative to z
#axisRatio_color_indices = [[0,9], [1,9], [2,9], [3,9], [4,9], [5,9]] # relative to k 
#axisRatio_color_indices = [[0,13], [1,13], [2,13], [3,13], [4,13], [5,13]] # relative to w4
axisRatio_color_indices = [[0,7], [1,7], [2,7], [3,7], [4,7], [5,7]] # relative to J
axisRatio_plot_names = ["FUV", "NUV", "u", "g", "r", "i"]

# Initialize DustPedia arrays
dp_flux = np.zeros((dp_num, len(band_names))) # first index specifies galaxy, second index specifies band 
dp_err = np.zeros((dp_num, len(band_names)))
dp_bool = np.zeros((dp_num, len(band_names))) # 1 indicates a negative flux (to be colored red) 
#dp_names = []
dp_axisRatio = np.zeros(dp_num)
dp_disk = np.zeros(dp_num) # 0 means not disk, 1 means disk (t >= 0)
dp_stellarMass = np.zeros(dp_num)
dp_dustMass = np.zeros(dp_num)

dp3_names = [] 
for i in range(dp_num):
    g3 = i+1 # skip first line of file (headers)
    dp3_names.append(dp_file3[g3].split(',')[0])
dp3_names = np.asarray(dp3_names) # doesn't include first line of dp3_file

start = timer()

# Fill DustPedia arrays
for i in range(dp_num):
    g = i+1 # index of current galaxy
    params = dp_file[g].split(',')
    params2 = dp_file2[g].split(',')
    g3 = np.where(dp3_names == params[0])[0][0] + 1 
    params3 = dp_file3[g3].split(',')
    t = float(params2[3])
    if t >= 0:
        dp_disk[i] = 1
    dp_axisRatio[i] = 1./float(params[4]) # change from a/b to b/a 
    dp_stellarMass[i] = params3[3] # in solar masses (CIGALE)
    dp_dustMass[i] = params3[17] # in solar masses (CIGALE)
    for j in range(len(band_names)):
        if params[int(dp_flux_index[j])]:
            dp_flux[i,j] = float(params[int(dp_flux_index[j])])
            dp_err[i,j] = float(params[int(dp_err_index[j])])
            if dp_flux[i,j] - 2*dp_err[i,j] <= 0:
                dp_flux[i,j] = 2 * dp_err[i,j] # if flux is consistent with 0, set to 2*sigma 
                dp_bool[i,j] = 1 # signals flux is consistent with 0 (upper limit)
        else:
            dp_flux[i,j] = float("NaN") # set flux and error to NaN if photometry not available 
            dp_err[i,j] = float("NaN")

diskMask = dp_disk == 1
faceOnMask = dp_axisRatio[diskMask] > 0.85 

end = timer()
time = end - start
time = str(datetime.timedelta(seconds=time))
print('Time to fill dp arrays:', time)  

# Final sample (66)
#galaxies = ['g1.88e10','g1.89e10','g1.90e10','g2.34e10','g2.63e10','g2.64e10','g2.80e10','g2.83e10','g2.94e10',
#            'g3.19e10','g3.44e10','g3.67e10','g3.93e10','g4.27e10','g4.48e10','g4.86e10','g4.94e10','g4.99e10',
#            'g5.05e10','g6.12e10','g6.37e10','g6.77e10','g6.91e10','g7.12e10','g8.89e10','g9.59e10','g1.05e11',
#            'g1.08e11','g1.37e11','g1.52e11','g1.57e11','g1.59e11','g1.64e11','g2.04e11','g2.19e11','g2.39e11',
#            'g2.41e11','g2.42e11','g2.54e11','g3.06e11','g3.21e11','g3.23e11','g3.49e11','g3.55e11','g3.59e11',
#            'g3.71e11','g4.90e11','g5.02e11','g5.31e11','g5.36e11','g5.38e11','g5.46e11','g5.55e11','g6.96e11',
#            'g7.08e11','g7.44e11','g7.55e11','g7.66e11','g8.06e11','g8.13e11','g8.26e11','g8.28e11','g1.12e12',
#            'g1.77e12','g1.92e12','g2.79e12']

# Final sample (65) excluding g3.19e10 (conatins two galaxies)
galaxies = ['g1.88e10','g1.89e10','g1.90e10','g2.34e10','g2.63e10','g2.64e10','g2.80e10','g2.83e10','g2.94e10',
            'g3.44e10','g3.67e10','g3.93e10','g4.27e10','g4.48e10','g4.86e10','g4.94e10','g4.99e10',
            'g5.05e10','g6.12e10','g6.37e10','g6.77e10','g6.91e10','g7.12e10','g8.89e10','g9.59e10','g1.05e11',
            'g1.08e11','g1.37e11','g1.52e11','g1.57e11','g1.59e11','g1.64e11','g2.04e11','g2.19e11','g2.39e11',
            'g2.41e11','g2.42e11','g2.54e11','g3.06e11','g3.21e11','g3.23e11','g3.49e11','g3.55e11','g3.59e11',
            'g3.71e11','g4.90e11','g5.02e11','g5.31e11','g5.36e11','g5.38e11','g5.46e11','g5.55e11','g6.96e11',
            'g7.08e11','g7.44e11','g7.55e11','g7.66e11','g8.06e11','g8.13e11','g8.26e11','g8.28e11','g1.12e12',
            'g1.77e12','g1.92e12','g2.79e12']

# Initialize arrays for SKIRT spatially integrated photometry and axis ratios for color plots
SKIRT_flux = np.zeros((len(galaxies), len(band_names)))
SKIRT_axisRatio = np.zeros(len(galaxies))
SKIRT_stellarMass = np.zeros(len(galaxies))
SKIRT_dustMass = np.zeros(len(galaxies))

start = timer()

os.system('mkdir -p '+plotPath+'SEDs/')
os.system('mkdir -p '+noDustPlotPath+'SEDs/')
os.system('mkdir -p '+plotPath+'zrgImages/')
os.system('mkdir -p '+plotPath+'AttenuationCurves/')
#os.system('rm -rf '+plotPath+'autoprof/') # remove directory if it exists
for i in range(len(galaxies)):
    print('starting '+galaxies[i])
    # stellar mass
    if os.path.isfile(massPath+galaxies[i]+'/stellarMass.npy'):
        SKIRT_stellarMass[i] = float(np.load(massPath+galaxies[i]+'/stellarMass.npy'))
    else:
        stars = np.loadtxt(textPath+galaxies[i]+'/stars.txt')
        youngStars = np.loadtxt(textPath+galaxies[i]+'/youngStars.txt')
        SKIRT_stellarMass[i] = np.sum(stars[:,7]) + np.sum(youngStars[:,7] * 1.e7) # in Msun
        os.system('mkdir -p '+massPath+galaxies[i]+'/')
        np.save(massPath+galaxies[i]+'/stellarMass.npy', SKIRT_stellarMass[i])
    # mass in metals
    if os.path.isfile(massPath+galaxies[i]+'/maxTemp'+args.maxTemp+'/metalMass.npy'):
        SKIRT_dustMass[i] = float(np.load(massPath+galaxies[i]+'/maxTemp'+args.maxTemp+'/metalMass.npy')) * float(args.dustFraction)
    else:
        gas = np.loadtxt(textPath+galaxies[i]+'/gas.txt')
        tempMask = gas[:,6] < float(args.maxTemp)
        metalMass = np.sum(gas[tempMask, 4] * gas[tempMask, 5]) # in Msun
        SKIRT_dustMass[i] = metalMass * float(args.dustFraction)
        os.system('mkdir -p '+massPath+galaxies[i]+'/maxTemp'+args.maxTemp+'/')
        np.save(massPath+galaxies[i]+'/maxTemp'+args.maxTemp+'/metalMass.npy', metalMass)
    bb = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_broadband_sed.dat', unpack = True)
    bb_wave = bb[0] # spatially integrated broadband wavelengths in microns
    bb_wave = bb_wave * 1e4 # convert to Angstroms
    bb_spec = bb[1] # spatially integrated broadband fluxes in Janskys
    SKIRT_flux[i,:] = bb_spec
    # skip this part of analysis if already done
    if os.path.isfile(plotPath+'autoprof/'+galaxies[i]+'/'+galaxies[i]+'-r-ba.csv'):
        print('skipping first analysis')
        axisRatioFile = np.loadtxt(plotPath+'autoprof/'+galaxies[i]+'/'+galaxies[i]+'-r-ba.csv', dtype=str)
        SKIRT_axisRatio[i] = float(axisRatioFile[1,1].split(',')[0])
        continue
    # SEDs 
    print('starting SEDs')
    sed = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_SED_sed.dat', unpack = True)
    wave = sed[0] # spatially integrated SED wavelengths in microns
    wave = wave * 1e4 # convert to Angstroms
    spec = sed[1] # spatially integrated SED fluxes in Janskys
    plt.figure(figsize=(10,8))
    mask = (wave >= 1e3) & (wave <= 1e7)
    plt.plot(wave[mask], spec[mask], color='k',alpha=1)
    plt.xlabel(r'$\lambda$'+' ['+r'$\AA$'+']', fontsize=16)
    plt.ylabel(r'$F_{\nu}$'+' [Jy]', fontsize=16)
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(plotPath+'SEDs/'+galaxies[i]+'_SED.png', dpi=300)
    plt.close()
    # noDust SEDs
    noDustSed = np.loadtxt(noDustSKIRTPath+galaxies[i]+'/sph_SED_sed.dat', unpack = True)
    noDustWave = noDustSed[0] # spatially integrated SED wavelengths in microns
    noDustWave = noDustWave * 1e4 # convert to Angstroms
    noDustSpec = noDustSed[1] # spatially integrated SED fluxes in Janskys
    plt.figure(figsize=(10,8))
    noDustMask = (noDustWave >= 1e3) & (noDustWave <= 1e7)
    plt.plot(noDustWave[noDustMask], noDustSpec[noDustMask], color='k', alpha=1)
    plt.xlabel(r'$\lambda$'+' ['+r'$\AA$'+']', fontsize=16)
    plt.ylabel(r'$F_{\nu}$'+' [Jy]', fontsize=16)
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(noDustPlotPath+'SEDs/'+galaxies[i]+'_SED.png', dpi=300)
    plt.close()
    # Attenuation Curves
    att_mask = (wave >= 912) & (wave <= 5e4)
    dustMags = 22.5 - 2.5*np.log10((spec/3631) * 1e9) # convert to Pogson magnitudes
    noDustMags = 22.5 - 2.5*np.log10((noDustSpec/3631) * 1e9) # convert to Pogson magnitudes
    plt.figure(figsize=(10,8))
    plt.plot(wave[att_mask], dustMags[att_mask] - noDustMags[att_mask], color='k', alpha=1)
    plt.xlabel(r'$\lambda$'+' ['+r'$\AA$'+']', fontsize=16)
    plt.ylabel(r'$A_{\lambda}$',fontsize=16)
    plt.xscale('log')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(plotPath+'AttenuationCurves/'+galaxies[i]+'.png', dpi=300)
    plt.close()
    # Optical images
    print('starting images')
    os.system('mkdir -p '+plotPath+'zrgImages/'+galaxies[i]+'/')
    cube_file = fits.open(SKIRTPath+galaxies[i]+'/sph_broadband_total.fits')
    cube = np.asarray(cube_file[0].data) # (2000, 2000, 20) cube of broadbands in MJy/sr 
    r_grid = np.asarray(cube[6,:,:]) # sdss_z band (see band_names for indicies)
    g_grid = np.asarray(cube[4,:,:]) # sdss_r band
    b_grid = np.asarray(cube[3,:,:]) # sdss_g band
    fig = plt.figure()
    sizes = np.shape(r_grid)
    fig.set_size_inches(1. * sizes[0] / sizes[1], 1, forward = False)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    tot_max = np.amax([np.amax(r_grid),np.amax(g_grid),np.amax(b_grid)])
    r_grid = 255 * r_grid / tot_max
    g_grid = 255 * g_grid / tot_max
    b_grid = 255 * b_grid / tot_max
    image = make_lupton_rgb(r_grid, g_grid, b_grid, stretch=3, Q=8) 
    ax.imshow(image,interpolation='none')
    plt.savefig(plotPath+'zrgImages/'+galaxies[i]+'_zrg.png',dpi=sizes[0])
    plt.close()
    # Axis Ratios
    print('starting autoprof')
    os.system('mkdir -p '+plotPath+'autoprof/'+galaxies[i]+'/')
    os.system('python '+codePath+'python/mkprof.py --galaxy='+galaxies[i]+
              ' --band=r --path='+SKIRTPath+galaxies[i]+' --save='+plotPath+'autoprof/'+galaxies[i]+'/')
    axisRatioFile = np.loadtxt(plotPath+'autoprof/'+galaxies[i]+'/'+galaxies[i]+'-r-ba.csv', dtype=str)
    SKIRT_axisRatio[i] = float(axisRatioFile[1,1].split(',')[0])
    os.system('rm '+plotPath+'autoprof/'+galaxies[i]+'/galaxy-r.fits')

end = timer()
time = end - start
time = str(datetime.timedelta(seconds=time))
print('Time for first analysis:', time)  

start = timer()

# Color-color Plots
print('starting color color plots')
os.system('mkdir -p '+plotPath+'colorPlots/')
for i in range(len(ratio_indicies)): # number of plots to make 
    if os.path.isfile(plotPath+'colorPlots/'+plot_names[i]+'.png'): # skip if already done
        print('color color plots already done')
        continue
    plt.figure(figsize=(10,8))
    # DustPedia propogation of errors
    xerr = np.absolute(dp_flux[:,ratio_indicies[i][0]] / dp_flux[:,ratio_indicies[i][1]]) * np.sqrt( 
    			      (dp_err[:,ratio_indicies[i][0]]/dp_flux[:,ratio_indicies[i][0]])**2 + 
    			      (dp_err[:,ratio_indicies[i][1]]/dp_flux[:,ratio_indicies[i][1]])**2 )
    yerr = np.absolute(dp_flux[:,ratio_indicies[i][2]] / dp_flux[:,ratio_indicies[i][3]]) * np.sqrt( 
    			      (dp_err[:,ratio_indicies[i][2]]/dp_flux[:,ratio_indicies[i][2]])**2 + 
    			      (dp_err[:,ratio_indicies[i][3]]/dp_flux[:,ratio_indicies[i][3]])**2 )
    colors = np.empty(dp_num,dtype=str)
    for d in range(dp_num):
        bad = False
        xuplims = False
        xlolims = False
        uplims = False
        lolims = False
        colors[d] = 'k'
        if dp_bool[d,ratio_indicies[i][0]] and dp_bool[d,ratio_indicies[i][1]]: # bad data point (0/0)
            bad = True # don't plot this point
        elif dp_bool[d,ratio_indicies[i][0]]: # x is an upper limit (numerator is upper limit)
            xuplims = True
            colors[d] = 'red'
        elif dp_bool[d,ratio_indicies[i][1]]: # x is a lower limit (denominator is upper limit)
            xlolims = True
            colors[d] = 'red'
        if dp_bool[d,ratio_indicies[i][2]] and dp_bool[d,ratio_indicies[i][3]]: # bad data point (0/0)
            bad = True # don't plot this point
        elif dp_bool[d,ratio_indicies[i][2]]: # y is an upper limit (numerator is upper limit)
            uplims = True
            colors[d] = 'red'
        elif dp_bool[d,ratio_indicies[i][3]]: # y is a lower limit (denominator is upper limit)
            lolims = True
            colors[d] = 'red'
        if bad:
            print('bad DustPedia data point')
        else:
            plt.errorbar(dp_flux[d,ratio_indicies[i][0]] / dp_flux[d,ratio_indicies[i][1]], 
                         dp_flux[d,ratio_indicies[i][2]] / dp_flux[d,ratio_indicies[i][3]], 
                         xerr=xerr[d], yerr=yerr[d], elinewidth=0.2, marker='o',
                         markersize=5, linewidth=0, color=colors[d], zorder=0, alpha=0.3,
                         xuplims=xuplims, xlolims=xlolims, uplims=uplims, lolims=lolims)
    plt.scatter(SKIRT_flux[:,ratio_indicies[i][0]] / SKIRT_flux[:,ratio_indicies[i][1]], 
                SKIRT_flux[:,ratio_indicies[i][2]] / SKIRT_flux[:,ratio_indicies[i][3]], 
                marker='o', s=20, zorder=10,alpha=0.7, c='blue')
    plt.xlabel(band_names[ratio_indicies[i][0]]+' / '+band_names[ratio_indicies[i][1]], fontsize=16)
    plt.ylabel(band_names[ratio_indicies[i][2]]+' / '+band_names[ratio_indicies[i][3]], fontsize=16)
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.xlim((xlims[i][0], xlims[i][1]))
    #plt.ylim((ylims[i][0], ylims[i][1]))
    plt.savefig(plotPath+'colorPlots/'+plot_names[i]+'.png',dpi=300)
    plt.close()

end = timer()
time = end - start
time = str(datetime.timedelta(seconds=time))
print('Time for color color plots:', time)  

# axis ratio vs. color plots
print('starting axis ratio color plots')
os.system('mkdir -p '+plotPath+'fitPlots/')
os.system('mkdir -p '+plotPath+'axisRatioColorPlots/')
for i in range(len(axisRatio_color_indices)):
    # Calculate average DustPedia face-on disk galaxy colors as function of stellar mass
    # polyfit does not handle NaNs, need to mask them out first
    dpFaceOnDustToStellar = np.log10(dp_dustMass[diskMask][faceOnMask] / dp_stellarMass[diskMask][faceOnMask])
    dpFaceOnFluxRatios = np.log10((dp_flux[diskMask][faceOnMask][:, axisRatio_color_indices[i][0]] / 
                         dp_flux[diskMask][faceOnMask][:, axisRatio_color_indices[i][1]]))
    nanMask = np.isfinite(dpFaceOnDustToStellar) & np.isfinite(dpFaceOnFluxRatios)
    linearFit = np.polyfit(dpFaceOnDustToStellar[nanMask], dpFaceOnFluxRatios[nanMask], 1)
    fit = np.poly1d(linearFit) # use fit(dp_stellarMass) to get average color
    # make plots showing fits
    plt.figure(figsize=(10,8))
    plt.scatter(dpFaceOnDustToStellar[nanMask], dpFaceOnFluxRatios[nanMask], marker='o', s=15)
    x_fit = np.linspace(np.amin(dpFaceOnDustToStellar[nanMask]), np.amax(dpFaceOnDustToStellar[nanMask]), num=2)
    y_fit = fit(x_fit)
    plt.plot(x_fit, y_fit, color='k')
    plt.xlabel('log(Dust Mass / Stellar Mass)', fontsize=16)
    plt.ylabel('log('+band_names[axisRatio_color_indices[i][0]]+' / '+
               band_names[axisRatio_color_indices[i][1]]+')', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(plotPath+'fitPlots/'+axisRatio_plot_names[i]+'.png',dpi=300)
    plt.close()
    # calculate errors
    x = dp_flux[:,axisRatio_color_indices[i][0]]
    y = dp_flux[:,axisRatio_color_indices[i][1]]
    x_err = dp_err[:,axisRatio_color_indices[i][0]]
    y_err = dp_err[:,axisRatio_color_indices[i][1]]
    yerr = np.sqrt((x_err**2 * (1/(x * np.log(10))**2)) + (y_err**2 * (1/(y * np.log(10)))**2))
    colors = np.empty(dp_num, dtype=str)
    plt.figure(figsize=(10,8))
    plt.scatter(SKIRT_axisRatio, 
                np.log10(SKIRT_flux[:,axisRatio_color_indices[i][0]] / 
                SKIRT_flux[:,axisRatio_color_indices[i][1]]) - 
                fit(np.log10(SKIRT_dustMass / SKIRT_stellarMass)), 
                marker='o', s=20, zorder=10, alpha=0.7, c='blue')
    for d in range(dp_num):
        if not diskMask[d]: # skip non-disk galaxies
            continue 
        bad = False
        xuplims = False
        xlolims = False
        uplims = False
        lolims = False
        colors[d] = 'k'
        if dp_bool[d,axisRatio_color_indices[i][0]] and dp_bool[d,axisRatio_color_indices[i][1]]: # bad data point (0/0)
            bad = True # don't plot this point
        elif dp_bool[d,axisRatio_color_indices[i][0]]: # y is an upper limit (numerator is upper limit)
            uplims = True
            colors[d] = 'red'
        elif dp_bool[d,axisRatio_color_indices[i][1]]: # y is a lower limit (denominator is upper limit)
            lolims = True
            colors[d] = 'red'
        if bad:
            print('bad DustPedia data point')
        else:
            plt.errorbar(dp_axisRatio[d], 
                         np.log10(dp_flux[d,axisRatio_color_indices[i][0]] / 
                         dp_flux[d,axisRatio_color_indices[i][1]]) - 
                         fit(np.log10(dp_dustMass[d] / dp_stellarMass[d])), 
                         xerr=0, yerr=yerr[d], elinewidth=0.2, marker='o',
                         markersize=5, linewidth=0, color=colors[d], zorder=0, alpha=0.3,
                         xuplims=xuplims, xlolims=xlolims, uplims=uplims, lolims=lolims)
    plt.xlabel('Axis Ratio', fontsize=16)
    plt.ylabel('log('+band_names[axisRatio_color_indices[i][0]]+' / '+band_names[axisRatio_color_indices[i][1]]+') - <log('+
               band_names[axisRatio_color_indices[i][0]]+' / '+band_names[axisRatio_color_indices[i][1]]+')>', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(plotPath+'axisRatioColorPlots/'+axisRatio_plot_names[i]+'.png',dpi=300)
    plt.close()
print('done')
