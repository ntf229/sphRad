import numpy as np
import os
from astropy.io import fits
import argparse
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb

parser = argparse.ArgumentParser()
parser.add_argument("--ageSmooth") # if True, smooth ages based on number of star particles (makeTextFiles parameter)
parser.add_argument("--SF") # if True, star particles younger than 10 Myrs are assigned MAPPINGS-III SEDs (makeTextFiles parameter)
parser.add_argument("--tauClear") # clearing time in Myrs for MAPPINGS-III f_PDR calculations (only matters if SF=True) (makeTextFiles parameter)
parser.add_argument("--clumps") # if True, add subgrid clumpiness to gas near MAPPINGS-III particles (makeTextFiles parameter)
parser.add_argument("--numCells") # number of cells along one dimension for subgrid clumping (only matters if clumps=True) (makeTextFiles parameter)
parser.add_argument("--numClumps") # number of clumps for subgrid clumping (only matters if clumps=True) (makeTextFiles parameter)
parser.add_argument("--dust") # include dust; True or False (RT parameter)
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
SKIRTPath = resultPath+'resources/SKIRT_views/'
plotPath = resultPath+'resources/Plots_views/'
if eval(args.ageSmooth):
    textPath += 'ageSmooth/'
    SKIRTPath += 'ageSmooth/'
    plotPath += 'ageSmooth/'
else:
    textPath += 'noAgeSmooth/'
    SKIRTPath += 'noAgeSmooth/'
    plotPath += 'noAgeSmooth/'
if eval(args.SF):
    textPath += 'SF/tauClear'+args.tauClear+'/'
    SKIRTPath += 'SF/tauClear'+args.tauClear+'/'
    plotPath += 'SF/tauClear'+args.tauClear+'/'
else:
    textPath += 'noSF/'
    SKIRTPath += 'noSF/'
    plotPath += 'noSF/'
if eval(args.clumps):
    textPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
    SKIRTPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
    plotPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
else:
    textPath += 'noClumps/'
    SKIRTPath += 'noClumps/'
    plotPath += 'noClumps/'
if eval(args.dust):
    SKIRTPath += 'dust/dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'
    plotPath += 'dust/dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'
else:
    SKIRTPath += 'noDust/'
    plotPath += 'noDust/'

SKIRTPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
plotPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'

# SKIRT broadband photometry names 
band_names = ['FUV', 'NUV', 'u', 'g', 'r', 'i', 'z', '2MASS_J', '2MASS_H', '2MASS_KS', 'W1', 
              'W2', 'W3', 'W4', 'PACS70', 'PACS100', 'PACS160', 'SPIRE250', 'SPIRE350', 'SPIRE500']

# DustPedia data
dp_file = np.loadtxt(codePath+'resources/DustPedia_Aperture_Photometry_2.2.csv',dtype=str)
dp_num = len(dp_file)
dp_flux_index = np.asarray([7,10,13,16,19,22,25,28,31,34,37,40,43,46,70,73,76,79,82,85]) # DP indicies corresponding to band_names
dp_err_index = dp_flux_index + 1 # index of DustPedia flux errors corresponding to band_names

#ratio_indicies = [[0,6,19,6], [0,6,15,6], [0,6,13,6], [0,6,12,6], [0,6,11,6], [0,6,10,6]]
#plot_names = ["spire500","pacs100","wise4","wise3","wise2","wise1"]
ratio_indicies = [[0,6,19,6], [0,6,15,6], [0,6,13,6], [0,6,11,6]]
plot_names = ["spire500","pacs100","wise4","wise2"]
xlims = [[2e-4, 0.6], [2e-4,1], [2e-4,1], [2e-4,1]]
ylims = [[1e-1, 50], [1e-2,5e2], [6e-3,7e1], [5e-2,5]]

# Initialize DustPedia arrays
dp_flux = np.zeros((dp_num, len(band_names))) # first index specifies galaxy, second index specifies band 
dp_err = np.zeros((dp_num, len(band_names)))
dp_bool = np.zeros((dp_num, len(band_names))) # 1 indicates a negative flux (to be colored red) 

# Fill DustPedia arrays
for i in range(dp_num-1):
    g = i+1 # index of current galaxy
    params = dp_file[g].split(',')
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

# Final sample (66)
#galaxies = ['g1.88e10','g1.89e10','g1.90e10','g2.34e10','g2.63e10','g2.64e10','g2.80e10','g2.83e10','g2.94e10',
#            'g3.19e10','g3.44e10','g3.67e10','g3.93e10','g4.27e10','g4.48e10','g4.86e10','g4.94e10','g4.99e10',
#            'g5.05e10','g6.12e10','g6.37e10','g6.77e10','g6.91e10','g7.12e10','g8.89e10','g9.59e10','g1.05e11',
#            'g1.08e11','g1.37e11','g1.52e11','g1.57e11','g1.59e11','g1.64e11','g2.04e11','g2.19e11','g2.39e11',
#            'g2.41e11','g2.42e11','g2.54e11','g3.06e11','g3.21e11','g3.23e11','g3.49e11','g3.55e11','g3.59e11',
#            'g3.71e11','g4.90e11','g5.02e11','g5.31e11','g5.36e11','g5.38e11','g5.46e11','g5.55e11','g6.96e11',
#            'g7.08e11','g7.44e11','g7.55e11','g7.66e11','g8.06e11','g8.13e11','g8.26e11','g8.28e11','g1.12e12',
#            'g1.77e12','g1.92e12','g2.79e12']

# test galaxy
galaxies = ['g7.55e11']

numViews = 100 # (named 0-99)

# Initialize array for SKIRT spatially integrated photometry for color-color plots
SKIRT_flux = np.zeros((len(galaxies), len(band_names)))

os.system('mkdir -p '+plotPath+'SEDs/')
os.system('mkdir -p '+plotPath+'zrgImages/')
for i in range(len(galaxies)):
    plt.figure(figsize=(10,8)) # each galaxy gets an SED plot
    for j in range(numViews):
	    #cube_file = fits.open(SKIRTPath+galaxies[i]+'/sph_broadband'+str(j)+'_total.fits')
	    #cube = np.asarray(cube_file[0].data) # (2000, 2000, 20) cube of broadbands in MJy/sr
	    sed = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_SED'+str(j)+'_sed.dat', unpack = True)
	    wave = sed[0] # spatially integrated SED wavelengths in microns
	    wave = wave * 1e4 # convert to Angstroms
	    spec = sed[1] # spatially integrated SED fluxes in Janskys
	    #bb = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_broadband'+str(j)+'_sed.dat', unpack = True)
	    #bb_wave = bb[0] # spatially integrated broadband wavelengths in microns
	    #bb_wave = bb_wave * 1e4 # convert to Angstroms
	    #bb_spec = bb[1] # spatially integrated broadband fluxes in Janskys
	    #SKIRT_flux[i,:] = bb_spec
	    # SEDs
	    #plt.figure(figsize=(10,8))
	    mask = (wave >= 1e3) & (wave <= 1e7)
	    plt.plot(wave[mask], spec[mask], color='k',alpha=0.5, linewidth=0.1)
	    #plt.xlabel('Wavelength ('+r'$\AA$'+')', fontsize=16)
	    #plt.ylabel('Flux (Janskys)', fontsize=16)
	    #plt.xscale('log')
	    #plt.yscale('log')
	    #plt.xticks(fontsize=16)
	    #plt.yticks(fontsize=16)
	    #plt.savefig(plotPath+'SEDs/'+galaxies[i]+'_SED.png', dpi=300)
	    #plt.close()
	    # Optical images
	    #r_grid = np.asarray(cube[6,:,:]) # sdss_z band (see band_names for indicies)
	    #g_grid = np.asarray(cube[4,:,:]) # sdss_r band
	    #b_grid = np.asarray(cube[3,:,:]) # sdss_g band
	    #fig = plt.figure()
	    #sizes = np.shape(r_grid)
	    #fig.set_size_inches(1. * sizes[0] / sizes[1], 1, forward = False)
	    #ax = plt.Axes(fig, [0., 0., 1., 1.])
	    #ax.set_axis_off()
	    #fig.add_axes(ax)
	    #tot_max = np.amax([np.amax(r_grid),np.amax(g_grid),np.amax(b_grid)])
	    #r_grid = 255 * r_grid / tot_max
	    #g_grid = 255 * g_grid / tot_max
	    #b_grid = 255 * b_grid / tot_max
	    #image = make_lupton_rgb(r_grid, g_grid, b_grid, stretch=3, Q=8) 
	    #ax.imshow(image,interpolation='none')
	    #plt.savefig(plotPath+'zrgImages/'+galaxies[i]+'_zrg'+str(j)+'.png',dpi=sizes[0])
	    #plt.close()
    plt.xlabel('Wavelength ('+r'$\AA$'+')', fontsize=16)
    plt.ylabel('Flux (Janskys)', fontsize=16)
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(plotPath+'SEDs/'+galaxies[i]+'_SED.png', dpi=600)
    plt.close()

for i in range(len(galaxies)):
    os.system('mkdir -p '+plotPath+'zrgImages/'+galaxies[i]+'/')
    for j in range(numViews):
	    cube_file = fits.open(SKIRTPath+galaxies[i]+'/sph_broadband'+str(j)+'_total.fits')
	    cube = np.asarray(cube_file[0].data) # (2000, 2000, 20) cube of broadbands in MJy/sr
	    #sed = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_SED'+str(j)+'_sed.dat', unpack = True)
	    #wave = sed[0] # spatially integrated SED wavelengths in microns
	    #wave = wave * 1e4 # convert to Angstroms
	    #spec = sed[1] # spatially integrated SED fluxes in Janskys
	    #bb = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_broadband'+str(j)+'_sed.dat', unpack = True)
	    #bb_wave = bb[0] # spatially integrated broadband wavelengths in microns
	    #bb_wave = bb_wave * 1e4 # convert to Angstroms
	    #bb_spec = bb[1] # spatially integrated broadband fluxes in Janskys
	    #SKIRT_flux[i,:] = bb_spec
	    # SEDs
	    #plt.figure(figsize=(10,8))
	    #mask = (wave >= 1e3) & (wave <= 1e7)
	    #plt.plot(wave[mask], spec[mask], color='k',alpha=1)
	    #plt.xlabel('Wavelength ('+r'$\AA$'+')', fontsize=16)
	    #plt.ylabel('Flux (Janskys)', fontsize=16)
	    #plt.xscale('log')
	    #plt.yscale('log')
	    #plt.xticks(fontsize=16)
	    #plt.yticks(fontsize=16)
	    #plt.savefig(plotPath+'SEDs/'+galaxies[i]+'_SED.png', dpi=300)
	    #plt.close()
	    # Optical images
	    r_grid = np.asarray(cube[6,:,:]) # sdss_z band (see band_names for indicies)
	    g_grid = np.asarray(cube[4,:,:]) # sdss_r band
	    b_grid = np.asarray(cube[3,:,:]) # sdss_g band
	    fig = plt.figure() # each galaxy and each orientation gets an image
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
	    plt.savefig(plotPath+'zrgImages/'+galaxies[i]+'/zrg'+str(j)+'.png',dpi=sizes[0])
	    plt.close()

exit()


# Color-color plots
for i in range(len(ratio_indicies)): # number of plots to make 
    os.system('mkdir -p '+plotPath+'colorPlots/')
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
    plt.xlim((xlims[i][0], xlims[i][1]))
    plt.ylim((ylims[i][0], ylims[i][1]))
    plt.savefig(plotPath+'colorPlots/'+plot_names[i]+'.png',dpi=600)
    plt.close()

print('done')
