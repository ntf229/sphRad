import numpy as np
import os
from astropy.io import fits
import argparse
import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from timeit import default_timer as timer
import datetime

def stellarMass(galaxy):
    if os.path.isfile(massPath+galaxy+'/stellarMass.npy'):
        SKIRT_stellarMass = float(np.load(massPath+galaxy+'/stellarMass.npy'))
    else:
        stars = np.loadtxt(textPath+galaxy+'/stars.txt')
        youngStars = np.loadtxt(textPath+galaxy+'/youngStars.txt')
        SKIRT_stellarMass = np.sum(stars[:,7]) + (np.sum(youngStars[:,7] * 1.e7)) # in Msun
        os.system('mkdir -p '+massPath+galaxy+'/')
        np.save(massPath+galaxy+'/stellarMass.npy', SKIRT_stellarMass)
    return SKIRT_stellarMass

def SFR(galaxy):
    if os.path.isfile(SFRPath+galaxy+'/SFR.npy'):
        SKIRT_SFR = float(np.load(SFRPath+galaxy+'/SFR.npy'))
    else:
        stars = np.loadtxt(textPath+galaxy+'/stars.txt')
        youngStars = np.loadtxt(textPath+galaxy+'/youngStars.txt')
        youngMask = stars[:,9] < 1e8 # younger than 100 Myrs
        SKIRT_SFR = np.sum(stars[youngMask,7]) + (np.sum(youngStars[:,7] * 1.e7)) / 1.e8 # in Msun per year
        os.system('mkdir -p '+SFRPath+galaxy+'/')
        np.save(SFRPath+galaxy+'/SFR.npy', SKIRT_SFR)
    return SKIRT_SFR

def dustMass(galaxy):
    if os.path.isfile(massPath+galaxy+'/maxTemp'+args.maxTemp+'/metalMass.npy'):
        SKIRT_dustMass = float(np.load(massPath+galaxy+'/maxTemp'+args.maxTemp+'/metalMass.npy')) * float(args.dustFraction)
    else:
        gas = np.loadtxt(textPath+galaxy+'/gas.txt')
        tempMask = gas[:,6] < float(args.maxTemp)
        ghostMask = np.asarray(gas[tempMask][:,4] > 0, dtype=bool) # mask out negative mass ghost particles
        metalMass = np.sum(gas[tempMask, 4][ghostMask] * gas[tempMask, 5][ghostMask]) # in Msun
        #metalMass += np.sum(-1 * gas[tempMask, 4][~ghostMask] * gas[tempMask, 5][~ghostMask]) # add SF particle dust mass
        SKIRT_dustMass = metalMass * float(args.dustFraction)
        os.system('mkdir -p '+massPath+galaxy+'/maxTemp'+args.maxTemp+'/')
        np.save(massPath+galaxy+'/maxTemp'+args.maxTemp+'/metalMass.npy', metalMass)
    return SKIRT_dustMass

def plotSEDs(galaxy, dust=True):
    if dust:
        os.system('mkdir -p '+plotPath+'SEDs/'+galaxy+'/')
        sed = np.loadtxt(SKIRTPath+galaxy+'/sph_SED_'+instName+'_sed.dat', unpack = True)
        tempPlotPath = plotPath
    else:
        os.system('mkdir -p '+noDustPlotPath+'SEDs/'+galaxy+'/')
        sed = np.loadtxt(noDustSKIRTPath+galaxy+'/sph_SED_'+instName+'_sed.dat', unpack = True)
        tempPlotPath = noDustPlotPath
    wave = sed[0] * 1e4 # spatially integrated SED wavelengths converted to Angstroms
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
    plt.savefig(tempPlotPath+'SEDs/'+galaxy+'/'+instName+'_SED.png', dpi=300)
    plt.close()

def plotBothSEDs(galaxy):
    os.system('mkdir -p '+plotPath+'bothSEDs/'+galaxy+'/')
    sed = np.loadtxt(SKIRTPath+galaxy+'/sph_SED_'+instName+'_sed.dat', unpack = True)
    noDustSed = np.loadtxt(noDustSKIRTPath+galaxy+'/sph_SED_'+instName+'_sed.dat', unpack = True)
    wave = sed[0] * 1e4 # spatially integrated SED wavelengths converted to Angstroms
    spec = sed[1] # spatially integrated SED fluxes in Janskys
    noDustSpec = noDustSed[1] # spatially integrated SED fluxes in Janskys
    plt.figure(figsize=(10,8))
    mask = (wave >= 1e3) & (wave <= 1e7)
    plt.plot(wave[mask], spec[mask], color='red',alpha=1)
    plt.plot(wave[mask], noDustSpec[mask], color='blue',alpha=1)
    plt.xlabel(r'$\lambda$'+' ['+r'$\AA$'+']', fontsize=16)
    plt.ylabel(r'$F_{\nu}$'+' [Jy]', fontsize=16)
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(plotPath+'bothSEDs/'+galaxy+'/'+instName+'_SED.png', dpi=300)
    plt.close()

def plotAttenuationCurves(galaxy):
    os.system('mkdir -p '+plotPath+'AttenuationCurves/'+galaxy+'/')
    sed = np.loadtxt(SKIRTPath+galaxy+'/sph_SED_'+instName+'_sed.dat', unpack = True)
    wave = sed[0] * 1e4 # spatially integrated SED wavelengths converted to Angstroms
    spec = sed[1] # spatially integrated SED fluxes in Janskys
    sed = np.loadtxt(noDustSKIRTPath+galaxy+'/sph_SED_'+instName+'_sed.dat', unpack = True)
    noDustSpec = sed[1] 
    att_mask = (wave >= 912) & (wave <= 2e4)
    dustMags = 22.5 - 2.5*np.log10((spec/3631) * 1e9) # convert to Pogson magnitudes
    noDustMags = 22.5 - 2.5*np.log10((noDustSpec/3631) * 1e9) # convert to Pogson magnitudes
    plt.figure(figsize=(10,8))
    plt.plot(wave[att_mask], dustMags[att_mask] - noDustMags[att_mask], color='k', alpha=1)
    plt.xlabel(r'$\lambda$'+' ['+r'$\AA$'+']', fontsize=16)
    plt.ylabel(r'$A_{\lambda}$',fontsize=16)
    plt.xscale('log')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(plotPath+'AttenuationCurves/'+galaxy+'/'+instName+'_attenuation.png', dpi=300)
    plt.close()

def getAvValues():
    AvValues = np.zeros((len(galaxies), numOrientations))
    for i in range(len(galaxies)):
        for j in range(numOrientations):
            instName = 'axisRatio'+str(np.round_(SKIRT_axisRatio[i,j], decimals = 4))
            sed = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_SED_'+instName+'_sed.dat', unpack = True)
            wave = sed[0] * 1e4 # spatially integrated SED wavelengths converted to Angstroms
            spec = sed[1] # spatially integrated SED fluxes in Janskys
            sed = np.loadtxt(noDustSKIRTPath+galaxies[i]+'/sph_SED_'+instName+'_sed.dat', unpack = True)
            noDustSpec = sed[1] 
            att_mask = (wave >= 912) & (wave <= 2e4)
            dustMags = 22.5 - 2.5*np.log10((spec/3631) * 1e9) # convert to Pogson magnitudes
            noDustMags = 22.5 - 2.5*np.log10((noDustSpec/3631) * 1e9) # convert to Pogson magnitudes
            attenuation = dustMags[att_mask] - noDustMags[att_mask]
            Av_index = np.abs(wave[att_mask] - 5500).argmin() # find wave index closes to 5500 angstroms (V)
            AvValues[i,j] =  attenuation[Av_index]
    return AvValues

def plotAllAttenuationCurves():
    os.system('mkdir -p '+plotPath+'AllAttenuationCurves/')
    plt.figure(figsize=(10,8))
    for i in range(len(galaxies)):
        if not sphMassMask[i]:
            continue
        for j in range(numOrientations):
            instName = 'axisRatio'+str(np.round_(SKIRT_axisRatio[i,j], decimals = 4))
            sed = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_SED_'+instName+'_sed.dat', unpack = True)
            wave = sed[0] * 1e4 # spatially integrated SED wavelengths converted to Angstroms
            spec = sed[1] # spatially integrated SED fluxes in Janskys
            sed = np.loadtxt(noDustSKIRTPath+galaxies[i]+'/sph_SED_'+instName+'_sed.dat', unpack = True)
            noDustSpec = sed[1] 
            att_mask = (wave >= 912) & (wave <= 2e4)
            dustMags = 22.5 - 2.5*np.log10((spec/3631) * 1e9) # convert to Pogson magnitudes
            noDustMags = 22.5 - 2.5*np.log10((noDustSpec/3631) * 1e9) # convert to Pogson magnitudes
            attenuation = dustMags[att_mask] - noDustMags[att_mask]
            plt.plot(wave[att_mask], attenuation, color='k', alpha=0.5, linewidth=0.2)
    plt.xlabel(r'$\lambda$'+' ['+r'$\AA$'+']', fontsize=16)
    plt.ylabel(r'$A_{\lambda}$',fontsize=16)
    plt.xscale('log')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(plotPath+'AllAttenuationCurves/attenuation.png', dpi=300)
    plt.close()

def plotAllAttenuationCurvesNorm():
    os.system('mkdir -p '+plotPath+'AllAttenuationCurves/')
    plt.figure(figsize=(10,8))
    for i in range(len(galaxies)):
        if not sphMassMask[i]:
            continue  
        for j in range(numOrientations):
            instName = 'axisRatio'+str(np.round_(SKIRT_axisRatio[i,j], decimals = 4))
            sed = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_SED_'+instName+'_sed.dat', unpack = True)
            wave = sed[0] * 1e4 # spatially integrated SED wavelengths converted to Angstroms
            spec = sed[1] # spatially integrated SED fluxes in Janskys
            sed = np.loadtxt(noDustSKIRTPath+galaxies[i]+'/sph_SED_'+instName+'_sed.dat', unpack = True)
            noDustSpec = sed[1] 
            att_mask = (wave >= 912) & (wave <= 2e4)
            dustMags = 22.5 - 2.5*np.log10((spec/3631) * 1e9) # convert to Pogson magnitudes
            noDustMags = 22.5 - 2.5*np.log10((noDustSpec/3631) * 1e9) # convert to Pogson magnitudes
            attenuation = dustMags[att_mask] - noDustMags[att_mask]
            plt.plot(wave[att_mask], attenuation / SKIRT_AvValues[i,j], color='k', alpha=0.5, linewidth=0.2)
    plt.xlabel(r'$\lambda$'+' ['+r'$\AA$'+']', fontsize=16)
    plt.ylabel(r'$A_{\lambda} / A_{V}$',fontsize=16)
    plt.xscale('log')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(plotPath+'AllAttenuationCurves/normAttenuation.png', dpi=300)
    plt.close()

def plotEnergyBalance(exclude=True):
    os.system('mkdir -p '+plotPath+'EnergyBalance/')
    c = 2.998e18 # speed of light in Anstroms per second
    plt.figure(figsize=(10,8))
    for i in range(len(galaxies)):
        if exclude:
            if not sphMassMask[i]:
                continue  
        for j in range(numOrientations):
            instName = 'axisRatio'+str(np.round_(SKIRT_axisRatio[i,j], decimals = 4))
            sed = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_SED_'+instName+'_sed.dat', unpack = True)
            wave = sed[0] * 1e4 # spatially integrated SED wavelengths converted to Angstroms
            freq = c / wave # in Hz
            spec = sed[1] # spatially integrated SED fluxes in Janskys
            sed = np.loadtxt(noDustSKIRTPath+galaxies[i]+'/sph_SED_'+instName+'_sed.dat', unpack = True)
            noDustSpec = sed[1] 
            att_mask = wave <= 2e4
            emit_mask = wave > 2e4
            attenuation = noDustSpec[att_mask] - spec[att_mask] # in Janskys
            emission = spec[emit_mask] - noDustSpec[emit_mask]
            attEnergy = np.trapz(attenuation, freq[att_mask])
            emitEnergy = np.trapz(emission, freq[emit_mask])
            plt.scatter(np.log10(SKIRT_stellarMass[i]), attEnergy / emitEnergy, color='k', alpha=0.5, s=15)
    plt.axhline(y=1, color='k')
    plt.xlabel('log(Stellar Mass / '+r'$M_{\odot}$'+')', fontsize=16)
    plt.ylabel('Attenuated Energy / Emitted Energy',fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    if exclude:
        plt.savefig(plotPath+'EnergyBalance/energyBalanceMassCut.png', dpi=300)
    else:
        plt.savefig(plotPath+'EnergyBalance/energyBalance.png', dpi=300)
    plt.close()

def plotAvAxisRatio(exclude=True):
    os.system('mkdir -p '+plotPath+'AvPlots/')
    plt.figure(figsize=(10,8))
    if exclude:
        plt.scatter(SKIRT_axisRatio[sphMassMask,:], SKIRT_AvValues[sphMassMask,:], color='k', alpha=0.5, s=15)
    else:
        plt.scatter(SKIRT_axisRatio, SKIRT_AvValues, color='k', alpha=0.5, s=15)
    plt.xlabel('Axis Ratio', fontsize=16)
    plt.ylabel(r'$A_{V}$',fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    if exclude:
        plt.savefig(plotPath+'AvPlots/AvAxisRatioMassCut.png', dpi=300)
    else:
        plt.savefig(plotPath+'AvPlots/AvAxisRatio.png', dpi=300)
    plt.close()

def plotAvStellarMass(exclude=True):
    os.system('mkdir -p '+plotPath+'AvPlots/')
    plt.figure(figsize=(10,8))
    for i in range(len(galaxies)):
        if exclude:
            if not sphMassMask[i]:
                continue
        for j in range(numOrientations):
            plt.scatter(np.log10(SKIRT_stellarMass[i]), SKIRT_AvValues[i,j], color='k', alpha=0.5, s=15)
    plt.xlabel('log(Stellar Mass / '+r'$M_{\odot}$'+')', fontsize=16)
    plt.ylabel(r'$A_{V}$',fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    if exclude:
        plt.savefig(plotPath+'AvPlots/AvStellarMassMassCut.png', dpi=300)
    else:
        plt.savefig(plotPath+'AvPlots/AvStellarMass.png', dpi=300)
    plt.close()

def plotAvsSFR(exclude=True):
    os.system('mkdir -p '+plotPath+'AvPlots/')
    plt.figure(figsize=(10,8))
    for i in range(len(galaxies)):
        if exclude:
            if not sphMassMask[i]:
                continue
        for j in range(numOrientations):
            plt.scatter(np.log10(SKIRT_sSFR[i]), SKIRT_AvValues[i,j], color='k', alpha=0.5, s=15)
    plt.xlabel('log(sSFR / '+r'$yr^{-1}$'+')', fontsize=16)
    plt.ylabel(r'$A_{V}$',fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    if exclude:
        plt.savefig(plotPath+'AvPlots/AvsSFRMassCut.png', dpi=300)
    else:
        plt.savefig(plotPath+'AvPlots/AvsSFR.png', dpi=300)
    plt.close()

def plotAvDustToStellar(exclude=True):
    os.system('mkdir -p '+plotPath+'AvPlots/')
    plt.figure(figsize=(10,8))
    for i in range(len(galaxies)):
        if exclude:
            if not sphMassMask[i]:
                continue
        for j in range(numOrientations):
            plt.scatter(np.log10(SKIRT_dustToStellar[i]), SKIRT_AvValues[i,j], color='k', alpha=0.5, s=15)
    plt.xlabel('log(Dust Mass / Stellar Mass)', fontsize=16)
    plt.ylabel(r'$A_{V}$',fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    if exclude:
        plt.savefig(plotPath+'AvPlots/AvDustToStellarMassCut.png', dpi=300)
    else:
        plt.savefig(plotPath+'AvPlots/AvDustToStellar.png', dpi=300)
    plt.close()

def makeImages(galaxy):
    os.system('mkdir -p '+plotPath+'zrgImages/'+galaxy+'/')
    cube_file = fits.open(SKIRTPath+galaxy+'/sph_broadband_'+instName+'_total.fits')
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
    ax.imshow(image, interpolation='none')
    plt.savefig(plotPath+'zrgImages/'+galaxy+'/'+instName+'_zrg.png', dpi=sizes[0])
    plt.close()

def colorColorPlots(i, exclude=True):
    #if os.path.isfile(plotPath+'colorPlots/'+plot_names[i]+'.png'): # skip if already done
    #    print('color color plots already done')
    #    return
    plt.figure(figsize=(10,8))
    # DustPedia propogation of errors
    xerr = np.absolute(dp_flux[:,ratio_indicies[i][0]] / dp_flux[:,ratio_indicies[i][1]]) * np.sqrt( 
                      (dp_err[:,ratio_indicies[i][0]] / dp_flux[:,ratio_indicies[i][0]])**2 + 
                      (dp_err[:,ratio_indicies[i][1]] / dp_flux[:,ratio_indicies[i][1]])**2 )
    yerr = np.absolute(dp_flux[:,ratio_indicies[i][2]] / dp_flux[:,ratio_indicies[i][3]]) * np.sqrt( 
                      (dp_err[:,ratio_indicies[i][2]] / dp_flux[:,ratio_indicies[i][2]])**2 + 
                      (dp_err[:,ratio_indicies[i][3]] / dp_flux[:,ratio_indicies[i][3]])**2 )
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
    #for j in range(numOrientations):
    #    plt.scatter(SKIRT_flux[:,j,ratio_indicies[i][0]] / SKIRT_flux[:,j,ratio_indicies[i][1]], 
    #                SKIRT_flux[:,j,ratio_indicies[i][2]] / SKIRT_flux[:,j,ratio_indicies[i][3]], 
    #                marker='o', s=20, zorder=10,alpha=0.7, c='blue')
    if exclude:
        plt.scatter(SKIRT_flux[sphMassMask,:,ratio_indicies[i][0]] / SKIRT_flux[sphMassMask,:,ratio_indicies[i][1]], 
                    SKIRT_flux[sphMassMask,:,ratio_indicies[i][2]] / SKIRT_flux[sphMassMask,:,ratio_indicies[i][3]], 
                    marker='o', s=20, zorder=10,alpha=0.7, c='blue')
    else:
        plt.scatter(SKIRT_flux[sphMassMask,:,ratio_indicies[i][0]] / SKIRT_flux[sphMassMask,:,ratio_indicies[i][1]], 
                    SKIRT_flux[sphMassMask,:,ratio_indicies[i][2]] / SKIRT_flux[sphMassMask,:,ratio_indicies[i][3]], 
                    marker='o', s=20, zorder=10,alpha=0.7, c='blue')
    plt.xlabel(band_names[ratio_indicies[i][0]]+' / '+band_names[ratio_indicies[i][1]], fontsize=16)
    plt.ylabel(band_names[ratio_indicies[i][2]]+' / '+band_names[ratio_indicies[i][3]], fontsize=16)
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.xlim((xlims[i][0], xlims[i][1]))
    #plt.ylim((ylims[i][0], ylims[i][1]))
    if exclude:
        plt.savefig(plotPath+'colorPlots/'+plot_names[i]+'MassCut.png',dpi=300)
    else:
        plt.savefig(plotPath+'colorPlots/'+plot_names[i]+'.png',dpi=300)
    plt.close()

def plotEachAxisRatioVsColor(i):
    for j in range(len(galaxies)):
        os.system('mkdir -p '+plotPath+'sphAxisRatioColorPlots/'+galaxies[j]+'/')
        plt.figure(figsize=(10,8))   
        plt.scatter(SKIRT_axisRatio[j,:], 
                    np.log10(SKIRT_flux[j,:,axisRatio_color_indices[i][0]] / 
                    SKIRT_flux[j,:,axisRatio_color_indices[i][1]]), 
                    marker='o', s=20, zorder=10, alpha=0.7, c='blue')
        plt.xlabel('Axis Ratio', fontsize=16)
        plt.ylabel('log('+band_names[axisRatio_color_indices[i][0]]+' / '+
                   band_names[axisRatio_color_indices[i][1]]+')', fontsize=16)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.savefig(plotPath+'sphAxisRatioColorPlots/'+galaxies[j]+'/'+axisRatio_plot_names[i]+'.png',dpi=300)
        plt.close()

def fitBysSFR(i):
    # polyfit does not handle NaNs, need to mask them out first
    dpPhotMask = ~np.asarray(dp_bool[diskMask][faceOnMask][:, axisRatio_color_indices[i][0]].tolist() or 
                 dp_bool[diskMask][faceOnMask][:, axisRatio_color_indices[i][1]].tolist())
    dpFaceOnsSFR = np.log10(dp_SFR[diskMask][faceOnMask][dpPhotMask] / 
                   dp_stellarMass[diskMask][faceOnMask][dpPhotMask])
    dpFaceOnFluxRatios = np.log10((dp_flux[diskMask][faceOnMask][dpPhotMask][:, axisRatio_color_indices[i][0]] / 
                         dp_flux[diskMask][faceOnMask][dpPhotMask][:, axisRatio_color_indices[i][1]]))
    nanMask = np.isfinite(dpFaceOnsSFR) & np.isfinite(dpFaceOnFluxRatios)
    dpLinearFit = np.polyfit(dpFaceOnsSFR[nanMask], dpFaceOnFluxRatios[nanMask], 1)
    dpFit = np.poly1d(dpLinearFit)
    # make plots showing fits
    plt.figure(figsize=(10,8))
    plt.scatter(dpFaceOnsSFR[nanMask], dpFaceOnFluxRatios[nanMask], marker='o', s=15, color='k')
    x_fit = np.linspace(np.amin(dpFaceOnsSFR[nanMask]), np.amax(dpFaceOnsSFR[nanMask]), num=2)
    y_fit = dpFit(x_fit)
    plt.plot(x_fit, y_fit, color='k')
    flat_sphFaceOnsSFR = np.asarray([])
    flat_sphFaceOnFluxRatios = np.asarray([])
    # flatten arrays (each galaxy may have multiple face-on orientations)
    for j in range(len(galaxies)):
        if not sphMassMask[j]: 
            continue
        current_sphFaceOnFluxRatios = np.log10(SKIRT_flux[j,sphFaceOnMask[j],axisRatio_color_indices[i][0]] /
                                      SKIRT_flux[j,sphFaceOnMask[j],axisRatio_color_indices[i][1]]) 
        current_sphFaceOnsSFR = np.repeat(np.log10(SKIRT_sSFR[j, np.newaxis]), 
                                len(current_sphFaceOnFluxRatios), axis=0)
        flat_sphFaceOnsSFR = np.append(flat_sphFaceOnsSFR, current_sphFaceOnsSFR)
        flat_sphFaceOnFluxRatios = np.append(flat_sphFaceOnFluxRatios, current_sphFaceOnFluxRatios)
    nanMask = np.isfinite(flat_sphFaceOnsSFR) & np.isfinite(flat_sphFaceOnFluxRatios)
    sphLinearFit = np.polyfit(flat_sphFaceOnsSFR[nanMask], flat_sphFaceOnFluxRatios[nanMask], 1)
    sphFit = np.poly1d(sphLinearFit) 
    plt.scatter(flat_sphFaceOnsSFR[nanMask], flat_sphFaceOnFluxRatios[nanMask], marker='o', s=15, color='blue')
    x_fit = np.linspace(np.amin(flat_sphFaceOnsSFR[nanMask]), np.amax(flat_sphFaceOnsSFR[nanMask]), num=2)
    y_fit = sphFit(x_fit)
    plt.plot(x_fit, y_fit, color='blue')
    plt.xlabel('log(sSFR)', fontsize=16)
    plt.ylabel('log('+band_names[axisRatio_color_indices[i][0]]+' / '+
               band_names[axisRatio_color_indices[i][1]]+')', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(plotPath+'fitPlots/'+axisRatio_plot_names[i]+'.png',dpi=300)
    plt.close()
    dpAvg = dpFit(np.log10(dp_SFR / dp_stellarMass))
    #sphAvg = np.repeat(sphFit(np.log10(SKIRT_sSFR))[:, np.newaxis], numOrientations, axis=1)
    sphAvg = np.repeat(sphFit(np.log10(SKIRT_sSFR))[sphMassMask, np.newaxis], numOrientations, axis=1)
    return dpAvg, sphAvg

def fitByDustToStellar(i):
    # polyfit does not handle NaNs, need to mask them out first
    dpPhotMask = ~np.asarray(dp_bool[diskMask][faceOnMask][:, axisRatio_color_indices[i][0]].tolist() or 
                 dp_bool[diskMask][faceOnMask][:, axisRatio_color_indices[i][1]].tolist())
    dpFaceOnDustToStellar = np.log10(dp_dustMass[diskMask][faceOnMask][dpPhotMask] / 
                            dp_stellarMass[diskMask][faceOnMask][dpPhotMask])
    dpFaceOnFluxRatios = np.log10((dp_flux[diskMask][faceOnMask][dpPhotMask][:, axisRatio_color_indices[i][0]] / 
                         dp_flux[diskMask][faceOnMask][dpPhotMask][:, axisRatio_color_indices[i][1]]))
    nanMask = np.isfinite(dpFaceOnDustToStellar) & np.isfinite(dpFaceOnFluxRatios)
    dpLinearFit = np.polyfit(dpFaceOnDustToStellar[nanMask], dpFaceOnFluxRatios[nanMask], 1)
    dpFit = np.poly1d(dpLinearFit) 
    # make plots showing fits
    plt.figure(figsize=(10,8))
    plt.scatter(dpFaceOnDustToStellar[nanMask], dpFaceOnFluxRatios[nanMask], marker='o', s=15, color='k')
    x_fit = np.linspace(np.amin(dpFaceOnDustToStellar[nanMask]), np.amax(dpFaceOnDustToStellar[nanMask]), num=2)
    y_fit = dpFit(x_fit)
    plt.plot(x_fit, y_fit, color='k')
    flat_sphFaceOnDustToStellar = np.asarray([])
    flat_sphFaceOnFluxRatios = np.asarray([])
    # flatten arrays (each galaxy may have multiple face-on orientations)
    for j in range(len(galaxies)):
        current_sphFaceOnFluxRatios = np.log10(SKIRT_flux[j,sphFaceOnMask[j],axisRatio_color_indices[i][0]] /
                                      SKIRT_flux[j,sphFaceOnMask[j],axisRatio_color_indices[i][1]]) 
        current_sphFaceOnDustToStellar = np.repeat(np.log10(SKIRT_dustToStellar[j, np.newaxis]), 
                                         len(current_sphFaceOnFluxRatios), axis=0)
        flat_sphFaceOnDustToStellar = np.append(flat_sphFaceOnDustToStellar, current_sphFaceOnDustToStellar)
        flat_sphFaceOnFluxRatios = np.append(flat_sphFaceOnFluxRatios, current_sphFaceOnFluxRatios)
    nanMask = np.isfinite(flat_sphFaceOnDustToStellar) & np.isfinite(flat_sphFaceOnFluxRatios)
    sphLinearFit = np.polyfit(flat_sphFaceOnDustToStellar[nanMask], flat_sphFaceOnFluxRatios[nanMask], 1)
    sphFit = np.poly1d(sphLinearFit) 
    plt.scatter(flat_sphFaceOnDustToStellar[nanMask], flat_sphFaceOnFluxRatios[nanMask], marker='o', s=15, color='blue')
    x_fit = np.linspace(np.amin(flat_sphFaceOnDustToStellar[nanMask]), np.amax(flat_sphFaceOnDustToStellar[nanMask]), num=2)
    y_fit = sphFit(x_fit)
    plt.plot(x_fit, y_fit, color='blue')
    plt.xlabel('log(Dust Mass / Stellar Mass)', fontsize=16)
    plt.ylabel('log('+band_names[axisRatio_color_indices[i][0]]+' / '+
               band_names[axisRatio_color_indices[i][1]]+')', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(plotPath+'fitPlots/'+axisRatio_plot_names[i]+'.png',dpi=300)
    plt.close()
    dpAvg = dpFit(np.log10(dp_dustMass / dp_stellarMass))
    sphAvg = np.repeat(sphFit(np.log10(SKIRT_dustToStellar))[:, np.newaxis], numOrientations, axis=1)
    return dpAvg, sphAvg

def fitByStellar(i):
    # polyfit does not handle NaNs, need to mask them out first
    dpPhotMask = ~np.asarray(dp_bool[diskMask][faceOnMask][:, axisRatio_color_indices[i][0]].tolist() or 
                 dp_bool[diskMask][faceOnMask][:, axisRatio_color_indices[i][1]].tolist())
    dpFaceOnStellar = np.log10(dp_stellarMass[diskMask][faceOnMask][dpPhotMask])
    dpFaceOnFluxRatios = np.log10((dp_flux[diskMask][faceOnMask][dpPhotMask][:, axisRatio_color_indices[i][0]] / 
                         dp_flux[diskMask][faceOnMask][dpPhotMask][:, axisRatio_color_indices[i][1]]))
    nanMask = np.isfinite(dpFaceOnStellar) & np.isfinite(dpFaceOnFluxRatios)
    dpLinearFit = np.polyfit(dpFaceOnStellar[nanMask], dpFaceOnFluxRatios[nanMask], 1)
    dpFit = np.poly1d(dpLinearFit) 
    # make plots showing fits
    plt.figure(figsize=(10,8))
    plt.scatter(dpFaceOnStellar[nanMask], dpFaceOnFluxRatios[nanMask], marker='o', s=15, color='k')
    x_fit = np.linspace(np.amin(dpFaceOnStellar[nanMask]), np.amax(dpFaceOnStellar[nanMask]), num=2)
    y_fit = dpFit(x_fit)
    plt.plot(x_fit, y_fit, color='k')
    flat_sphFaceOnStellar = np.asarray([])
    flat_sphFaceOnFluxRatios = np.asarray([])
    # flatten arrays (each galaxy may have multiple face-on orientations)
    for j in range(len(galaxies)):
        current_sphFaceOnFluxRatios = np.log10(SKIRT_flux[j,sphFaceOnMask[j],axisRatio_color_indices[i][0]] /
                                      SKIRT_flux[j,sphFaceOnMask[j],axisRatio_color_indices[i][1]]) 
        current_sphFaceOnStellar = np.repeat(np.log10(SKIRT_stellarMass[j, np.newaxis]), 
                                         len(current_sphFaceOnFluxRatios), axis=0)
        flat_sphFaceOnStellar = np.append(flat_sphFaceOnStellar, current_sphFaceOnStellar)
        flat_sphFaceOnFluxRatios = np.append(flat_sphFaceOnFluxRatios, current_sphFaceOnFluxRatios)
    nanMask = np.isfinite(flat_sphFaceOnStellar) & np.isfinite(flat_sphFaceOnFluxRatios)
    sphLinearFit = np.polyfit(flat_sphFaceOnStellar[nanMask], flat_sphFaceOnFluxRatios[nanMask], 1)
    sphFit = np.poly1d(sphLinearFit) 
    plt.scatter(flat_sphFaceOnStellar[nanMask], flat_sphFaceOnFluxRatios[nanMask], marker='o', s=15, color='blue')
    x_fit = np.linspace(np.amin(flat_sphFaceOnStellar[nanMask]), np.amax(flat_sphFaceOnStellar[nanMask]), num=2)
    y_fit = sphFit(x_fit)
    plt.plot(x_fit, y_fit, color='blue')
    plt.xlabel('log(Stellar Mass)', fontsize=16)
    plt.ylabel('log('+band_names[axisRatio_color_indices[i][0]]+' / '+
               band_names[axisRatio_color_indices[i][1]]+')', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(plotPath+'fitPlots/'+axisRatio_plot_names[i]+'.png',dpi=300)
    plt.close()
    dpAvg = dpFit(np.log10(dp_stellarMass))
    sphAvg = np.repeat(sphFit(np.log10(SKIRT_stellarMass))[:, np.newaxis], numOrientations, axis=1)
    return dpAvg, sphAvg

def fitByFluxRatios(i):
    # polyfit does not handle NaNs, need to mask them out first 
    dpPhotMask1 = ~np.asarray(dp_bool[diskMask][faceOnMask][:, axisRatio_color_indices[i][0]].tolist() or 
                  dp_bool[diskMask][faceOnMask][:, axisRatio_color_indices[i][1]].tolist())
    dpPhotMask2 = ~np.asarray(dp_bool[diskMask][faceOnMask][dpPhotMask1][:, fitFluxes[0]].tolist() or 
                  dp_bool[diskMask][faceOnMask][dpPhotMask1][:, fitFluxes[1]].tolist())
    dpFaceOnFitRatio = np.log10((dp_flux[diskMask][faceOnMask][dpPhotMask1][dpPhotMask2][:, fitFluxes[0]] / 
                       dp_flux[diskMask][faceOnMask][dpPhotMask1][dpPhotMask2][:, fitFluxes[1]]))
    dpFaceOnFluxRatios = np.log10((dp_flux[diskMask][faceOnMask][dpPhotMask1][dpPhotMask2][:, axisRatio_color_indices[i][0]] / 
                         dp_flux[diskMask][faceOnMask][dpPhotMask1][dpPhotMask2][:, axisRatio_color_indices[i][1]]))
    nanMask = np.isfinite(dpFaceOnFitRatio) & np.isfinite(dpFaceOnFluxRatios)
    dpLinearFit = np.polyfit(dpFaceOnFitRatio[nanMask], dpFaceOnFluxRatios[nanMask], 1)
    dpFit = np.poly1d(dpLinearFit) 
    # make plots showing fits
    plt.figure(figsize=(10,8))
    plt.scatter(dpFaceOnFitRatio[nanMask], dpFaceOnFluxRatios[nanMask], marker='o', s=15, color='k')
    x_fit = np.linspace(np.amin(dpFaceOnFitRatio[nanMask]), np.amax(dpFaceOnFitRatio[nanMask]), num=2)
    y_fit = dpFit(x_fit)
    plt.plot(x_fit, y_fit, color='k')
    flat_sphFaceOnFitRatio = np.asarray([])
    flat_sphFaceOnFluxRatios = np.asarray([])
    # flatten arrays (each galaxy may have multiple face-on orientations)
    for j in range(len(galaxies)):
        current_sphFaceOnFitRatio = np.log10(SKIRT_flux[j,sphFaceOnMask[j],fitFluxes[0]] /
                                    SKIRT_flux[j,sphFaceOnMask[j],fitFluxes[1]])
        current_sphFaceOnFluxRatios = np.log10(SKIRT_flux[j,sphFaceOnMask[j],axisRatio_color_indices[i][0]] /
                                      SKIRT_flux[j,sphFaceOnMask[j],axisRatio_color_indices[i][1]]) 
        flat_sphFaceOnFitRatio = np.append(flat_sphFaceOnFitRatio, current_sphFaceOnFitRatio)
        flat_sphFaceOnFluxRatios = np.append(flat_sphFaceOnFluxRatios, current_sphFaceOnFluxRatios)
    nanMask = np.isfinite(flat_sphFaceOnFitRatio) & np.isfinite(flat_sphFaceOnFluxRatios)
    sphLinearFit = np.polyfit(flat_sphFaceOnFitRatio[nanMask], flat_sphFaceOnFluxRatios[nanMask], 1)
    sphFit = np.poly1d(sphLinearFit) 
    plt.scatter(flat_sphFaceOnFitRatio[nanMask], flat_sphFaceOnFluxRatios[nanMask], marker='o', s=15, color='blue')
    x_fit = np.linspace(np.amin(flat_sphFaceOnFitRatio[nanMask]), np.amax(flat_sphFaceOnFitRatio[nanMask]), num=2)
    y_fit = sphFit(x_fit)
    plt.plot(x_fit, y_fit, color='blue')
    plt.xlabel('log('+band_names[fitFluxes[0]]+' / '+band_names[fitFluxes[1]]+')', fontsize=16)
    plt.ylabel('log('+band_names[axisRatio_color_indices[i][0]]+' / '+
               band_names[axisRatio_color_indices[i][1]]+')', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(plotPath+'fitPlots/'+axisRatio_plot_names[i]+'.png',dpi=300)
    plt.close()
    dpAvg = dpFit(np.log10(dp_flux[:,fitFluxes[0]] / dp_flux[:,fitFluxes[1]]))
    sphAvg = sphFit(np.log10(SKIRT_flux[:,:,fitFluxes[0]] / SKIRT_flux[:,:,fitFluxes[1]]))
    return dpAvg, sphAvg

def sSFRnoDust(i):
    # polyfit does not handle NaNs, need to mask them out first
    dpPhotMask = ~np.asarray(dp_bool[diskMask][faceOnMask][:, axisRatio_color_indices[i][0]].tolist() or 
                 dp_bool[diskMask][faceOnMask][:, axisRatio_color_indices[i][1]].tolist())
    dpFaceOnsSFR = np.log10(dp_SFR[diskMask][faceOnMask][dpPhotMask] / 
                   dp_stellarMass[diskMask][faceOnMask][dpPhotMask])
    dpFaceOnFluxRatios = np.log10((dp_flux[diskMask][faceOnMask][dpPhotMask][:, axisRatio_color_indices[i][0]] / 
                         dp_flux[diskMask][faceOnMask][dpPhotMask][:, axisRatio_color_indices[i][1]]))
    nanMask = np.isfinite(dpFaceOnsSFR) & np.isfinite(dpFaceOnFluxRatios)
    dpLinearFit = np.polyfit(dpFaceOnsSFR[nanMask], dpFaceOnFluxRatios[nanMask], 1)
    dpFit = np.poly1d(dpLinearFit)
    # make plots showing fits
    plt.figure(figsize=(10,8))
    plt.scatter(dpFaceOnsSFR[nanMask], dpFaceOnFluxRatios[nanMask], marker='o', s=15, color='k')
    x_fit = np.linspace(np.amin(dpFaceOnsSFR[nanMask]), np.amax(dpFaceOnsSFR[nanMask]), num=2)
    y_fit = dpFit(x_fit)
    plt.plot(x_fit, y_fit, color='k')
    flat_sphFaceOnsSFR = np.asarray([])
    flat_sphFaceOnFluxRatios = np.asarray([])
    # flatten arrays (each galaxy may have multiple face-on orientations)
    for j in range(len(galaxies)):
        current_sphFaceOnFluxRatios = np.log10(noDustSKIRT_flux[j,sphFaceOnMask[j],axisRatio_color_indices[i][0]] /
                                      noDustSKIRT_flux[j,sphFaceOnMask[j],axisRatio_color_indices[i][1]]) 
        current_sphFaceOnsSFR = np.repeat(np.log10(SKIRT_sSFR[j, np.newaxis]), 
                                len(current_sphFaceOnFluxRatios), axis=0)
        flat_sphFaceOnsSFR = np.append(flat_sphFaceOnsSFR, current_sphFaceOnsSFR)
        flat_sphFaceOnFluxRatios = np.append(flat_sphFaceOnFluxRatios, current_sphFaceOnFluxRatios)
    nanMask = np.isfinite(flat_sphFaceOnsSFR) & np.isfinite(flat_sphFaceOnFluxRatios)
    sphLinearFit = np.polyfit(flat_sphFaceOnsSFR[nanMask], flat_sphFaceOnFluxRatios[nanMask], 1)
    sphFit = np.poly1d(sphLinearFit) 
    plt.scatter(flat_sphFaceOnsSFR[nanMask], flat_sphFaceOnFluxRatios[nanMask], marker='o', s=15, color='blue')
    x_fit = np.linspace(np.amin(flat_sphFaceOnsSFR[nanMask]), np.amax(flat_sphFaceOnsSFR[nanMask]), num=2)
    y_fit = sphFit(x_fit)
    plt.plot(x_fit, y_fit, color='blue')
    plt.xlabel('log(sSFR)', fontsize=16)
    plt.ylabel('log('+band_names[axisRatio_color_indices[i][0]]+' / '+
               band_names[axisRatio_color_indices[i][1]]+')', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(noDustPlotPath+'sSFRPlots/'+axisRatio_plot_names[i]+'.png',dpi=300)
    plt.close()

parser = argparse.ArgumentParser()
parser.add_argument("--ageSmooth") # if True, smooth ages based on number of star particles (makeTextFiles parameter)
parser.add_argument("--SF") # if True, star particles younger than 10 Myrs are assigned MAPPINGS-III SEDs (makeTextFiles parameter)
parser.add_argument("--tauClear") # clearing time in Myrs for MAPPINGS-III f_PDR calculations (only matters if SF=True) (makeTextFiles parameter)
parser.add_argument("--clumps") # if True, add subgrid clumpiness to gas near MAPPINGS-III particles (makeTextFiles parameter)
parser.add_argument("--numCells") # number of cells along one dimension for subgrid clumping (only matters if clumps=True) (makeTextFiles parameter)
parser.add_argument("--numClumps") # number of clumps for subgrid clumping (only matters if clumps=True) (makeTextFiles parameter)
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
SKIRTPath = resultPath+'resources/selectedOrientations_SKIRT/'
noDustSKIRTPath = resultPath+'resources/selectedOrientations_SKIRT/'
plotPath = resultPath+'resources/selectedPlots/'
massPath = resultPath+'resources/NIHAO/GlobalProps/stellarMasses/'
SFRPath = resultPath+'resources/NIHAO/GlobalProps/SFR/'
selectedPath = resultPath+'resources/selectedOrientations/'
if eval(args.ageSmooth):
    textPath += 'ageSmooth/'
    SKIRTPath += 'ageSmooth/'
    noDustSKIRTPath += 'ageSmooth/'
    plotPath += 'ageSmooth/'
    SFRPath += 'ageSmooth/'
else:
    textPath += 'noAgeSmooth/'
    SKIRTPath += 'noAgeSmooth/'
    noDustSKIRTPath += 'noAgeSmooth/'
    plotPath += 'noAgeSmooth/'
    SFRPath += 'noAgeSmooth/'
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
noDustPlotPath = plotPath
SKIRTPath += 'dust/dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'
plotPath += 'dust/dustFraction'+args.dustFraction+'/maxTemp'+args.maxTemp+'/'
noDustSKIRTPath += 'noDust/'
noDustPlotPath += 'noDust/'
if eval(args.clumps):
    textPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
    SKIRTPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
    plotPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
else:
    textPath += 'noClumps/'
    SKIRTPath += 'noClumps/'
    plotPath += 'noClumps/'

SKIRTPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
plotPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
noDustSKIRTPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'
noDustPlotPath += 'numPhotons'+args.numPhotons+'/'+args.SSP+'/'

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
ratio_indicies = [[0,6,7,6], [0,6,8,6], [0,6,9,6], [0,6,10,6], [0,6,11,6], [0,6,12,6], [0,6,13,6], 
                  [0,6,14,6], [0,6,15,6], [0,6,16,6], [0,6,17,6], [0,6,18,6], [0,6,19,6]]
plot_names = ["J","H","K","W1","W2","W3","W4","pacs70","pacs100","pacs160","spire250","spire350","spire500"]

# axis ratio color plot info
axisRatio_color_indices = [[0,7], [5,7]] # relative to J
#axisRatio_color_indices = [[0,6], [5,6]] # relative to z
axisRatio_plot_names = ["FUV", "i"]

# flux ratio to fit average axis ratio colors 
fitFluxes = [15, 7] # PACS100 / J 
#fitFluxes = [15, 6] # PACS100 / z
#fitFluxes = [13, 7] # W4 / J 
#fitFluxes = [11, 7] # W2 / J 

# Initialize DustPedia arrays
dp_flux = np.zeros((dp_num, len(band_names))) # first index specifies galaxy, second index specifies band 
dp_err = np.zeros((dp_num, len(band_names)))
dp_bool = np.zeros((dp_num, len(band_names)), dtype=bool) # True means flux consistent with 0 (to be colored red) 
dp_axisRatio = np.zeros(dp_num)
dp_disk = np.zeros(dp_num) # 0 means not disk, 1 means disk (t >= 0)
dp_stellarMass = np.zeros(dp_num)
dp_dustMass = np.zeros(dp_num)
dp_SFR = np.zeros(dp_num)

dp3_names = [] 
for i in range(dp_num):
    g3 = i+1 # skip first line of file (headers)
    dp3_names.append(dp_file3[g3].split(',')[0])
dp3_names = np.asarray(dp3_names) # doesn't include first line of dp3_file

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
    dp_SFR[i] = params3[1] # in solar masses per year (CIGALE)
    for j in range(len(band_names)):
        if params[int(dp_flux_index[j])]:
            dp_flux[i,j] = float(params[int(dp_flux_index[j])])
            dp_err[i,j] = float(params[int(dp_err_index[j])])
            if dp_flux[i,j] - 2*dp_err[i,j] <= 0:
                dp_flux[i,j] = 2 * dp_err[i,j] # if flux is consistent with 0, set to 2*sigma 
                dp_bool[i,j] = True # signals flux is consistent with 0 (upper limit)
            else:
                dp_bool[i,j] = False
        else:
            dp_flux[i,j] = float("NaN") # set flux and error to NaN if photometry not available 
            dp_err[i,j] = float("NaN")
            dp_bool[i,j] = True 

diskMask = dp_disk == 1
faceOnMask = dp_axisRatio[diskMask] > 0.85 

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
numOrientations = 10
SKIRT_flux = np.zeros((len(galaxies), numOrientations, len(band_names)))
noDustSKIRT_flux = np.zeros((len(galaxies), numOrientations, len(band_names)))
SKIRT_axisRatio = np.zeros((len(galaxies), numOrientations))
SKIRT_stellarMass = np.zeros(len(galaxies))
SKIRT_SFR = np.zeros(len(galaxies))
SKIRT_dustMass = np.zeros(len(galaxies))
sphFaceOnMask = np.zeros((len(galaxies), numOrientations), dtype=bool)
sphAvMask = np.zeros(len(galaxies), dtype=bool)

start = timer()

os.system('mkdir -p '+plotPath+'SEDs/')
os.system('mkdir -p '+noDustPlotPath+'SEDs/')
os.system('mkdir -p '+plotPath+'zrgImages/')
os.system('mkdir -p '+plotPath+'AttenuationCurves/')
for i in range(len(galaxies)):
    print('starting '+galaxies[i])
    SKIRT_stellarMass[i] = stellarMass(galaxies[i])
    SKIRT_SFR[i] = SFR(galaxies[i])
    SKIRT_dustMass[i] = dustMass(galaxies[i])
    # import axis ratios from table
    selections = np.load(selectedPath+galaxies[i]+'/selectedAxisRatios.npy') # [inc, az, axisRatio]
    SKIRT_axisRatio[i,:] = selections[:,2]
    sphFaceOnMask[i,:] = SKIRT_axisRatio[i,:] > 0.85
    for j in range(numOrientations):
        instName = 'axisRatio'+str(np.round_(SKIRT_axisRatio[i,j], decimals = 4))
        bb = np.loadtxt(SKIRTPath+galaxies[i]+'/sph_broadband_'+instName+'_sed.dat', unpack = True)
        SKIRT_flux[i,j,:] = bb[1] # spatially integrated broadband fluxes in Janskys
        noDustbb = np.loadtxt(noDustSKIRTPath+galaxies[i]+'/sph_broadband_'+instName+'_sed.dat', unpack = True)
        noDustSKIRT_flux[i,j,:] = noDustbb[1] # spatially integrated broadband fluxes in Janskys
        # skip this part of analysis if already done
        if os.path.isfile(plotPath+'zrgImages/'+galaxies[i]+'/'+instName+'_zrg.png'):
            print('skipping first analysis')
            continue
        #plotSEDs(galaxies[i])
        #plotSEDs(galaxies[i], dust=False)
        #plotAttenuationCurves(galaxies[i])
        #makeImages(galaxies[i]) 
        #plotBothSEDs(galaxies[i]) # dust and noDust

SKIRT_sSFR = SKIRT_SFR / SKIRT_stellarMass
SKIRT_dustToStellar = SKIRT_dustMass / SKIRT_stellarMass
SKIRT_AvValues = getAvValues()
sphMassMask = np.log10(SKIRT_stellarMass) > 9.5 

for i in range(len(galaxies)):
    if np.amax(SKIRT_AvValues[i,:]) < 0.1:
        sphAvMask[i] = False # don't include galaxies with Av values that don't go above 0.1
    else:
        sphAvMask[i] = True 

plotAllAttenuationCurves()
plotAllAttenuationCurvesNorm()
plotEnergyBalance()
plotAvAxisRatio()
plotAvStellarMass()
plotAvsSFR()
plotAvDustToStellar()

end = timer()
time = end - start
time = str(datetime.timedelta(seconds=time))
print('Time for first analysis:', time)  

# Color-color Plots
print('starting color color plots')
os.system('mkdir -p '+plotPath+'colorPlots/')
for i in range(len(ratio_indicies)): # number of plots to make 
    colorColorPlots(i)

# axis ratio vs. color plots
print('starting axis ratio color plots')
os.system('mkdir -p '+plotPath+'fitPlots/')
os.system('mkdir -p '+plotPath+'axisRatioColorPlots/')
os.system('mkdir -p '+noDustPlotPath+'sSFRPlots/')
for i in range(len(axisRatio_color_indices)):
    #sSFRnoDust(i)
    #plotEachAxisRatioVsColor(i) # one plot for each NIHAO galaxy 
    dpAvg, sphAvg = fitBysSFR(i)
    #dpAvg, sphAvg = fitByDustToStellar(i)
    #dpAvg, sphAvg = fitByStellar(i)
    #dpAvg, sphAvg = fitByFluxRatios(i) # specify with fitFluxes 
    # calculate dp errors
    x = dp_flux[:,axisRatio_color_indices[i][0]]
    y = dp_flux[:,axisRatio_color_indices[i][1]]
    x_err = dp_err[:,axisRatio_color_indices[i][0]]
    y_err = dp_err[:,axisRatio_color_indices[i][1]]
    yerr = np.sqrt((x_err**2 * (1/(x * np.log(10))**2)) + (y_err**2 * (1/(y * np.log(10)))**2))
    colors = np.empty(dp_num, dtype=str)
    plt.figure(figsize=(10,8))   
    #plt.scatter(SKIRT_axisRatio, 
    #            np.log10(SKIRT_flux[:,:,axisRatio_color_indices[i][0]] / 
    #            SKIRT_flux[:,:,axisRatio_color_indices[i][1]]) - sphAvg, 
    #            marker='o', s=20, zorder=10, alpha=0.7, c='blue')
    plt.scatter(SKIRT_axisRatio[sphMassMask,:], 
                np.log10(SKIRT_flux[sphMassMask,:,axisRatio_color_indices[i][0]] / 
                SKIRT_flux[sphMassMask,:,axisRatio_color_indices[i][1]]) - sphAvg, 
                marker='o', s=20, zorder=10, alpha=0.7, c='blue')
    good_dp_x = np.asarray([]) # includes only disks with reliable photometry
    good_dp_y = np.asarray([])
    dpFitMask = np.asarray([], dtype=bool) # don't include red points in fit
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
            dpFitMask = np.append(dpFitMask, False)
        elif dp_bool[d,axisRatio_color_indices[i][1]]: # y is a lower limit (denominator is upper limit)
            lolims = True
            colors[d] = 'red'
            dpFitMask = np.append(dpFitMask, False)
        else:
            dpFitMask = np.append(dpFitMask, True)
        if bad:
            print('bad DustPedia data point')
        else:
            good_dp_x = np.append(good_dp_x, dp_axisRatio[d])
            good_dp_y = np.append(good_dp_y, np.log10(dp_flux[d,axisRatio_color_indices[i][0]] / 
                        dp_flux[d,axisRatio_color_indices[i][1]]) - dpAvg[d])
            plt.errorbar(good_dp_x[-1], good_dp_y[-1], 
                         xerr=0, yerr=yerr[d], elinewidth=0.2, marker='o',
                         markersize=5, linewidth=0, color=colors[d], zorder=0, alpha=0.3,
                         xuplims=xuplims, xlolims=xlolims, uplims=uplims, lolims=lolims)
    # linear fits for dp galaxies
    nanMask = np.isfinite(good_dp_x[dpFitMask]) & np.isfinite(good_dp_y[dpFitMask])
    dpAxisRatioColorLinearFit = np.polyfit(good_dp_x[dpFitMask][nanMask], good_dp_y[dpFitMask][nanMask], 1)
    dpAxisRatioColorFit = np.poly1d(dpAxisRatioColorLinearFit) 
    # save fit coefficients as numpy arrays
    np.save(plotPath+'axisRatioColorPlots/dp_'+axisRatio_plot_names[i]+'_fit.npy', dpAxisRatioColorLinearFit)
    dp_x_fit = np.linspace(np.amin(good_dp_x[dpFitMask][nanMask]), np.amax(good_dp_x[dpFitMask][nanMask]), num=2)
    dp_y_fit = dpAxisRatioColorFit(dp_x_fit)
    plt.plot(dp_x_fit, dp_y_fit, color='black', linewidth=2)
    # linear fits for NIHAO galaxies
    #sphAxisRatioColorLinearFit = np.polyfit(SKIRT_axisRatio.flatten(), 
    #                             np.log10(SKIRT_flux[:,:,axisRatio_color_indices[i][0]] / 
    #                             SKIRT_flux[:,:,axisRatio_color_indices[i][1]]).flatten() - 
    #                             sphAvg.flatten(), 1)
    sphAxisRatioColorLinearFit = np.polyfit(SKIRT_axisRatio[sphMassMask,:].flatten(), 
                                 np.log10(SKIRT_flux[sphMassMask,:,axisRatio_color_indices[i][0]] / 
                                 SKIRT_flux[sphMassMask,:,axisRatio_color_indices[i][1]]).flatten() - 
                                 sphAvg.flatten(), 1)
    sphAxisRatioColorFit = np.poly1d(sphAxisRatioColorLinearFit) 
    # save fit coefficients as numpy arrays
    np.save(plotPath+'axisRatioColorPlots/sph_'+axisRatio_plot_names[i]+'_fit.npy', sphAxisRatioColorLinearFit)
    sph_x_fit = np.linspace(np.amin(good_dp_x[dpFitMask]), np.amax(good_dp_x[dpFitMask]), num=2)
    sph_y_fit = sphAxisRatioColorFit(sph_x_fit)
    plt.plot(sph_x_fit, sph_y_fit, color='blue', linewidth=2)
    # format and save 
    plt.xlabel('Axis Ratio', fontsize=16)
    plt.ylabel('log('+band_names[axisRatio_color_indices[i][0]]+' / '+
               band_names[axisRatio_color_indices[i][1]]+') - <log('+
               band_names[axisRatio_color_indices[i][0]]+' / '+
               band_names[axisRatio_color_indices[i][1]]+')>', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.savefig(plotPath+'axisRatioColorPlots/'+axisRatio_plot_names[i]+'.png',dpi=300)
    plt.close()
print('done')
