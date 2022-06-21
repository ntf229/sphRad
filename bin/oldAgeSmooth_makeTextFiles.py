# Creates text files for importing SPH data into SKIRT

import pynbody
import numpy as np
import math
from timeit import default_timer as timer
from os.path import expanduser
import os
import argparse
from scipy import spatial
from scipy.spatial.transform import Rotation as R

class galaxy:
    def __init__(self,name):
        filePath = '/scratch/ntf229/NIHAO/'+name+'/'+name+'.01024'
        
        # Load NIHAO data
        self.data = pynbody.load(filePath)
        
        self.full_x_pos = np.float32(self.data.star['pos'].in_units('pc')[:][:,0])
        self.full_y_pos = np.float32(self.data.star['pos'].in_units('pc')[:][:,1])
        self.full_z_pos = np.float32(self.data.star['pos'].in_units('pc')[:][:,2])
        self.full_mass = np.float32(self.data.star['massform'].in_units('Msol')[:]) # in solar masses	
        self.full_current_mass = np.float32(self.data.star['mass'].in_units('Msol')[:])
        self.full_metals = np.float32(self.data.star['metals'][:])
        self.full_age = np.float32(self.data.star['age'].in_units('yr')[:])
        self.full_x_vel = np.float32(self.data.star['vel'].in_units('km s**-1')[:][:,0])
        self.full_y_vel = np.float32(self.data.star['vel'].in_units('km s**-1')[:][:,1])
        self.full_z_vel = np.float32(self.data.star['vel'].in_units('km s**-1')[:][:,2])
        self.full_smooth = 2*np.float32(self.data.star['smooth'].in_units('pc')[:]) # 2 times gravitational softening length
        
        self.full_x_pos_dust = np.float32(self.data.gas['pos'].in_units('pc')[:][:,0])
        self.full_y_pos_dust = np.float32(self.data.gas['pos'].in_units('pc')[:][:,1])
        self.full_z_pos_dust = np.float32(self.data.gas['pos'].in_units('pc')[:][:,2])
        self.full_smooth_dust = 2*np.float32(self.data.gas['smooth'].in_units('pc')[:]) # 2 times gravitational softening length
        self.full_mass_dust = np.float32(self.data.gas['mass'].in_units('Msol')[:]) # in solar masses	
        self.full_metals_dust = np.float32(self.data.gas['metals'][:])
        self.full_temp_dust = np.float32(self.data.gas['temp'][:])
        self.full_density_dust = np.float32(self.data.gas['rho'].in_units('Msol pc**-3')[:])

        print('gas density units:',self.data.gas['rho'].units)
        print('max stellar particle mass:',np.amax(self.full_mass))
        print('min stellar particle mass:',np.amin(self.full_mass))
        
        self.full_length_star = len(self.full_x_pos) # starting length
        print(self.full_length_star, 'full star') 
        
        self.full_length_dust = len(self.full_x_pos_dust) # starting length
        print(self.full_length_dust, 'full dust') 
                
        # Halo catalogue
        h = self.data.halos() # ordered by number of particles (starts at 1)
        
        #vir_radius = h[1].properties['Rvir']
        
        if name == 'g3.49e11':
            haloNum = int(2) # first halo is off center for this galaxy
        else:
            haloNum = int(1)
        
        xMin = np.amin(h[haloNum]['x'].in_units('pc'))
        xMax = np.amax(h[haloNum]['x'].in_units('pc'))
        yMin = np.amin(h[haloNum]['y'].in_units('pc'))
        yMax = np.amax(h[haloNum]['y'].in_units('pc'))
        zMin = np.amin(h[haloNum]['z'].in_units('pc'))
        zMax = np.amax(h[haloNum]['z'].in_units('pc'))
        
        xLength = abs(xMax - xMin)
        yLength = abs(yMax - yMin)
        zLength = abs(zMax - zMin)
        
        diameter = np.amax([xLength,yLength,zLength]) / 6 # division by 6 here scales all galaxies (calibrated from images)
        
        # Custom changes based on optical (zrg) images
        if name == 'g1.05e11':
            diameter = diameter * 1.5
        elif name == 'g1.08e11':
            diameter = diameter * 1.5
        elif name == 'g1.52e11':
            diameter = diameter * 2
        elif name == 'g1.57e11':
            diameter = diameter * 1.5
        elif name == 'g1.59e11':
            diameter = diameter * 1.5
        elif name == 'g1.64e11':
            diameter = diameter * 1.5
        elif name == 'g3.19e10':
        	diameter = diameter * 4
        elif name == 'g3.23e11':
            diameter = diameter * 1.5
        elif name == 'g3.44e10':
            diameter = diameter * 1.5
        elif name == 'g3.49e11':
            diameter = diameter * 1.5
        elif name == 'g3.55e11':
            diameter = diameter * 1.5
        elif name == 'g3.59e11':
            diameter = diameter * 1.5
        elif name == 'g3.93e10':
            diameter = diameter * 2
        elif name == 'g6.37e10': 
            diameter = diameter * 2
        elif name == 'g2.83e10':
            diameter = diameter * 1.5
        elif name == 'g2.64e10':
            diameter = diameter * 2
        elif name == 'g2.80e10':
            diameter = diameter * 2
        elif name == 'g2.19e11':
            diameter = diameter * 2
        elif name == 'g2.34e10':
            diameter = diameter * 1.5
        elif name == 'g2.94e10':
            diameter = diameter * 1.5
        elif name == 'g3.71e11':
            diameter = diameter * 2
        elif name == 'g4.48e10':
            diameter = diameter * 2
        elif name == 'g4.27e10':
            diameter = diameter * 2
        elif name == 'g4.90e11':
            diameter = diameter * 1.5
        elif name == 'g4.94e10':
            diameter = diameter * 1.5
        elif name == 'g5.02e11':
            diameter = diameter * 1.5
        elif name == 'g5.05e10':
            diameter = diameter * 2
        elif name == 'g6.12e10':
            diameter = diameter * 1.5
        elif name == 'g6.37e10':
            diameter = diameter * 2
        elif name == 'g6.91e10':
            diameter = diameter * 1.5
        elif name == 'g6.77e10':
            diameter = diameter * 1.5
        elif name == 'g7.44e11':
            diameter = diameter * 1.5
        elif name == 'g8.06e11':
            diameter = diameter * 2
        elif name == 'g8.89e10':
            diameter = diameter * 1.5
        elif name == 'g9.59e10':
            diameter = diameter * 1.5
        elif name == 'g1.12e12':
            diameter = diameter * 0.75
        elif name == 'g1.92e12':
            diameter = diameter * 0.75
        elif name == 'g7.66e11':
            diameter = diameter * 0.75 
        elif name == 'g8.13e11':
            diameter = diameter * 0.75      
        
        xCenter = (xMax + xMin)/2
        yCenter = (yMax + yMin)/2
        zCenter = (zMax + zMin)/2
        
        # Make cuts based on position
        self.xAbsMin = xCenter - (diameter/2)															
        self.xAbsMax = xCenter + (diameter/2)
        self.yAbsMin = yCenter - (diameter/2)
        self.yAbsMax = yCenter + (diameter/2)
        self.zAbsMin = zCenter - (diameter/2)
        self.zAbsMax = zCenter + (diameter/2)
        
    def starCut(self):
        xMin = self.full_x_pos > self.xAbsMin
        xMax = self.full_x_pos[xMin] < self.xAbsMax
        yMin = self.full_y_pos[xMin][xMax] > self.yAbsMin
        yMax = self.full_y_pos[xMin][xMax][yMin] < self.yAbsMax
        zMin = self.full_z_pos[xMin][xMax][yMin][yMax] > self.zAbsMin
        zMax = self.full_z_pos[xMin][xMax][yMin][yMax][zMin] < self.zAbsMax
        
        self.x_pos = self.full_x_pos[xMin][xMax][yMin][yMax][zMin][zMax]
        self.y_pos = self.full_y_pos[xMin][xMax][yMin][yMax][zMin][zMax]
        self.z_pos = self.full_z_pos[xMin][xMax][yMin][yMax][zMin][zMax]
        self.smooth = self.full_smooth[xMin][xMax][yMin][yMax][zMin][zMax]
        self.x_vel = self.full_x_vel[xMin][xMax][yMin][yMax][zMin][zMax]
        self.y_vel = self.full_y_vel[xMin][xMax][yMin][yMax][zMin][zMax]
        self.z_vel = self.full_z_vel[xMin][xMax][yMin][yMax][zMin][zMax]
        self.mass = self.full_mass[xMin][xMax][yMin][yMax][zMin][zMax]
        self.current_mass = self.full_current_mass[xMin][xMax][yMin][yMax][zMin][zMax]
        self.metals = self.full_metals[xMin][xMax][yMin][yMax][zMin][zMax]
        self.age = self.full_age[xMin][xMax][yMin][yMax][zMin][zMax]
        
        self.starLength = len(self.x_pos)
        print(self.starLength, 'stars')
        
    def dustCut(self):    
        xMin = self.full_x_pos_dust > self.xAbsMin
        xMax = self.full_x_pos_dust[xMin] < self.xAbsMax
        yMin = self.full_y_pos_dust[xMin][xMax] > self.yAbsMin
        yMax = self.full_y_pos_dust[xMin][xMax][yMin] < self.yAbsMax
        zMin = self.full_z_pos_dust[xMin][xMax][yMin][yMax] > self.zAbsMin
        zMax = self.full_z_pos_dust[xMin][xMax][yMin][yMax][zMin] < self.zAbsMax
        
        self.x_pos_dust = self.full_x_pos_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        self.y_pos_dust = self.full_y_pos_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        self.z_pos_dust = self.full_z_pos_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        self.smooth_dust = self.full_smooth_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        self.mass_dust = self.full_mass_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        self.metals_dust = self.full_metals_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        self.temp_dust = self.full_temp_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        self.density_dust = self.full_density_dust[xMin][xMax][yMin][yMax][zMin][zMax]
        
        self.dustLength = len(self.x_pos_dust)
        print(self.dustLength, 'dust')
        
    def shift(self):    
        xCenter = (self.xAbsMax + self.xAbsMin)/2
        yCenter = (self.yAbsMax + self.yAbsMin)/2
        zCenter = (self.zAbsMax + self.zAbsMin)/2
        
        self.x_pos = self.x_pos - xCenter
        self.y_pos = self.y_pos - yCenter
        self.z_pos = self.z_pos - zCenter
        
        self.x_pos_dust = self.x_pos_dust - xCenter
        self.y_pos_dust = self.y_pos_dust - yCenter
        self.z_pos_dust = self.z_pos_dust - zCenter

    def youngStars(self, ageSmooth, SF, tauClear):
        tau_clear = float(tauClear) * 1.e6 # convert from Myrs to years
        if eval(ageSmooth):
            # star-formation lasts smooth_time years
            if  self.starLength > 2e5: # ~ 35 galaxies
                smooth_time = 0.
            elif self.starLength > 1e5:
                smooth_time = 1.e8
            elif self.starLength > 5e4:
                smooth_time = 2.e8
            elif self.starLength > 4e4:
                smooth_time = 3.e8
            elif self.starLength > 3e4:
                smooth_time = 4.e8
            elif self.starLength > 2e4:
                smooth_time = 5.e8
            else:
                smooth_time = 6.e8
        else:
            smooth_time = 0.
        k = 1.3806485*10**(-19) # boltzmann constant in cm**2 kg s**-2 K**-1
        youngStarMask = self.age < (smooth_time + 1.e7) # Mask for young stars (smooth_time + 10 Myrs)
        youngStarIndex = [] # Indices of young stars 
        for i in range(self.starLength):
            if youngStarMask[i]:
                youngStarIndex.append(i)
        print('number of young star particles before smoothing:',len(youngStarIndex))
        self.young_x_pos = []
        self.young_y_pos = []
        self.young_z_pos = []
        self.young_mass = []
        self.young_current_mass = []
        self.young_metals = []
        self.young_age = []
        self.young_x_vel = []
        self.young_y_vel = []
        self.young_z_vel = []
        self.young_smooth = []
        self.young_SFR = []
        self.young_logC = []
        self.young_p = []
        self.young_f_PDR = []
        def mass_prob(x): # star cluster mass probability distribution
            return x**(-1.8)
        mass_min = 700
        mass_max = 1e6
        prob_max = 7.57*10**(-6) # slightly larger than 700**(-1.8)
        N_MC = 10000000 # masked distribution is well resolved and still runs fast
        # accept / reject Monte Carlo:
        mass = np.random.uniform(mass_min,mass_max,N_MC)  # get uniform temporary mass values
        prob = np.random.uniform(0,prob_max,N_MC)  # get uniform random probability values
        mask = prob < mass_prob(mass) # accept / reject
        sampled_masses = mass[mask] # sample of star cluster masses following the desired distribution
        for i in range(len(youngStarIndex)):
            parent_index = youngStarIndex[i]
            parent_mass = self.mass[parent_index]
            parent_age = self.age[parent_index]
            # parent_mass_not_formed is the fraction of mass that hasn't formed yet
            if smooth_time == 0 and eval(SF):
                parent_mass_not_formed = 0. # no smoothing
                ind = len(self.young_x_pos) # don't subtract 1, haven't appended yet
                # create MAPPINGS-III particle
                self.young_x_pos.append(self.x_pos[parent_index])
                self.young_y_pos.append(self.y_pos[parent_index])
                self.young_z_pos.append(self.z_pos[parent_index])
                self.young_current_mass.append(self.mass[parent_index]) # assume no mass loss 
                self.young_metals.append(self.metals[parent_index])
                self.young_age.append(self.age[parent_index])
                self.young_x_vel.append(self.x_vel[parent_index])
                self.young_y_vel.append(self.y_vel[parent_index])
                self.young_z_vel.append(self.z_vel[parent_index])
                self.young_smooth.append(self.smooth[parent_index])
                self.young_f_PDR.append(np.exp(-1*parent_age / tau_clear)) # Groves 2008 equation 16
                self.young_mass.append(parent_mass * (1. - parent_mass_not_formed))
                self.young_SFR.append(1.e-7 * self.young_mass[ind]) # (units: yr**-1) assumes constant SFR over the last 10 Myrs
                self.young_logC.append(np.random.normal(5., 0.4))
                self.young_p.append(k*10**((5/2)*(self.young_logC[ind]-(3/5)*np.log10(np.random.choice(sampled_masses))))) # in kg cm**-1 s**-2
                self.young_p[ind] *= 100 # convert to Pascals
                # create ghost particle
                self.x_pos_dust = np.append(self.x_pos_dust, self.x_pos[parent_index])
                self.y_pos_dust = np.append(self.y_pos_dust, self.y_pos[parent_index])
                self.z_pos_dust = np.append(self.z_pos_dust, self.z_pos[parent_index])
                self.smooth_dust = np.append(self.smooth_dust, 3*self.smooth[parent_index]) # 3 times smoothing length of parent star
                self.mass_dust = np.append(self.mass_dust, -10.*self.young_mass[ind]*self.young_f_PDR[ind] + parent_mass*parent_mass_not_formed) # 10 times larger negative mass in PDR plus mass not yet formed 
                self.metals_dust = np.append(self.metals_dust, self.metals[parent_index]) 
                self.temp_dust = np.append(self.temp_dust, 7000.) # assume 7000K (doesn't make a difference as long as lower than maxTemp)
                self.density_dust = np.append(self.density_dust, 0.) # don't need density, set to 0
            elif smooth_time != 0: # some age smoothing
                if (parent_age > 1.e7) and (parent_age <= smooth_time): # 10 Myrs < age < smooth_time
                    parent_mass_not_formed = (smooth_time - parent_age) / smooth_time
                elif parent_age <= 1.e7: # younger than 10 Myrs
                    parent_mass_not_formed = (smooth_time - parent_age) / smooth_time
                else: # between smooth_time and smooth_time + 10 Myrs
                    parent_mass_not_formed = 0.
                # create particles every 1 Myrs 
                num = math.floor(parent_age*1.e-6) + 1
                new_mass = (parent_mass*(1-parent_mass_not_formed)) / num 
                ghost_mass = 0.
                for j in range(num):
                    # create new particle and track ghost mass
                    new_age = (parent_age*1.e-6 - j)*1.e6
                    if new_age < 1.e7 and eval(SF): # less than 10 Myrs
                        ind = len(self.young_x_pos) # don't subtract 1, haven't appended yet
                        # create MAPPINGS-III particle 
                        self.young_x_pos.append(self.x_pos[parent_index])
                        self.young_y_pos.append(self.y_pos[parent_index])
                        self.young_z_pos.append(self.z_pos[parent_index])
                        self.young_current_mass.append(new_mass) # assume no mass loss 
                        self.young_metals.append(self.metals[parent_index])
                        self.young_age.append(new_age)
                        self.young_x_vel.append(self.x_vel[parent_index])
                        self.young_y_vel.append(self.y_vel[parent_index])
                        self.young_z_vel.append(self.z_vel[parent_index])
                        self.young_smooth.append(self.smooth[parent_index])
                        self.young_f_PDR.append(np.exp(-1*new_age / tau_clear)) # Groves 2008 equation 16
                        self.young_mass.append(new_mass)
                        self.young_SFR.append(1.e-7 * self.young_mass[ind]) # (units: yr**-1) assumes constant SFR over the last 10 Myrs
                        self.young_logC.append(np.random.normal(5., 0.4))
                        self.young_p.append(k*10**((5/2)*(self.young_logC[ind]-(3/5)*np.log10(np.random.choice(sampled_masses))))) # in kg cm**-1 s**-2
                        self.young_p[ind] *= 100 # convert to Pascals
                        ghost_mass += -10.*self.young_mass[ind]*self.young_f_PDR[ind]
                    else:
                        # create star particle
                        self.x_pos = np.append(self.x_pos, self.x_pos[parent_index])
                        self.y_pos = np.append(self.y_pos, self.y_pos[parent_index])
                        self.z_pos = np.append(self.z_pos, self.z_pos[parent_index])
                        self.mass = np.append(self.mass, new_mass)
                        self.current_mass = np.append(self.current_mass, new_mass) # assume no mass loss
                        self.metals = np.append(self.metals, self.metals[parent_index])
                        self.age = np.append(self.age, new_age)
                        self.x_vel = np.append(self.x_vel, self.x_vel[parent_index])
                        self.y_vel = np.append(self.y_vel, self.y_vel[parent_index])
                        self.z_vel = np.append(self.z_vel, self.z_vel[parent_index])
                        self.smooth = np.append(self.smooth, self.smooth[parent_index])
                # create ghost particle (considering mass not yet formed)
                # if SF == False, these will have positive mass (mass not yet formed)
                self.x_pos_dust = np.append(self.x_pos_dust, self.x_pos[parent_index])
                self.y_pos_dust = np.append(self.y_pos_dust, self.y_pos[parent_index])
                self.z_pos_dust = np.append(self.z_pos_dust, self.z_pos[parent_index])
                self.smooth_dust = np.append(self.smooth_dust, 3*self.smooth[parent_index]) # 3 times smoothing length of parent star
                self.mass_dust = np.append(self.mass_dust, ghost_mass + parent_mass*parent_mass_not_formed) # 10 times larger negative mass in PDR plus mass not yet formed 
                self.metals_dust = np.append(self.metals_dust, self.metals[parent_index]) 
                self.temp_dust = np.append(self.temp_dust, 7000.) # assume 7000K (doesn't make a difference as long as lower than maxTemp)
                self.density_dust = np.append(self.density_dust, 0.) # don't need density, set to 0
        # change from lists to numpy arrays
        self.young_x_pos = np.asarray(self.young_x_pos)
        self.young_y_pos = np.asarray(self.young_y_pos)
        self.young_z_pos = np.asarray(self.young_z_pos)
        self.young_mass = np.asarray(self.young_mass)
        self.young_current_mass = np.asarray(self.young_current_mass)
        self.young_metals = np.asarray(self.young_metals)
        self.young_age = np.asarray(self.young_age)
        self.young_x_vel = np.asarray(self.young_x_vel)
        self.young_y_vel = np.asarray(self.young_y_vel)
        self.young_z_vel = np.asarray(self.young_z_vel)
        self.young_smooth = np.asarray(self.young_smooth)
        self.young_SFR = np.asarray(self.young_SFR)
        self.young_logC = np.asarray(self.young_logC)
        self.young_p = np.asarray(self.young_p)
        self.young_f_PDR = np.asarray(self.young_f_PDR)
        # delete parent particles
        if smooth_time == 0 and not eval(SF):
            print('Keeping parent particles')
        else:
            self.x_pos = np.delete(self.x_pos, youngStarIndex)
            self.y_pos = np.delete(self.y_pos, youngStarIndex)
            self.z_pos = np.delete(self.z_pos, youngStarIndex)
            self.mass = np.delete(self.mass, youngStarIndex)
            self.current_mass = np.delete(self.current_mass, youngStarIndex)
            self.metals = np.delete(self.metals, youngStarIndex)
            self.age = np.delete(self.age, youngStarIndex)
            self.x_vel = np.delete(self.x_vel, youngStarIndex)
            self.y_vel = np.delete(self.y_vel, youngStarIndex)
            self.z_vel = np.delete(self.z_vel, youngStarIndex)
            self.smooth = np.delete(self.smooth, youngStarIndex)
        
    def clumps(self, numCells, numClumps):
        # replace particles that enclose one or more SF regions by a set of smaller particles with the same total mass
        # use gridded phantom particles with negligable mass to force small grid cells
        if len(self.young_x_pos) == 0: # no extra clumpiness
            print('no young star particles so no clumps')
            self.ag_x = self.x_pos_dust
            self.ag_y = self.y_pos_dust
            self.ag_z = self.z_pos_dust
            self.ag_h = self.smooth_dust
            self.ag_M = self.mass_dust
            self.ag_Z = self.metals_dust
            self.ag_T = self.temp_dust
        else:    
            # pour star formation region positions into a search tree
            points = np.stack((self.young_x_pos,self.young_y_pos,self.young_z_pos)).T
            tree = spatial.KDTree(points)
            # new list of gas particles, each particle is a tuple
            aglist = []
            phantom_w = 1.e-6 # fraction of mass to use for phantom particles
            numCells = int(numCells)
            numClumps = int(numClumps)
            clump_w = (1 - (phantom_w * (numCells**3 - numClumps)) ) / numClumps # fraction of mass to use for clump particles
            grid = np.zeros((numCells**3, 3)) # each cell gets 3 values (x,y,z) indicies (0 - numCells-1)
            count = 0
            stepSizeUnscaled = 1./numCells # total range of grid is 1
            left_edge = -1*stepSizeUnscaled*numCells/2
            edges_1d = np.linspace(left_edge, left_edge + (numCells-1)*stepSizeUnscaled, num=numCells)
            for i in range(numCells):
            	for j in range(numCells):
            		for k in range(numCells):
            		    # position of left bin (centered at 0)
            			grid[count, 0] = edges_1d[i]
            			grid[count, 1] = edges_1d[j]
            			grid[count, 2] = edges_1d[k]
            			count += 1
            for x,y,z,h,M,Z,T in zip(self.x_pos_dust,self.y_pos_dust,self.z_pos_dust,self.smooth_dust,self.mass_dust,self.metals_dust,self.temp_dust):
                neighbors = tree.query_ball_point((x,y,z), r=h)
                if M > 0 and len(neighbors) > 0: # no subgrid clumping for negative mass ghost particles
                    stepSize = stepSizeUnscaled*h
                    self.young_smooth[neighbors] = stepSize # change smoothing length of neighboring young stars to fit within grid cells
                    grid_centers = h*grid + stepSize/2
                    # apply random rotations
                    rot = R.random()
                    grid_centers = rot.apply(grid_centers)
                    # shift positions
                    grid_centers[:,0] += x
                    grid_centers[:,1] += y
                    grid_centers[:,2] += z
                    # create array which only stores available grid cells 
                    grid_a = grid_centers
                    # find index corresponding to position of young stars
                    for n in range(len(neighbors)):
                        minDist = 1e9 # start with very high number (1 Mpc)
                        for i in range(numCells**3 - n): # minus n because each n loop deletes a particle
                            dist = np.sqrt( (self.young_x_pos[neighbors[n]]-grid_a[i,0])**2 +
                                            (self.young_y_pos[neighbors[n]]-grid_a[i,1])**2 + 
                                            (self.young_z_pos[neighbors[n]]-grid_a[i,2])**2 )
                            if dist < minDist:
                                minDist = dist
                                index = i
                        # put phantom particle in center of cell containing young star, remove cell from list
                        aglist.append((grid_a[index,0], grid_a[index,1], grid_a[index,2], stepSize, M*phantom_w, Z, T))
                        grid_a = np.delete(grid_a, index, axis=0) 
                    # sample clump positions (no repeats)
                    for i in range(numClumps):
                    	grid_len = len(grid_a[:,0])
                    	clumpIndex = int(np.random.choice(np.linspace(0, grid_len-1, num=grid_len)))
                    	aglist.append((grid_a[clumpIndex,0], grid_a[clumpIndex,1], grid_a[clumpIndex,2], stepSize, M*clump_w, Z, T))
                    	grid_a = np.delete(grid_a, clumpIndex, axis=0)
                    # phantom particles
                    for i in range(len(grid_a[:,0])):
                        aglist.append((grid_a[i,0], grid_a[i,1], grid_a[i,2], stepSize, M*phantom_w, Z, T))
                else:
                    aglist.append((x,y,z,h,M,Z,T))
            # convert list to arrays
            self.ag_x, self.ag_y, self.ag_z, self.ag_h, self.ag_M, self.ag_Z, self.ag_T = np.array(aglist).T

if __name__=='__main__':
    
    start = timer()
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--ageSmooth") # if True, smooth ages based on number of star particles (see galaxy.youngStars()) 
    parser.add_argument("--SF") # if True, star particles younger than 10 Myrs are assigned MAPPINGS-III SEDs
    parser.add_argument("--tauClear") # clearing time in Myrs for MAPPINGS-III f_PDR calculations (only matters if SF=True)
    parser.add_argument("--clumps") # if True, add subgrid clumpiness to gas near MAPPINGS-III particles
    parser.add_argument("--numCells") # number of cells along one dimension for subgrid clumping (only matters if clumps=True)
    parser.add_argument("--numClumps") # number of clumps for subgrid clumping (only matters if clumps=True)
    parser.add_argument("--galaxy") # name of galaxy 
    args = parser.parse_args()

    # Directory structure stores important parameters
    textPath = '/scratch/ntf229/sphRad/resources/NIHAO/TextFiles/'
    if eval(args.ageSmooth):
        textPath += 'ageSmooth/'
    else:
        textPath += 'noAgeSmooth/'
    if eval(args.SF):
        textPath += 'SF/tauClear'+args.tauClear+'/'
    else:
        textPath += 'noSF/'
    if eval(args.clumps):
        textPath += 'clumps/numCells'+args.numCells+'/numClumps'+args.numClumps+'/'
    else:
        textPath += 'noClumps/'

    textPath += args.galaxy+'/'

    os.system('mkdir -p '+textPath)

    star_header = 'Column 1: position x (pc)\nColumn 2: position y (pc)\nColumn 3: position z (pc)\nColumn 4: smoothing length (pc)\nColumn 5: v_x (km/s)\nColumn 6: v_y (km/s)\nColumn 7: v_z (km/s)\nColumn 8: mass (Msun)\nColumn 9: metallicity ()\nColumn 10: age (yr)'
    gas_header = 'Column 1: position x (pc)\nColumn 2: position y (pc)\nColumn 3: position z (pc)\nColumn 4: smoothing length (pc)\nColumn 5: mass (Msun)\nColumn 6: metallicity ()\nColumn 7: temperature (K)'
        
    
    print('starting ', args.galaxy)
    g = galaxy(args.galaxy)
    
    g.starCut()
    g.dustCut()
    g.shift()
    print('stellar mass before MAPPINGSIII:', np.sum(g.mass))
    g.youngStars(args.ageSmooth, args.SF, args.tauClear)
    print('FSPS stellar mass:', np.sum(g.mass))
    print('MAPPINGSIII stellar mass:', np.sum(g.young_mass))
    if eval(args.clumps):
        g.clumps(args.numCells, args.numClumps)
        print('gas mass after clumps:', np.sum(g.ag_M))
        print('fractional increase in number of gas particles from clumps:', float(len(g.ag_x)) / float(len(g.x_pos_dust)))
        np.savetxt(textPath+'gas.txt',np.float32(np.c_[g.ag_x, g.ag_y, g.ag_z, g.ag_h, g.ag_M, g.ag_Z, g.ag_T]),header=gas_header)
    else:
        np.savetxt(textPath+'gas.txt',np.float32(np.c_[g.x_pos_dust, g.y_pos_dust, g.z_pos_dust, g.smooth_dust, g.mass_dust, g.metals_dust, g.temp_dust]),header=gas_header)
    
    np.savetxt(textPath+'stars.txt',np.float32(np.c_[g.x_pos, g.y_pos, g.z_pos, g.smooth, g.x_vel, g.y_vel, g.z_vel, g.mass, g.metals, g.age]),header=star_header)
    if eval(args.SF):
        np.savetxt(textPath+'youngStars.txt',np.float32(np.c_[g.young_x_pos, g.young_y_pos, g.young_z_pos, g.young_smooth, g.young_x_vel, g.young_y_vel, g.young_z_vel, g.young_SFR, g.young_metals, g.young_logC, g.young_p, g.young_f_PDR]))
    
    np.savetxt(textPath+'current_mass_stars.txt',np.float32(np.c_[g.current_mass]))
    np.savetxt(textPath+'gas_density.txt',np.float32(np.c_[g.density_dust]))
    
    end = timer()
    print('time: ', end - start)

