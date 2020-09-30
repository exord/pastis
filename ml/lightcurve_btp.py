import os
import sys
import numpy as np
from numpy import savetxt

# from matplotlib import pyplot as plt

import pastis
from pastis.MCMC import priors
from pastis import ObjectBuilder as ob
import pastis.models as mod

import toi_dist as td
import tools

# >>> CORRER SOLO UNA VEZ

# Prepare TOI dist
home = os.getenv('HOME')
# csvdir = os.path.join(home, 'ExP/pastisML')
csvdir = os.path.join(home, 'rocky/pastis/ml')
csvfile = 'csv-file-toi-catalog.csv'
csvfullpath = os.path.join(csvdir, csvfile)

toidist = td.TOIdist(csvfullpath)

# Append configfiles to searchpath
sys.path.append('../examples/configfiles/')

# Read import dict
from example_BTP import input_dict

# Get priors from configuration file
priordict = priors.prior_constructor(input_dict, {})

out_file = open("./simulations/btp-index.txt", "w")
# <<<

for simu_number in range(5): #las n veces
    not_passed = True
    i = 0
    out_file_line=[]
    
    while not_passed:
        # Randomly draw from prior a value for each parameter with
        # flag > 0.
        for obj in input_dict:
            pd = input_dict[obj]
            for par in pd:
                if isinstance(pd[par], list) and pd[par][1] > 0:
                    pd[par][0] = priordict[obj+'_'+par].rvs()
                    out_file_line.append((par,pd[par][0]))
    
        # Sample from TOI list
        (lp, ld, ldur), toidict = toidist.sample(1)
        depth = 10**ld[0] * 1e-6 # originally in ppm
        
        # Fix period
        input_dict['FitPlanet1']['P'][0] = 10**lp[0]
        
        # Foreground flux based on depth
        
# =============================================================================
#         # Option 1; fix foreflux based on kr
#         kr = input_dict['FitPlanet1']['kr'][0]
#         
#         if kr**2 < depth:
#             print('Cannot dilute depth; already too shallow')
#             i+=1
#             continue
#         
#         input_dict['FitPlanet1']['foreflux'][0] = kr**2/depth
#         
# =============================================================================

        # Option 2; fix kr based on foreflux        
        foreflux = input_dict['FitPlanet1']['foreflux'][0] 
        
        kr2 = foreflux * depth
        
        if kr2 > 0.1:
            print('Unreallistically large kr')
            i+=1
            continue
        
        input_dict['FitPlanet1']['kr'][0] = np.sqrt(kr2)
        
        # set a/R* based on transit duration
        ecc = input_dict['FitPlanet1']['ecc'][0]
        omega_deg = input_dict['FitPlanet1']['omega'][0]                         
        b = input_dict['FitPlanet1']['b'][0]
        kr = input_dict['FitPlanet1']['kr'][0]
                          
        ar = tools.invert_duration(10**ldur[0], 10**lp[0], ecc, 
                                   omega_deg, kr, b)
        input_dict['FitPlanet1']['ar'][0] = ar
                
        # Instantiate binary and foreground star
        try:
            objs = ob.ObjectBuilder(input_dict)
            not_passed = False
        except pastis.exceptions.EvolTrackError:
            print('Objects could not be built.\n'
                  'Try again ({})'.format(i+1))
            i+=1
            continue
         
        try:
            tc = objs[0].orbital_parameters.T0
            per = objs[0].orbital_parameters.P
            dt = 2.0 / (24.0 * 60.0)
            
            # Sampling array, centred in Tc, width P, sampling; 2 min)
            t = np.arange(tc - per/2, tc + per/2, dt) 
    
            f = mod.PHOT.PASTIS_PHOT(t, 'Kepler', False, 0.0, 1, 0.0,
                                     *objs)
            not_passed = False
        except pastis.exceptions.EBOPparamError:
            print('Light curve could not be computed.\n'
                  'Try again ({})'.format(i+1))
            i+=1
            continue
    
    # Dilute binary LC using element in input_dict
    f3 = input_dict['FitPlanet1']['foreflux'][0]
    lc = (f/f3 + 1)/(1 + 1/f3)
    
    #save simulation and values            
    simu_name = './simulations/btp-simu-'+str(simu_number)+'.csv'   
    savetxt(simu_name, lc, delimiter=',') #as np array

    for tuple in out_file_line:
        out_file.write(str(tuple[0]) + " "+ str(tuple[1]) + ",")
    
    out_file.write(simu_name + "\n")
    

out_file.close()
