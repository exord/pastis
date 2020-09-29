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

# >>> CORRER SOLO UNA VEZ

# Prepare TOI dist
home = os.getenv('HOME')
csvdir = os.path.join(home, 'rocky/pastis/ml')
csvfile = 'csv-file-toi-catalog.csv'
csvfullpath = os.path.join(csvdir, csvfile)

toidist = td.TOIdist(csvfullpath)

# Append configfiles to searchpath
sys.path.append('../examples/configfiles/')

# Read import dict
from example_BEBfit import input_dict

# Get priors from configuration file
priordict = priors.prior_constructor(input_dict, {})

out_file = open("./simulations/beb-index.txt", "w")
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
        depth = 10**ld * 1e-6 # originally in ppm
        
        # Fix period
        input_dict['FitBinary1']['P'][0] = 10**lp
        
        # Foreground flux based on depth
        kr = input_dict['FitBinary1']['kr'][0]
        input_dict['FitBinary1']['foreflux'][0] = kr**2/depth
        
        ## TODO: include duration info    
                
        # Instantiate binary and foreground star
        try:
            objs = ob.ObjectBuilder(input_dict)
            not_passed = False
        except pastis.exceptions.EvolTrackError:
            print('Objects could not be built.\n'
                  'Try again ({})'.format(i+1))
            i+=1
            break
         
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
            break
    
    # Dilute planet LC using element in input_dict
    f3 = input_dict['FitBinary1']['foreflux'][0]
    lc = (f/f3 + 1)/(1 + 1/f3)
    
    simu_name = './simulations/beb-simu-'+str(simu_number)+'.csv'
    out_file_line.append(simu_name)
    
    #  write simulation to csv file
    savetxt(simu_name, lc, delimiter=',')
    # write file
    out_file.write(str(out_file_line)+'\n')

out_file.close()
