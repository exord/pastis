import os
import sys
import numpy as np
# from matplotlib import pyplot as plt

import pastis
from pastis.MCMC import priors
from pastis import ObjectBuilder as ob
import pastis.models as mod

import toi_dist as td

# Prepare TOI dist
home = os.getenv('HOME')
csvdir = os.path.join(home, 'ExP/pastisML')
csvfile = 'csv-file-toi-catalog.csv'
csvfullpath = os.path.join(csvdir, csvfile)

toidist = td.prepare_toi_dist(csvfullpath)

# Append configfiles to searchpath
sys.path.append('examples/configfiles/')

# Read import dict
from example_BTP import input_dict

# Get priors from configuration file
priordict = priors.prior_constructor(input_dict, {})

not_passed = True
i = 0
while not_passed:
    # Randomly draw from prior a value for each parameter with
    # flag > 0.
    for obj in input_dict:
        pd = input_dict[obj]
        for par in pd:
            if isinstance(pd[par], list) and pd[par][1] > 0:
                pd[par][0] = priordict[obj+'_'+par].rvs()

    # Sample from TOI list
    lp, ld, ldur = toidist.resample(1)[:, 0]
    depth = 10**ld * 1e-6 # originally in ppm
    
    # Fix period
    input_dict['FitPlanet1']['P'][0] = 10**lp
    
    # Foreground flux based on depth
    kr = input_dict['FitPlanet1']['kr'][0]
    
    if kr**2 < depth:
        print('Cannot dilute depth; already too shallow')
        i+=1
        continue
    
    input_dict['FitPlanet1']['foreflux'][0] = kr**2/depth
    
    ## TODO: include duration info    
            
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