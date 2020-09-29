import os
import sys
import numpy as np

import pastis
# pastis.initialize()

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
# sys.path.append('../examples/configfiles/')
sys.path.append('../../examples/configfiles/')

# Read import dict
from example_PIB import input_dict

# Get priors from configuration file
priordict = priors.prior_constructor(input_dict, {})

not_passed = True
i = 0

# Esta parte hay que repetirla N veces
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
    input_dict['FitBinary1']['P'][0] = 10**lp
    
    # Foreground flux based on depth
    kr = input_dict['FitPlanet1']['kr'][0]
    
    if kr**2 < depth:
        print('Cannot dilute depth; already too shallow')
        i+=1
        break
    
    input_dict['FitPlanet1']['foreflux'][0] = kr**2/depth
    
    # Instantiate binary and foreground star
    try:
        objs = ob.ObjectBuilder(input_dict)
        not_passed = False
    except pastis.exceptions.EvolTrackError:
        print('Try again ({})'.format(i+1))
        i+=1
        pass

tc = objs[0].object2.planets[0].orbital_parameters.T0
per = objs[0].object2.planets[0].orbital_parameters.P
dt = 2.0 / (24.0 * 60.0)

# Sampling array, centred in Tc, width P, sampling; 2 min)
t = np.arange(tc - per/2, tc + per/2, dt) 

lc = mod.PHOT.PASTIS_PHOT(t, 'Kepler', False, 0.0, 1, 0.0,
                         *objs)
