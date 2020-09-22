import sys
import numpy as np
from matplotlib import pyplot as plt

import pastis
pastis.initialize()

from pastis.MCMC import priors
from pastis import ObjectBuilder as ob
import pastis.models as mod

# Append configfiles to searchpath
sys.path.append('configfiles/')

# Read import dict
from example_TRIPLE import input_dict

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
            
    # Instantiate binary and foreground star
    try:
        objs = ob.ObjectBuilder(input_dict)
        not_passed = False
    except pastis.exceptions.EvolTrackError:
        print('Try again ({})'.format(i+1))
        i+=1
        pass

# =============================================================================
# # Set zero redenning for all stars
# beb.star1.ebmv = 0.0
# beb.star2.ebmv = 0.0
# star.ebmv = 0.0 
# =============================================================================

# Print stellar properties
# =============================================================================
# beb.star1.R
# beb.star2.R
# beb.dist
# star.dist
# star.get_LC(t)
# beb.star1.logL
# beb.star1.L
# beb.star2.L
#
# =============================================================================
# t = np.linspace(beb.orbital_parameters.T0 - beb.orbital_parameters.P,
#                 brb.orbital_parameters.T0 + beb.orbital_parameters.P, 1000)

tc = objs[0].object2.orbital_parameters.T0
per = objs[0].object2.orbital_parameters.P
dt = 2.0 / (24.0 * 60.0)

# Sampling array, centred in Tc, width P, sampling; 2 min)
t = np.arange(tc - per/2, tc + per/2, dt) 

lc = mod.PHOT.PASTIS_PHOT(t, 'Kepler', False, 0.0, 1, 0.0,
                         *objs)

ax = plt.gca()
ax.plot(t, f)
# Next: write results to file
