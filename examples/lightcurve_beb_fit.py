import sys
import numpy as np
from matplotlib import pyplot as plt

import pastis
from pastis.MCMC import priors
from pastis import ObjectBuilder as ob
import pastis.models as mod

# Append configfiles to searchpath
sys.path.append('configfiles/')

# Read import dict
from example_BEBfit import input_dict

# Get priors from configuration file
priordict = priors.prior_constructor(input_dict, {})

# Time array
t = np.linspace(-0.2, 1.2, 1000)

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
        print('Objects could not be built.\n'
              'Try again ({})'.format(i+1))
        i+=1
        break
     
    try:
        f = mod.PHOT.PASTIS_PHOT(t, 'Kepler', True, 0.0, 1, 0.0,
                                 *objs)
        not_passed = False
    except pastis.exceptions.EBOPparamError:
        print('Light curve could not be computed.\n'
              'Try again ({})'.format(i+1))
        i+=1
        break

# Dilute binary LC using element in input_dict
f3 = input_dict['FitBinary1']['foreflux'][0]
lc = (f/f3 + 1)/(1 + 1/f3)

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

ax = plt.gca()
ax.plot(t, lc)
# Next: write results to file
