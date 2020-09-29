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
from example_PIB import input_dict

# Get priors from configuration file
priordict = priors.prior_constructor(input_dict, {})

not_passed = True
i = 0

# Esta parte hay que repetirla N veces
out_file = open("test.txt","w") 
line=[]
while not_passed:
    # Randomly draw from prior a value for each parameter with
    # flag > 0.
    for obj in input_dict:
        pd = input_dict[obj]
        for par in pd:
            if isinstance(pd[par], list) and pd[par][1] > 0:
                pd[par][0] = priordict[obj+'_'+par].rvs()
                line.append((obj+'_'+par,pd[par][0]))
            
    # Instantiate binary and foreground star
    try:
        objs = ob.ObjectBuilder(input_dict)
        not_passed = False
    except pastis.exceptions.EvolTrackError:
        print('Try again ({})'.format(i+1))
        i+=1
        pass

t = np.linspace(-0.2, 1.2, 1000)
f = mod.PHOT.PASTIS_PHOT(t, 'Kepler', True, 0.0, 1, 0.0,
                         *objs)
ax = plt.gca()
ax.plot(t, f)
# Next: write results to file

out_file.write(str(line))    
print(line)

