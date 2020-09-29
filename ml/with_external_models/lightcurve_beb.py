import sys
import numpy as np
from numpy import savetxt


import pastis
pastis.initialize()

from pastis.MCMC import priors
from pastis import ObjectBuilder as ob
import pastis.models as mod

# Append configfiles to searchpath
#sys.path.append('configfiles/')
#tuve que cambiar para que me lo tome
sys.path.append('../../examples/configfiles/')



# Read import dict
from example_BEB import input_dict

# Get priors from configuration file
priordict = priors.prior_constructor(input_dict, {})

#not_passed = True
#i = 0

# por probar, lista de tuplas 

out_file = open("beb-index.txt", "w")


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
                    
                
        # Instantiate binary and foreground star
        try:
            beb, star = ob.ObjectBuilder(input_dict)
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
    tc = beb.orbital_parameters.T0
    per = beb.orbital_parameters.P
    dt = 2.0 / (24.0 * 60.0)
    
    # Sampling array, centred in Tc, width P, sampling; 2 min)
    t = np.arange(tc - per/2, tc + per/2, dt) 
    
    lc = mod.PHOT.PASTIS_PHOT(t, 'Kepler', False, 0.0, 1, 0.0,
                             *[star, beb])
    
    simu_name = 'beb-simu-'+str(simu_number)+'.csv'
    out_file_line.append(simu_name)
    
    
    #  write simulation to csv file
    savetxt(simu_name, lc, delimiter=',')
    # write file
    out_file.write(str(out_file_line)+'\n')


out_file.close()
# ax = plt.gca()
# ax.plot(t, f)
# Next: write results to file
