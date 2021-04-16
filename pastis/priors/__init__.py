import os
import pickle
import numpy as np
from math import pi, log10
from scipy import interpolate
import scipy.stats as st

# import ephem
from .. import *

## Absorption ratio in a given band, relative to visual V (for Rv = 3.1)
## TAken from table 21.6 of Allen's astrophysical quantities.
A_Av = {'L' : 0.051,
        'K' : 0.108,
        'H' : 0.176,
        'J' : 0.282,
        'I' : 0.479,
        'R' : 0.749,
        'V' : 1.00,
        'B' : 1.31,
        'U' : 1.56}
        
def run_prior(pastisfile):
    """
    Run computation of hypotheses priors.
    """

    print("Welcome to the prior calculation model")

    initialize_priors(pastisfile)

    import moduleprior as mp
    import download_model as dm
    
    f = open(pastisfile, 'r')
    dd = pickle.load(f)
    f.close()

    infodict = dd[0]
    datadict = dd[2]
    hyppriordict = dd[-1]

    ## Transform coordinate to galactic frame
    ec = ephem.Equatorial(infodict['alpha']*pi/180.0,
                          infodict['delta']*pi/180.0)
    gc = ephem.Galactic(ec)

    latit = gc.lat*180.0/pi
    longit = gc.lon*180.0/pi
    ## 

    ## Compute size of simulation from minsizefile
    if 'simsize' in hyppriordict:
        solid_angle_arcsec = hyppriordict['simsize'] * 3600**2.0
        solid_angle_degree_squared = hyppriordict['simsize']

    else:
        x, y = np.loadtxt(minsizefile, unpack = True)
        solid_angle_arcsec =  interpolate.interp1d(x, y)(latit)
        solid_angle_degree_squared = solid_angle_arcsec / 3600**2

    ## Compute maximum magnitude to simulate
    target_magnitude = hyppriordict.pop('target_magnitude')
    max_magnitude = target_magnitude - 2.5 * log10(hyppriordict['transit_depth'])
    ## Get name of contrast curve
    sensi_name = datadict['HYPpriors']['sensifile']

    ## Get name of simulation
    file_name = datadict['HYPpriors']['simulfile']

    ## If simulation file is missing, run Besancon model.
    if not os.path.exists(file_name):
        dm.besancon_model_form(name_target = infodict['name'],
                               model = hyppriordict['model'],
                               minimum_distance = hyppriordict['MinDist'],
                               maximum_distance = infodict['MaxDist'],
                               longit=longit, latit=latit, 
                               solid_angle=solid_angle_degree_squared,
                               maximum_magnitude = max_magnitude,
                               minimum_magnitude = target_magnitude - 3,
                               email = infodict['email'],
                               outfile = file_name)

    return mp.compute_prior(file_name, sensi_name,
                            model = hyppriordict.pop('model'),
                            mag_target = target_magnitude,
                            mag_max = max_magnitude,
                            solid_angle_arc = solid_angle_arcsec,
                            **hyppriordict)


def initialize_priors(pastisfile):

    f = open(pastisfile, 'r')
    dd = pickle.load(f)
    f.close()
    
    infodict = dd[0]
    datadict = dd[2]
    hyppriordict = dd[-1]


    ##################
    global lowmass_period_dist, giant_period_dist_FGK
    global giant_period_dist_OBA, giant_period_dist_M
    global planet_period_dist
    global binary_period_dist_FGK, multiplicity_fraction_FGK
    global binary_period_dist_OBA, multiplicity_fraction_OBA
    global binary_period_dist_M, multiplicity_fraction_M


    ######################
    ## Planet period distributions
    ######################
    
    ### Mayor statistics ####
    if infodict['PS'] == 'Mayor':

        # Interpolation of the frequency of giant planets for
        # FGK stars from Mayor et al 2011
        x, y = np.loadtxt(PRIORdict[infodict['PS']][1], unpack = True)
        y = y/100
        giant_period_dist_FGK = interpolate.interp1d(x, y)


        ## Read period distribution for small planets and compute CDF
        x, y = np.loadtxt(PRIORdict[infodict['PS']][0], unpack = True)
        total = sum(y)
        for i in range(len(y)):
            if i == 0 :
                y[i] = y[i]/total
            if i > 0 :
                y[i] = (y[i]/ total) + y[i-1]
    
        # Interpolate obtained CDF
        lowmass_period_dist = interpolate.interp1d(x, y)

        ## Compatibility with Fressin
        planet_period_dist = {}
        for ptype in ['Earth', 'SuperEarth', 'SmallNeptune', 'LargeNeptune',
                      'Jupiter']:
            if ptype != 'Jupiter':
                planet_period_dict[ptype] = lowmass_period_dist
            else:
                planet_period_dict[ptype] = giant_period_dist_FGK
        


    ### Fressin statistics ####
    elif infodict['PS'] == 'Fressin':

        # Read cumulative data
        G, LN, SN, SE, E = n.loadtxt(PRIORdict[infodict['PS']][0],
                                     skiprows = 2, unpack = True,
                                     usecols = (0, 2, 4, 6, 8))

        # Define period array (SEE THIS IN DETAIL!!)
        peredges = n.loadtxt(PRIORdict[infodict['PS']][1], skiprows = 2)
        P = 0.5*(peredges[1:] + peredges[:-1])
        
        ## Interpolate fraction of stars that have planets for different
        ## size ranges
        planet_period_dist = {}

        for ptype, pdist in zip(('Earth', 'SuperEarth', 'SmallNeptune',
                                 'LargeNeptune', 'Jupiter'),
                                (E, SE, SN, LN, G)):
            ##
            cond = n.array(map(n.isfinite, pdist))
            planet_period_dist[ptype] = interpolate.interp1d(P[cond],
                                                             pdist[cond]*1e-2)
            
        giant_period_dist_FGK = planet_period_dist['Jupiter']
        """
        ## Interpolate fraction of stars that have small planets
        ## A dictionary is used to manage different sizes of planets.
        lowmass_period_dist = interpolate.interp1d(P[:-3], E[:-3])
        """
        
    # Planet fraction of OBA stars from  
    giant_period_dist_OBA = linear(0.14/1000, 0) 

    # Planet fraction of M stars from Bonfils
    giant_period_dist_M = linear(0.05/1000, 0) 

    
    ######################
    ## Binary period distributions
    ######################
    # Binary fraction of FGK stars from  
    multiplicity_fraction_FGK = 1 - 0.56 # Including all multiple stars
    # Period distribution (in logP!)
    binary_period_dist_FGK = st.norm(5.03, 2.28)

    # Binary fraction of M stars from  
    multiplicity_fraction_M = 0.35
    # Period distribution (in logP!)
    # Fischer & Marcy 92 say: "period distribution of G and M star is not
    # greatly different within the errors."
    binary_period_dist_M = st.norm(5.03, 2.28) 

    # Binary fraction of OBA stars from  
    multiplicity_fraction_OBA = 0.725
    # Period distribution (in logP!)
    binary_period_dist_OBA = st.norm(5.03, 2.28)

    #####
    # Sensitivity curve and extinction curve
    #####
   
    # Interpolate contrast curve
    global contrast_curve
    rdist, deltamag = np.loadtxt(datadict['HYPpriors']['sensifile'],
                                 unpack = True, skiprows = 2)
    contrast_curve = interpolate.interp1d(deltamag, rdist,
                                          fill_value = n.max(rdist),
                                          bounds_error = False)

    # Interpolate extinction curve
    global ext_curve

    # Create file name based on input parameters
    extpath = os.path.join(libpath, 'GIE')
    EXTfile = os.path.join(extpath, 'ext_ra%.3f_dec%.3f_max%.d_step%.1f.txt'%\
                           (infodict['alpha'], infodict['delta'],
                            infodict['MaxDist'], infodict['EXTstep'])
                           )
    dist, redening, Av = np.loadtxt(EXTfile, unpack = True, skiprows = 2)
    ext_curve = interpolate.interp1d(dist, redening, bounds_error = False,
                                     fill_value = n.max(redening))


def linear(a, b):
    '''
    Function that allows the creation of a linear function
    '''
    def result(x):
        return a * x + b
    return result
