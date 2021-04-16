import numpy as n
import scipy as sp
import math
from scipy import interpolate

from .. import *
from . import *

SPTdict = {'O' : 1, 'B' : 2, 'A' : 3, 'F' : 4, 'G' : 5, 'K' : 6, 'M' : 7}

def compute_prior(model_file, sensi_name, model, mag_target, mag_max,
                  solid_angle_arc, niterations = 100,
                  **kwargs):

    mag_target = float(mag_target)
    photband = kwargs.pop('photband', 'J')
    
    ## Read keyword parameters
    period = kwargs.pop('period')
    error_period = kwargs.pop('error_period')
    #
    M1 = kwargs.pop('Mstar') # In Msun
    R1 = kwargs.pop('Rstar')  # In Rsun
    spectral_type = kwargs.pop('spectral_type')    
    spectral_number = SPTdict[spectral_type]

    # Radius of transiting planet candidate.
    radius_planet = kwargs.pop('Rplanet', 0.0) # Rearth

    # Radius of false positive blended planet.
    radius_blended_planet = kwargs.pop('Rbplanet', 2.0) # Rjup

    # Minimum magnitud difference with target star for background blends
    deltamag = kwargs.pop('deltamag', 1.0)

    # Size of magniude bins over which to compute probabilities
    binsize = kwargs.pop('binsize', 0.5)
    
    # Get column-matching dictionary for galactic models and make it global
    global cols
    cols = parameter_column(model)
    
    ## Read results from Galactic model simulation
    results = n.genfromtxt(model_file, skip_header = 83, skip_footer = 4)

    print('Running %d iterations.'%niterations)
    
    #######
    ### BEB
    #######
    pBEB, pBB, nBin, nEB, nStars, nRegions = proba_BEB(results,
                                                       mag_target + deltamag,
                                                       mag_max, mag_target,
                                                       solid_angle_arc,
                                                       period, error_period,
                                                       binsize = binsize,
                                                       model = model,
                                                       photband = photband,
                                                       niterations = niterations)


    #######
    ### BTP
    #######
    pBTP, pBP, nPla, nTP, nStars, nRegions = proba_BTP(results,
                                                       mag_target + deltamag,
                                                       mag_max, mag_target,
                                                       solid_angle_arc,
                                                       period, error_period,
                                                       binsize = binsize,
                                                       photband = photband,
                                                       model = model)


    ##########
    ### TRIPLE
    ###########
    pTRIPLE = proba_TRI(spectral_number, M1, period, error_period)


    #######
    ### PIB
    #######
    pPIB, pPl, pTP = proba_PIB(spectral_number, M1, period, error_period,
                               Rp = radius_blended_planet,
                               niterations = niterations)

    ##########
    ### PLANET
    ###########
    pPLA = proba_PLA(M1, R1, radius_planet, period, error_period)

    
    ## Print results
    print('='*80)
    print('BB: background binary')
    print('BEB: background eclipsing binary')
    
    print('p(BB): %e' % pBB)
    print('p(BEB): %e +/- %e' % (n.mean(pBEB), n.std(pBEB)))

    print('='*80)
    print('BP: background planetary system')
    print('BTP: background transiting planetary system')

    print('p(BP): %e' % pBP)
    print('p(BTP): %e +/- %e' % (n.mean(pBTP), n.std(pBTP)))

    print('='*80)
    print('PiB : transiting planetary system in binary system')
    print('p(PiB): %e +/- %e' %(n.mean(pPIB), n.std(pPIB)))
    print('p(TRIPLE): %e' %pTRIPLE)
    print('p(PLANET): %e' %pPLA)

    print('='*80)
    print('Total number of stars used for simulation: %d'%nStars)
    print('Number of binary stars in simulation: %e' %nBin)
    print('Number of eclipsing binary stars in simulation: %e +/- %e' %(n.mean(nEB), n.std(nEB)))
    print('Number of planetary systems in simulation: %e' %nPla)
    print('Number of transiting planetary systems in simulation: %e' % n.mean(nTP))
    return nStars, nBin, nEB, nRegions


def parameter_column(model):
    '''
    Function that associates the column number to the paremeter for the two models of Besancon and Trilegal
    '''
    if model == 'TRI':
        param = {'logAge': 1, 'MH': 2, 'Mass': 3, 'logL': 4, 'logTe': 5,
                 'logg': 6, 'Av': 8, 'mbol': 10, 'mag': 17}
        return param
        
    elif model == 'BES':
        param = {'dist': 0, 'Mv': 1, 'lumclass': 2, 'spectral_number': 3, 
                 'logTe': 4, 'logg': 5, 'Age': 6, 'Mass': 7, 'MH': 13,   
                  'Av': 16, 'Mbol': 17, 'mag': 12}
        return param

    else:
        raise NameError('Model parameter not recognised')
        



def proba_BEB(data, mag_min, mag_max, mag_target, solid_angle_arc, period, error_period, binsize = 0.5, model = 'BES', photband = 'J', niterations = 100):
    '''
    Function that separates the data into bins of magnitudes and then
    calculates the probability of having a blended eclipsing binary
    for the given magnitude.
    Returns the sum of the probability of blended eclipsing binary,
    blended binary, the number of binary stars and the number of eclipsing
    stars, the probability of eclipse and the total number of stars used
    for the calculations.
    '''
    magbins = n.arange(mag_min, mag_max + binsize, binsize)

    ## If the column-matching dicitonary is not defined, define it
    if 'cols' not in locals():
        cols = parameter_column(model)
    
    ### For all modeled stars, compute mass and radius of a potential companion
    ### and the corresponding eclipse probability
    
    # Compute extinguished magnitude for all modeled stars
    mag_model = data[:, cols['mag']]
    distance = data[:, cols['dist']]*1e3 ## in pc

    ## Mo = M + Rv * Alambda / Av * E(B-V)(dist) = M + Alambda / Av * Av(dist)
    mag_ext = mag_model + 3.1*A_Av[photband]*ext_curve(distance)

    # Compute luminosity (using model, unextinguished bolometric magnitude)
    Mbol = data[:, cols['Mbol']]
    L = 10**(-0.4*(Mbol - Mbolsun)) # In solar units

    # Compute radius using Teff and luminosity of all stars
    Teff = 10**(data[:, cols['logTe']])
    R1 = n.sqrt(L*(5777.0/Teff)**4.0) # In solar units
    M1 = data[:, cols['Mass']] # In solar units

    # Obtain binary probability in complete model based on spectral type of
    # stars. Use proba_binary_spt function.
    sp_type = data[:, cols['spectral_number']]
    bin_prob = proba_binary_spt(sp_type, period, error_period)
    
    #### Prepare arrays containing final result
    nBin = n.zeros(len(magbins) - 1)
    pBB = n.zeros(len(magbins) - 1)

    nEB = n.zeros((niterations, len(magbins) - 1))
    pBEB = n.zeros((niterations, len(magbins) - 1))

    nstars = n.zeros(len(magbins) - 1)
    nRegions =  n.zeros(len(magbins) - 1)

    # Iterate over magnitude bin
    for j, mag in enumerate(magbins[:-1]):
        # Write the condition for a star to be in this bin
        cond =  (mag_ext > mag) * (mag_ext <= mag + binsize)

        # Make list with stars in this magnitude bin
        data_mag = data[cond, :]

        nstars[j] = len(data_mag)

        # Compute deltaMagnitude with respecto to target, and corresponding
        # number of regions.
        delta_mag = (mag + binsize/2.0) - mag_target
        area = (2*contrast_curve(delta_mag))**2 # Area of square around target
        nRegions[j] = float(solid_angle_arc)/area
    
        # For this bin, the number of binaries is the sum of the probabilities
        # of each star.
        nBin[j] = n.sum(bin_prob[cond])

        ## Do iteration to consider different possible binaries
        for i in range(niterations):
    
            # Draw a random mass ratio and get radius of secondary
            # from mass-radius relationship (for main-sequence stars).
            M2 = M1 * q_def()
            R2 = massradiusrel(M2)

            # Compute semi-major axis and eclipse probability
            a = semi_major_axis(M1, M2, period) # In Rsun
            eclip_prob = (R1 + R2)/a

            # Skip bins if no model stars fall in it.
            #if len(data_mag) == 0: continue

            # Multiply binary probability with eclip_prob. The sum over all
            # modeled star in this bin is the number
            # of EBs in this bin and iteration.
            nEB[i, j] = n.sum(bin_prob[cond] * eclip_prob[cond])

            # Compute probability for final number of eclipsing binaries
            # for this bin and iteration
            # see Blend Probabilities by Rodrigo F. Diaz, 2011
            pBEB[i, j] = 1.0 - ((nRegions[j] - 1.0)/nRegions[j])**nEB[i, j]

        # Compute probability of blended binary
        pBB[j] = 1.0 - ((nRegions[j] - 1.0)/nRegions[j])**nBin[j]


    ## Return sum over bins
    return map(n.sum, (pBEB, pBB, nBin, nEB, nstars, nRegions),
               (1, 0, 0, 1, 0, 0))
    """
    return mag_ext, M1, R1, bin_prob, nBin, pBB, nEB, pBEB, nRegions, magbins, nstars
    """
    
def proba_BTP(data, mag_min, mag_max, mag_target, solid_angle_arc, period, error_period, binsize = 0.5, photband = 'J', model = 'BES'):
    '''
    Function that separates the data into bins of magnitudes and then
    calculates the probability of having a blended transiting planet for
    the given magnitude, it returns the sum of the probability of
    blended transiting planet, blended planet and the probability of eclipse.
    '''    
    magbins = n.arange(mag_min, mag_max + binsize, binsize)

    ## If the column-matching dicitonary is not defined, define it
    if 'cols' not in locals():
        cols = parameter_column(model)

    ### For all modeled stars, compute mass and radius of a potential companion
    ### and the corresponding eclipse probability
    
    # Compute extinguished magnitude for all modeled stars
    mag_model = data[:, cols['mag']]
    distance = data[:, cols['dist']]

    ## Mo = M + Rv * Alambda / Av * E(B-V)(dist) = M + Alambda / Av * Av(dist)
    mag_ext = mag_model + 3.1*A_Av[photband]*ext_curve(distance)

    # Compute luminosity (using model, unextinguished bolometric magnitude)
    Mbol = data[:, cols['Mbol']]
    L = 10**(-0.4*(Mbol - Mbolsun)) # In solar units

    # Compute radius using Teff and luminosity
    Teff = 10**(data[:, cols['logTe']])
    R1 = n.sqrt(L*(5777.0/Teff)**4.0) # In solar units
    M1 = data[:, cols['Mass']] # In solar units

    # Compute transit probability by a zero-mass, 2-Rjup planet.
    R2 = 2*Rjup2Rsun

    # Compute semi-major axis and eclipse probability
    a = semi_major_axis(M1, 0.0, period) # In Rsun; mass of the planet neglected
    transit_prob = (R1 + R2)/a

    # Obtain the probability of having a planet for the complete model
    # based on spectral type of
    # stars. Use proba_planet_spt function.
    sp_type = data[:, cols['spectral_number']]
    planet_prob = proba_planet_spt(sp_type, period, error_period)
    
    #### Prepare arrays containing final result
    nPlanet = n.zeros(len(magbins) - 1)
    pBP = n.zeros(len(magbins) - 1)
    
    nTP = n.zeros(len(magbins) - 1)
    pBTP = n.zeros(len(magbins) - 1)

    nstars = n.zeros(len(magbins) - 1)
    nRegions =  n.zeros(len(magbins) - 1)

    # Iterate over magnitude bin
    for j, mag in enumerate(magbins[:-1]):
        # Write the condition for a star to be in this bin
        cond =  (mag_ext > mag) * (mag_ext <= mag + binsize)

        # Make list with stars in this magnitude bin
        data_mag = data[cond, :]

        nstars[j] = len(data_mag)

        # Compute deltaMagnitude with respecto to target, and corresponding
        # number of regions.
        delta_mag = (mag + binsize/2.0) - mag_target
        area = (2*contrast_curve(delta_mag))**2 # Area of square around target
        nRegions[j] = float(solid_angle_arc)/area
        
        # For this bin, the number of planets ins the sum of the probabilities
        # of individual stars
        nPlanet[j] = n.sum(planet_prob[cond])

        # Compute the probability of blended planet
        # (both transiting and not transiting).
        pBP[j] = 1.0 - ((nRegions[j] - 1.0)/nRegions[j])**nPlanet[j]

        ## Compute the number of transiting planet as the sum of
        ## the product of planet_prob and transit_prob
        nTP[j] = n.sum(planet_prob[cond] * transit_prob[cond])
        
        # Compute probability for final number of transiting planets
        # for each bin and iteration
        # see Blend Probabilities by Rodrigo F. Diaz, 2011
        pBTP[j]  = 1.0 - ((nRegions[j] - 1.0)/nRegions[j])**nTP[j]

    ## Return sum over bins
    return map(n.sum, (pBTP, pBP, nPlanet, nTP, nstars, nRegions))


def proba_TRI(sp_type, M1, period, error_period):
    """
    Function that calculates the probability of having a triple system,
    given the spectral type, and mass of the primary.
    Then compute the probability that the binary eclipses
    """
    # We first calculate the probability of the target being a binary.
    prob_bin  = proba_binary_spt(sp_type, period, error_period)

    # 20 % of binaries are in fact triple systems from Rappaport+ 2013 ApJ.768.33
    prob_triple = prob_bin * 0.2 

    # Compute mass and radius of other objects in the system
    M2 = M1 * q_def()
    M3 = M2 * q_def()
    R2 = massradiusrel(M2)
    R3 = massradiusrel(M3)
    
    # Compute semi-major axis and eclipse probability
    a = semi_major_axis(M2, M3, period) # In Rsun
    p_eclip = (R2 + R3)/a
    return prob_triple * p_eclip


def proba_PIB(sp_type, M1, period, error_period, Rp = 2, niterations = 100):

    """
    Function that calculates the probability of having a planetary system in a binary system, given spectral type and mass of target star.
    """

    # Calculate the that the star has a binary, with a period at least 
    # 10 times larger than the candidate period.
    bin_prob = proba_binary_spt(sp_type, period*10.0, error_period = 0.0)
    
    # Prepare arrays
    pPla = n.zeros(niterations)
    pTP = n.zeros(niterations)
    M2 = n.zeros(niterations)
    
    for i in range(niterations):
        # The compute the probability that the secondary star has
        # a planetary system

        # Draw a random mass ratio and get radius of secondary
        # from mass-radius relationship (for main-sequence stars).
        M2[i] = M1 * q_def()
        R2 = massradiusrel(M2[i])

        # Compute the probability that the secondary hosts a giant planet
        # at the period of the candidate.
        planet_prob = proba_planet_mass(M2[i], period, error_period)

        pPla[i] = planet_prob

        a = semi_major_axis(M2[i], 0, period)  ## Neglects planet mass
        transit_prob = (R2 + Rp*Rjup2Rsun)/a

        pTP[i] = planet_prob * transit_prob

    ## Return the final probability, the product of all three probabilities;
    ## this assumes these probabilities are independent, which might not be
    ## the case between prob_bin and planet_prob
    return bin_prob * pTP, pPla, pTP

    
def proba_PLA(mass1, radius1, radius2, period, error_period):
    
    """
    Function that calculates the probability of having a small-mass
    planetary system.
    """

    if radius2 >= 6 and radius2 < 22:
        lowmass_period_dist = planet_period_dist['Jupiter']
        print('Jupiter')
    elif radius2 >= 4 and radius2 < 6:
        lowmass_period_dist = planet_period_dist['LargeNeptune']
        print('LargeNeptune')
    elif radius2 >= 2 and radius2 < 4:
        lowmass_period_dist = planet_period_dist['SmallNeptune']
        print('SmallNeptune')
    elif radius2 >= 1.25 and radius2 < 2:
        lowmass_period_dist = planet_period_dist['SuperEarth']
        print('SuperEarth')
    elif radius2 < 1.25:
        lowmass_period_dist = planet_period_dist['Earth']
        print('Earth')
    else:
        raise ValueError('Radius of transiting planet too large; has to be smaller than 22 Rearth')

    p_inf = lowmass_period_dist(period - error_period)
    p_sup = lowmass_period_dist(period + error_period)
    p_fin = (p_sup - p_inf)
    
    ## Compute semi-major axis (planet mass neglected next to stellar mass).
    a = semi_major_axis(mass1, 0.0, period) # In Rsun

    # We then calculate the probability that the planet transits
    return p_fin* (radius1 + radius2*Rearth2Rjup*Rjup2Rsun)/a


def proba_binary_spt(sp_type, period, error_period):
    """
    For a list of modeled stars, compute the probability that they are
    binaries, based in their spectral type.
    """
    if not n.iterable(sp_type):
        sp_type = n.array([sp_type,])
    else:
        sp_type = n.array(sp_type)
        
    bin_prob = n.zeros(len(sp_type))
    for j, spt in enumerate(sp_type):

        if spt < 4.0:
            # If the star is a OBA we calculate the associated probability
            p_inf = binary_period_dist_OBA.cdf(log10(period - error_period))
            p_sup = binary_period_dist_OBA.cdf(log10(period + error_period))
            normalization = multiplicity_fraction_OBA

        elif spt >= 4.0 and spt < 7.0:
            # If the star is a FGK we calculate the associated probabiliy
            p_inf = binary_period_dist_FGK.cdf(log10(period - error_period))
            p_sup = binary_period_dist_FGK.cdf(log10(period + error_period))
            normalization = multiplicity_fraction_FGK
                

        elif spt >= 7.0 and spt < 8.0:
            # If the star is a M we calculate the associated probability
            p_inf = binary_period_dist_M.cdf(log10(period - error_period))
            p_sup = binary_period_dist_M.cdf(log10(period + error_period))
            normalization = multiplicity_fraction_M

        else:
            ## For weird stars, set everything to zero
            p_inf = 0.0
            p_sup = 0.0
            normalization = 0.0

        ## If the period error is zero, assume the given period is the lower
        ## limit of the required period range
        if error_period == 0.0:
            p_fin  = (1 - p_inf) * normalization
        else:
            p_fin = (p_sup - p_inf) * normalization

        bin_prob[j] = p_fin

    return bin_prob

def proba_planet_spt(sp_type, period, error_period = 0.0):

    # Check is sp_type is iterable
    if not n.iterable(sp_type):
        sp_type = n.array([sp_type,])
    else:
        sp_type = n.array(sp_type)

    planet_prob = n.zeros(len(sp_type))
    ## Iterate for all stars
    for j, spt in enumerate(sp_type):
        
        if spt < 4.0:
            # If the star is a OBA we calculate the associated probability
            p_inf = giant_period_dist_OBA(period - error_period)
            p_sup = giant_period_dist_OBA(period + error_period)

        elif spt >= 4.0 and spt < 7.0:
            # If the star is a FGK we calculate the associated probability
            p_inf = giant_period_dist_FGK(period - error_period)
            p_sup = giant_period_dist_FGK(period + error_period)

        elif spt >= 7.0 and spt < 8.0:
            # If the star is a M we calculate the associated probability
            p_inf = giant_period_dist_M(period - error_period)
            p_sup = giant_period_dist_M(period + error_period)

        else:
            ## For weird stars, set everything to zero
            p_inf = 0.0
            p_sup = 0.0
            
        ## If the period error is zero, assume the given period is the lower
        ## limit of the required period range
        if error_period == 0.0:
            ### THIS IS WRONG!!! INSTEAD OF ONE, WE SHOULD HAVE THE
            ### PLANET FRACTION FOR EACH TYPE OF STAR; SEE BINARY EQUIVALENT
            p_fin  = 1 - p_inf
        else:
            p_fin = p_sup - p_inf

        # Add the probability to the number of planets in this bin
        planet_prob[j] = p_fin

    return planet_prob

def proba_planet_mass(mass, period, error_period = 0.0):
    """
    Planet probability as a function of host star mass.
    """
    # Check is mass is iterable
    if not n.iterable(mass):
        mass = n.array([mass,])
    else:
        mass = n.array(mass)

    planet_prob = n.zeros(len(mass))

    ## Iterate for all stars in this bin
    for j, mm in enumerate(mass):

        ## Depending on the mass, compute planet probability
        if mm > 0.45 and mm < 1.4:
            p_inf = giant_period_dist_FGK(period - error_period)
            p_sup = giant_period_dist_FGK(period + error_period)

        elif mm > 1.4:
            p_inf = giant_period_dist_OBA(period - error_period)
            p_sup = giant_period_dist_OBA(period + error_period)

        elif mm < 0.45:
            p_inf = giant_period_dist_M(period - error_period)
            p_sup = giant_period_dist_M(period + error_period)

        else :
            ### CHECK IF THIS IS REASONABLE!!!!
            planet_prob[j] = 0

        ## If the period error is zero, assume the given period is the lower
        ## limit of the required period range
        if error_period == 0.0:
            ### THIS IS WRONG!!! INSTEAD OF ONE, WE SHOULD HAVE THE
            ### PLANET FRACTION FOR EACH TYPE OF STAR; SEE BINARY EQUIVALENT
            planet_prob[j]  = 1 - p_inf
        else:
            planet_prob[j] = p_sup - p_inf

    return planet_prob

        
    
def semi_major_axis(mass1, mass2, period):
    a = ((G*((mass1+mass2)*Msun)*(period*24*3600)**2)/(4.0* n.pi**2))**(1.0/3.0)
    return a/Rsun

def q_def():
    """
    Draw random mass ratio from realistic distribution.
    
    This function that gives the ratio mass of a binary system based on the 
    statistics from Rghavan et al 2010 (Fig. 16, left panel).
    
    The distribution is modelled using a piece-wise constant distribution.
    
    A number is first drawn to choose which part of the distribution the mass 
    ratio will be drawn from and then a 
    """
    rand = sp.rand()
    if rand < 0.0454:
        # Lower end of distribution
        q =  sp.rand()*0.15
    elif rand > 0.882:
        # Upper end of distribution
        q = ((sp.rand() * 0.05)+0.95)
    else:
        # Middle part of distribution
        q =  ((sp.rand() * 0.8)+0.15)
    return q

def massradiusrel(M):
    """
    Compute radius of dwarf star.
    
    Use mass-radius relationship in Allen's Astrophysical Quantities.
    """
    if (n.log10(M) > 1.3).any() or (n.log10(M) < -1.0).any():
        #print('Warning! Mass-radius relation not calibrated for this mass range')
        pass
        
    ## Compute standard M-R relationship
    radius =  n.where(n.log10(M) > 0.12, 10**(0.640*n.log10(M) + 0.011),
                      10**(0.917*n.log10(M) - 0.020))

    return radius
    
"""
def rayondef(mass2):
    # Calculates the radius of a dwarf star given it's mass from Allens astrophysical quantities
    if (math.log1p(mass2) < 1.3) and (math.log1p(mass2)) > 0.12:
        rayon2 = Rsun*10**((0.640*math.log1p(mass2))+0.011)
    elif (math.log1p(mass2) < 0.12) and (math.log1p(mass2) < 0.12) > -1.0:
        rayon2 = Rsun*10**((0.917*math.log1p(mass2))-0.020)
    else:
        rayon2 = 0
    return rayon2

def binary_probability(type_spectral):
    #Funtion that gives the probability of a star being a binary depending on it's spectral type, these values represent minimum values and a more precise function would take in account the uper values also
    
    proba = {'O': 0.75, 'B': 0.70, 'A': 0.70, 'F': 0.70, 'G': 0.45, 'M': 0.35, 'L': 0.2, 'T': 0.2, 'logTe': 4, 'logg': 5, 'Av': 16, 'mbol': 17, 'j': 12}
    return proba[type_spectral]

"""
