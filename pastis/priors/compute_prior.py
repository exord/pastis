#! /usr/bin/python2.7

import moduleprior as mp  # module with all the functions used to calculate the pobability
import numpy as np
from scipy import interpolate

from .. import *
from . import *

SPTdict = {'M' : 3, 'K' : 4, 'G' : 5, 'F' : 7, 'A' : 8, 'B' : 9, 'A' : 10}

def prior(model_name, sensi_name, model, target_magnitude, magni_max, solid_angle_arc, number_iterations = 100, sensitivity_curve=True, parameter_name='Mass', parameter_limit=0, parameter=False, **kwargs):

    period = kwargs.pop('period')
    error_period = kwargs.pop('error_period')
    #
    mass_init = kwargs.pop('Mstar')
    rayon_etoile = kwargs.pop('Rstar')*Rsun
    spectral_type = kwargs.pop('spectral_type')    
    spectral_number = SPTdict[spectral_type]
    
    rayon_planete = kwargs.pop('Rplanet')*Rjup

    
    results = np.genfromtxt(model_name, skip_header = 84, skip_footer = 8)

    array_proba_blended_binary = []
    array_proba_blended_ecplising_binary = []
    array_number_of_binary_stars = []
    array_number_of_blended_binary_stars = []
    array_proba_blended_planet = []
    array_proba_blended_transiting_planet = []
    array_probability_eclipsing_binary = []
    array_probability_planet_number = []
    array_stars_used_in_calculations = []
    array_triple = []
    array_planet = []
    array_planet_in_binary =[]

    
    for i in range(number_iterations):
        percentage = number_iterations / 100.0
        print (str(i/percentage)+'%')
        target_magnitude = float(target_magnitude)
        proba_bb, proba_beb, binary, eclipsing_binary, probability_eclipse, number_stars = mp.proba_beb(results, target_magnitude, magni_max, target_magnitude, solid_angle_arc, model, period, error_period, fct_ext)

        array_proba_blended_binary.append(proba_bb)
        array_proba_blended_ecplising_binary.append(proba_beb)
        array_number_of_binary_stars.append(binary)
        array_probability_eclipsing_binary.append(eclipsing_binary)
        array_number_of_blended_binary_stars.append(probability_eclipse)
        array_stars_used_in_calculations.append(number_stars)

        proba_bp, proba_btp, planet_number = mp.proba_btp(results, target_magnitude, magni_max, target_magnitude, solid_angle_arc, model, period, error_period, fct_ext)
        array_proba_blended_planet.append(proba_bp)
        array_proba_blended_transiting_planet.append(proba_btp)
        array_probability_planet_number.append(planet_number)

        
        array_triple.append(mp.proba_triple(spectral_number, mass_init, period, error_period))
        array_planet.append(mp.proba_planete(mass_init, rayon_etoile, rayon_planete, period, error_period))
        array_planet_in_binary.append(mp.proba_planetinbinary(spectral_number, spectral_type, mass_init, rayon_planete, period, error_period))

    #proba_bin_2 = np.mean(array_bin) / np.mean(array_bin_cut)


    print('probabilite bb : %e +/- %e ' % (np.mean(array_proba_blended_binary), np.std(array_proba_blended_binary)))
    print('probabilite beb : %e +/- %e ' % (np.mean(array_proba_blended_ecplising_binary), np.std(array_proba_blended_ecplising_binary)))

    print('======================================================================' )

    print('probabilite bp: %e +/- %e' % (np.mean(array_proba_blended_planet), np.std(array_proba_blended_planet)))
    print('probabilite btp: %e +/- %e' % (np.mean(array_proba_blended_transiting_planet), np.std(array_proba_blended_transiting_planet)))

    print('======================================================================' )

    print('probabilite pib: %e +/- %e' %(np.mean(array_planet_in_binary), np.std(array_planet_in_binary)))
    print('probabilite triple: %e +/- %e' %(np.mean(array_triple), np.std(array_triple)))
    print('probabilite planet: %e +/- %e' %(np.mean(array_planet), np.std(array_planet)))

    print('Number of binary stars : %e' % np.mean(array_number_of_binary_stars))
    print('Number of planetary systems : %e' % np.mean(array_probability_planet_number))





