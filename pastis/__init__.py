import os
import numpy as n

from .paths import libpath, filterpath, zeromagfile
from .extlib import EMdict
from . import photometry
from . import limbdarkening as ld
from . import extinction as ext
from . import velocimetry as vel
from . import isochrones as iso

from .extlib import SAMdict, LDdict

# Define current version and print
__version__ = 'PASTIS'
version = __version__

hostname = os.getenv('HOSTNAME')
if hostname is None:
    import socket
    hostname = socket.gethostname()
    
print('Running {} on {}'.format(version, hostname))

# List of Classes and parameters
ClassesList = {'Star': {'achromatic': ['minit', 'z', 'logage', 'dens', 'teff',
                                       'logg', 'mact', 'R', 'dist', 'vsini',
                                       'v0', 'ebmv'],
                        'chromatic': ['albedo', 'ua', 'ub', 'gd', 'B',
                                      'alphaS'],
                        'iterable': {}},

               'PlanetHost': {'achromatic': ['teff', 'dens', 'z',  'dist',
                                             'vsini', 'v0', 'ebmv'],
                              'chromatic': ['albedo', 'ua', 'ub', 'gd', 'B',
                                            'alphaS'],
                              'iterable': {}},
               
               'Blend': {'achromatic': ['minit', 'z', 'logage', 'dist', 'vsini',
                                        'v0', 'ebmv'],
                         'chromatic': ['albedo', 'ua', 'ub', 'gd', 'B',
                                       'alphaS'],
                         'iterable': {}},

               'Target': {'achromatic': ['teff', 'logg', 'z', 'dist', 'vsini',
                                         'v0', 'ebmv'],
                          'chromatic': ['albedo', 'ua', 'ub', 'gd', 'B',
                                        'alphaS'],
                          'iterable': {}
                          },
                           
               'WhiteDwarf': {'achromatic': ['teff', 'logg', 'z', 'dist',
                                             'vsini', 'v0', 'ebmv'],
                              'chromatic': ['albedo', 'ua', 'ub', 'gd', 'B',
                                            'alphaS'],
                              'iterable': {}},

               'Planet': {'achromatic': ['Mp', 'Rp', 'orbital_parameters'],
                          'chromatic': ['albedo', 'f'],
                          'iterable': {}},

               'FitBinary': {'achromatic': ['kr', 'sumr', 'K1', 'q', 'v0',
                                            'vsin1', 'vsin2', 'spinorbit',
                                            'orbital_parameters', 'ar'],
                             'chromatic': ['sbr', 'ua1', 'ub1', 'ua2', 'ub2',
                                           'gd1', 'gd2', 'B1', 'B2', 'albedo1',
                                           'albedo2'],
                             'iterable': {}},

               'FitPlanet': {'achromatic': ['kr', 'sumr', 'K1', 'q', 'v0',
                                            'vsin1', 'vsin2', 'spinorbit',
                                            'orbital_parameters', 'ar'],
                             'chromatic': ['ua1', 'ub1', 'gd1', 'B1',
                                           'albedo2'],
                             'iterable': {}},

               'IsoBinary': {'achromatic': ['orbital_parameters', ],
                             'chromatic': ['sbr', 'ua1', 'ub1', 'ua2', 'ub2',
                                           'gd1', 'gd2', 'B1', 'B2', 'albedo1',
                                           'albedo2'],
                             'iterable': {'star': 2}},

               'qBinary': {'achromatic': ['orbital_parameters', 'q'],
                           'chromatic': ['sbr', 'ua1', 'ub1', 'ua2', 'ub2',
                                         'gd1', 'gd2', 'B1', 'B2', 'albedo1',
                                         'albedo2'],
                           'iterable': {'star': 1}},

               'PlanSys': {'achromatic': [],
                           'chromatic': [],
                           'iterable': {'star': 1,
                                        'planet': n.inf}},

               'Triple': {'achromatic': ['orbital_parameters', ],
                          'chromatic': [],
                          'iterable': {'object': 2}},

               'orbital_parameters': ['P', 'ecc', 'omega', 'ecos', 'esin',
                                      'incl', 'b', 'bprime', 'T0', 'Tp',
                                      'spinorbit']
               }
               
    
### List of Classes per Type
TypeList = {'binary': ['FitBinary', 'IsoBinary', 'qBinary', 'FitPlanet'],
            'plansys': ['PlanSys', ],
            'star': ['Star', 'Target', 'Blend', 'PlanetHost', 'WhiteDwarf']}

# Prior type list
PriorTypeList = ['Normal', 'TruncatedUNormal', 'TruncatedJNormal', 'Binormal',
                 'Uniform', 'LnNormal', 'LogNormal', 'Jeffreys', 'PowerLaw',
                 'DoublePowerLaw', 'AssymetricNormal', 'Sine', 'Beta',
                 'LogBinormal', 'File']

beaming = False
spotmodel = None
Nmax = None
TrefRV = 0.0
checkblendmag = False

###
# Function to initialize PASTIS
###
def initialize(*args):
    """Initialise PASTIS."""
    global inputdicts
    global beaming, Nmax, TrefRV, spotmodel, checkblendmag

    ## If no input is given, initialize using default options.
    ## No initilization of extinction
    if len(args) == 0:
        pbands = ['Kepler', 'CoRoT-W', 'Johnson-B', 'Johnson-V']
        pbands_limbdarkening = ['Kepler', 'CoRoT-W']
        photometry.initialize_phot(pbands, zeromagfile, filterpath,
                                   AMmodel=SAMdict['BT-settl'])
        #photometry.initialize_phot_WD()

        iso.interpol_tracks(EMdict['Dartmouth'])
        iso.prepare_tracks_target(EMdict['Dartmouth'])
        #isochrones.interpol_WD(os.path.join(libpath, 'AM', 'WD', 'Table_DA'))

        ld.initialize_limbdarkening(pbands_limbdarkening,
                                    ATMmodel='A',
                                    LDCfile=LDdict['Claret2011_wTESS'])

        inputdicts = [{}, {}, {}]

        beaming = False
        Nmax = None
        TrefRV = 0.0
        spotmodel = None
        checkblendmag = False
        
        return
    
    inputdicts = [args[0], args[1], args[2]]

    infodict = args[0]
    datadict = args[1]
    objectdict = args[2]

    # Compute beaming?
    try:
        beaming = infodict['beaming']
    except KeyError:
        beaming = False

    # Maximun number of points computed by EBOP
    try:
        Nmax = infodict['Nmax']
    except KeyError:
        Nmax = None

    # Tiempo de referencia para RV
    try:
        TrefRV = infodict['TrefRV']
    except KeyError:
        TrefRV = 2450000.0

    # Spot model
    try:
        spotmodel = infodict['spotmodel']
    except KeyError:
        spotmodel = None

    #HARD CODED magnitude LIMITATION TO BLENDS 
    try:
        checkblendmag = infodict['checkblendmag']
    except KeyError:
        checkblendmag = False

    # Check if a target star is in objectdict
    contains_target = False
    contains_blend = False
    contains_wd = False
    for obj in objectdict.keys():
        if 'Target' in obj or 'PlanetHost' in obj or 'Star' in obj:
            contains_target = True
        if 'Blend' in obj or 'qBinary' in obj:
            contains_blend = True
        if 'WhiteDwarf' in obj:
            contains_wd = True
         
    # Construct list of photometric bands needed for PASTIS run
    corotcolor = False
    pbands_limbdarkening = []
    for kk in datadict.keys():
        if 'filter' in datadict[kk]:
            if datadict[kk]['filter'] == 'CoRoT-R' or \
                    datadict[kk]['filter'] == 'CoRoT-G' or \
                    datadict[kk]['filter'] == 'CoRoT-B':
                corotcolor = True
            else:       
                pbands_limbdarkening.append(datadict[kk]['filter'])

    if corotcolor: 
        # To compute LD for CoRoT colors, sloan LD coefficients and bandpass
        # are used
        if 'SDSS-U' not in pbands_limbdarkening:
            pbands_limbdarkening.append('SDSS-U')
        if 'SDSS-G' not in pbands_limbdarkening:
            pbands_limbdarkening.append('SDSS-G')
        if 'SDSS-R' not in pbands_limbdarkening:
            pbands_limbdarkening.append('SDSS-R')
        if 'SDSS-I' not in pbands_limbdarkening:
            pbands_limbdarkening.append('SDSS-I')
        if 'SDSS-Z' not in pbands_limbdarkening:
            pbands_limbdarkening.append('SDSS-Z')

    # In order to copy pbands_limbdarkening to a _different_ list.
    pbands = pbands_limbdarkening[:]

    # Add B and V
    if 'Johnson-B'not in pbands:
        pbands.append('Johnson-B')
    if 'Johnson-V'not in pbands:
        pbands.append('Johnson-V')
    if corotcolor:
        if 'CoRoT-W' not in pbands:
            pbands.append('CoRoT-W')

    # Add all bands needed for SED
    if 'SED' in datadict:
        pbands.extend(datadict['SED']['data']['band'])

    ## Read BT spectra if needed
    if 'SED' in datadict or contains_target or contains_blend:
        photometry.initialize_phot(pbands, zeromagfile,
                                   filterpath, AMmodel=SAMdict[infodict['SAM']])

    ## initialize RV if needed
    contains_rv = False
    for kk in datadict.keys():
        if datadict[kk]['type'] == 'RV':
            contains_rv = True
    
    if contains_rv:
        vel.initialize_RV(datadict)

    if contains_wd:
        iso.interpol_WD(os.path.join(libpath, 'AM', 'WD', 'Table_DA'))
        photometry.initialize_phot_WD()

    try:
        if infodict['EvolModel'] == 'Padova':

            print('Using Padova Stellar Evolution Models.')
            print('Imposible to use Padova tracks: they\'re not yet available.')
            return
            ###
            ## Padova
            ###
            #isochrones.initialize_isochrones(EMdict['Padova'])

        elif infodict['EvolModel'] == 'Dartmouth':

            print('Using Dartmouth Stellar Evolution Models.')
            ###
            ## Dartmouth
            ###
            if contains_blend:
                iso.interpol_tracks(EMdict['Dartmouth'])
            if contains_target:
                iso.prepare_tracks_target(EMdict['Dartmouth'])

        elif infodict['EvolModel'] == 'Geneva':

            print('Using Geneva Stellar Evolution Models.')
            ###
            ## Geneva
            ###
            if contains_blend:
                iso.interpol_tracks(EMdict['Geneva'])
            if contains_target:
                iso.prepare_tracks_target(EMdict['Geneva'])

        elif infodict['EvolModel'] == 'StarEvol':

            print('Using StarEvol Stellar Evolution Models.')
            ###
            ## StarEvol
            ###
            if contains_blend:
                iso.interpol_tracks(EMdict['StarEvol'])
            if contains_target:
                iso.prepare_tracks_target(EMdict['StarEvol'])

        elif infodict['EvolModel'] == 'Parsec':

            print('Using Parsec Stellar Evolution Models.')
            ###
            ## StarEvol
            ###
            if contains_blend:
                iso.interpol_tracks(EMdict['Parsec'])
            if contains_target:
                iso.prepare_tracks_target(EMdict['Parsec'])

        else:
            print('Evolution model not recongnized.')
            return

    except KeyError:
        print('Evolution model not specified.')
        return
    
    ra = infodict['alpha']
    dec = infodict['delta']
    ext.initialize_extinction(ra, dec, infodict['MaxDist'],
                              infodict['EXTstep'], Rv=3.1)

    if len(pbands_limbdarkening) > 0:

        if SAMdict[infodict['SAM']] == 'CK':
            atmmodel = 'A'
        elif SAMdict[infodict['SAM']] == 'BT':
            atmmodel = 'A'
            ## 'P' only z=0 for PHOENIX model, interpolation gives error
        else:
            raise NameError('Atmospheric model not recognized.')

        ld.initialize_limbdarkening(pbands_limbdarkening,
                                    ATMmodel=atmmodel,
                                    LDCfile=LDdict[infodict['LDC']])
    
    return
