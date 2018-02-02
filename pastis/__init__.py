import os
import sys
import numpy as n
from .DataTools import readdata

# Define current version and print
version = 'PASTIS_NM'
hostname = os.getenv('HOSTNAME')
print('Running {} on {}'.format(version, hostname))

# Read environment variables
PASTISroot = os.getenv('PASTISPATH')
libpath = os.getenv('PASTISLIB')
idlexec = os.getenv('IDLEXEC')
widgetpath = os.path.join(PASTISroot, 'src', version, 'widget')

# Define other useful paths
datapath = os.path.join(PASTISroot, 'datafiles')
configpath = os.path.join(PASTISroot, 'configfiles')
resultpath = os.path.join(PASTISroot, 'resultfiles')
runpath = os.path.join(PASTISroot, 'run')

# Add directory to pythonpath to read modules from fortran
sys.path.append(os.path.join(os.getenv('PASTISLIB'), 'fortran'))

# Define physical constants
G = 6.67428e-11  # m3 kg-1 s-2
c = 299792458.0  # m/s
SB = 5.67051e-8  # W m-2 K-4

# Astronomical constants
pc = 3.08568025e16  # parsec [m]
# au = 149597870691  # astronomical unit [m]
au = 149597870700  # New definition (IAU 2012) of the astronomical unit [m]
Msun = 1.98842e30  # kg
Mjup = 1.89852e27  # kg
Mearth = 5.9736e24  # kg
Rsun = 6.95508e8  # equatorial radius [m]
Rjup = 71492000.  # equatorial radius [m]
Rearth = 6378137.0  # equatorial radius [m]
R_V = 3.1  # normal extinction of the galaxy

# Bolometric solar magnitude for prior computation
Mbolsun = 4.74

# Define useful unit conversions
# Msun2Mjup = 1.047348644e3 # From IAU
# (http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
Mjup2Mearth = Mjup/Mearth
Mearth2Mjup = Mearth/Mjup
Msun2Mjup = Msun/Mjup
Mjup2Msun = Mjup/Msun
Rjup2Rearth = Rjup/Rearth
Rearth2Rjup = Rearth/Rjup
Rsun2Rjup = Rsun/Rjup
Rjup2Rsun = Rjup/Rsun
pc2Rsun = pc/Rsun

# Define PATHS
filterpath = os.path.join(libpath, 'Filters/')
zeromagfile = os.path.join(filterpath, 'ZeroFlux.dat')
setpath = os.path.join(libpath, 'SET')
ldpath = os.path.join(libpath, 'LD')
priorspath = os.path.join(libpath, 'Priors')
minsizefile = os.path.join(priorspath, 'MinimumSize_arcsec2.dat')

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


# Define Exception classes
class EvolTrackError(Exception):
    pass


class OutofIsochroneError(Exception):
    pass


class SpectrumInterpolError(Exception):
    def __init__(self, mesg, indz, indteff, indlogg, outside=False):
        self.mesg = mesg
        self.indz = indz
        self.indteff = indteff
        self.indlogg = indlogg
        self.outside = outside

    def __str__(self):
        return self.mesg

"""
class SpectrumGridError(Exception):
    def __init__(self, mesg, indz, indteff, indlogg):
        self.mesg = mesg
        self.indz = indz
        self.indteff = indteff
        self.indlogg = indlogg
"""


class EBOPparamError(Exception):
    pass


class GlobalSpectrumError(Exception):
    pass

from . import photometry
#should be imported after the SpectrumInterpolError definition

from . import velocimetry


# Dictionaries
LDdict = {'Claret2011': os.path.join(ldpath, 'claret2011ab.txt')}
EMdict = {'Dartmouth': os.path.join(setpath, 'Dartmouth.trk'),
          'Geneva': os.path.join(setpath, 'Geneva.trk'),
          'StarEvol': os.path.join(setpath, 'Starevol.trk'),
          'Parsec': os.path.join(setpath, 'Parsec.trk')}
SAMdict = {'CastelliKurucz': 'CK', 'BT-settl': 'BT'}
PRIORdict = {'Mayor': [os.path.join(priorspath, 'PDF_SmallPlanets_Mayor.dat'),
                       os.path.join(priorspath, 'CDF_GiantPlanets_Mayor.dat')],
             'Fressin': [os.path.join(priorspath,
                                      'Fressin_cumulative_occurrence.txt'),
                         os.path.join(priorspath,
                                      'Fressin_period_bin_edges.txt')]}
# GIEdict = {'AmoresLepine2005': ''}

beaming = False
Nmax = None
TrefRV = 0.0
spotmodel = 'Macula'
checkblendmag = False


###
# Function to initialize PASTIS
###
def initialize(*args):
    """
    Function to initialise PASTIS
    """
    
    global inputdicts
    global beaming, Nmax, TrefRV, spotmodel, checkblendmag

    ## If no input is given, initialize using default options.
    ## No initilization of extinction
    if len(args) == 0:
        pbands = ['Kepler', 'CoRoT-W', 'Johnson-B', 'Johnson-V']
        pbands_limbdarkening = ['Kepler', 'CoRoT-W']
        photometry.initialize_phot(pbands, zeromagfile, filterpath,
                                   AMmodel=SAMdict['BT-settl'])
        photometry.initialize_phot_WD()

        import isochrones

        isochrones.interpol_tracks(EMdict['Dartmouth'])
        isochrones.prepare_tracks_target(EMdict['Dartmouth'])
        isochrones.interpol_WD(os.path.join(libpath, 'AM', 'WD', 'Table_DA'))
        import limbdarkening

        limbdarkening.initialize_limbdarkening(pbands_limbdarkening,
                                               ATMmodel='A',
                                               LDCfile=LDdict['Claret2011'])

        inputdicts = [{}, {}, {}]

        beaming = False
        Nmax = None
        TrefRV = 0.0
        spotmodel = 'Macula'
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
        spotmodel = 'Macula'

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
        velocimetry.initialize_RV(datadict)

    import isochrones

    if contains_wd:
        isochrones.interpol_WD(os.path.join(libpath, 'AM', 'WD', 'Table_DA'))
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
                isochrones.interpol_tracks(EMdict['Dartmouth'])
            if contains_target:
                isochrones.prepare_tracks_target(EMdict['Dartmouth'])

        elif infodict['EvolModel'] == 'Geneva':

            print('Using Geneva Stellar Evolution Models.')
            ###
            ## Geneva
            ###
            if contains_blend:
                isochrones.interpol_tracks(EMdict['Geneva'])
            if contains_target:
                isochrones.prepare_tracks_target(EMdict['Geneva'])

        elif infodict['EvolModel'] == 'StarEvol':

            print('Using StarEvol Stellar Evolution Models.')
            ###
            ## StarEvol
            ###
            if contains_blend:
                isochrones.interpol_tracks(EMdict['StarEvol'])
            if contains_target:
                isochrones.prepare_tracks_target(EMdict['StarEvol'])

        elif infodict['EvolModel'] == 'Parsec':

            print('Using Parsec Stellar Evolution Models.')
            ###
            ## StarEvol
            ###
            if contains_blend:
                isochrones.interpol_tracks(EMdict['Parsec'])
            if contains_target:
                isochrones.prepare_tracks_target(EMdict['Parsec'])

        else:
            print('Evolution model not recongnized.')
            return

    except KeyError:
        print('Evolution model not specified.')
        return
    
    import extinction
    ra = infodict['alpha']
    dec = infodict['delta']
    extinction.initialize_extinction(ra, dec, infodict['MaxDist'],
                                     infodict['EXTstep'], Rv=3.1)

    if len(pbands_limbdarkening) > 0:
        import limbdarkening

        if SAMdict[infodict['SAM']] == 'CK':
            atmmodel = 'A'
        elif SAMdict[infodict['SAM']] == 'BT':
            atmmodel = 'A'
            ## 'P' only z=0 for PHOENIX model, interpolation gives error
        else:
            raise NameError('Atmospheric model not recognized.')

        limbdarkening.initialize_limbdarkening(pbands_limbdarkening,
                                               ATMmodel=atmmodel,
                                               LDCfile=LDdict[infodict['LDC']])
    
    return
