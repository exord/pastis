import sys
import warnings
import numpy as n
from scipy import interpolate

from .extlib import LDdict
from . import photometry as phot


def initialize_limbdarkening(pbands,
                             ATMmodel = 'A',
                             LDCfile = LDdict['Claret2011_wTESS']):

    print('Loading LimbDarkening coefficients...')

    LDCinterpol = interpol_LD(LDCfile, pbands, ATMmodel)

    global LDCs
    def LDCs(logg, teff, z, band):
        return n.array([LDCinterpol[band][0](logg, teff, z),
                        LDCinterpol[band][1](logg, teff, z)])

    print('... DONE! \n')
    return


def interpol_LD(LDCfile, photbands = ['Kepler', 'CoRoT'], ATMmodel = 'A'):
    """
    Interpolate limb darkening tables.
    
    Read limb darkenning file, interpolate over log(g), Teff, 
    and metallicity (z).

    If a combination of log(g), Teff, and z is not in the tables but within its
    limits, return the closest value in the exisiting table.

    Model: A or P for Atlas (Castelli Kurucz) or Phoenix (BT-Settl)
    Return interpolated function that returns an array with quadratic LDC for a given photmetric bandpass.
    """
    # Define dictionary that translates photbands to Claret notation
    trad_to_claret = {'CoRoT-W' : 'C', 'Kepler' : 'Kp', 'IRAC-I1' : 'S1',
                      'IRAC-I2' : 'S2', 'IRAC-I3' : 'S3', 'IRAC-I4' : 'S4',
                      'SDSS-U' : 'u*', 'SDSS-G' : 'g*', 'SDSS-R' : 'r*',
                      'SDSS-I' : 'i*', 'SDSS-Z' : 'z*', 'Johnson-U' : 'U',
                      'Johnson-B' :'B', 'Johnson-V' : 'V', 'Johnson-R' : 'R' ,
                      'Johnson-I' : 'I', '2MASS-J' : 'J', '2MASS-H' : 'H',
                      '2MASS-Ks' : 'K', 'STROMGREN-u' : 'u',
                      'STROMGREN-v' : 'v', 'STROMGREN-b' : 'b',
                      'STROMGREN-y' : 'y', 'TESS': 'T'}

    # Read data from file
    f = open(LDCfile, 'r')
    lines = f.readlines()
    f.close()

    # Number of records in LDile
    nrecord = len(lines[32: -2])

    # Initialize float arrays
    logg = n.zeros((nrecord,), float)
    Teff = n.zeros((nrecord,), float)
    z = n.zeros((nrecord,), float)
    microturb = n.zeros((nrecord,), float)
    a = n.zeros((nrecord,), float)
    b = n.zeros((nrecord,), float)

    # Initialize string arrays
    band = n.zeros((nrecord,), 'U8')
    method = n.zeros((nrecord,), 'U1')
    model = n.zeros((nrecord,), 'U1')

    ## Read variables
    for i, line in enumerate(lines[32:-2]):
        try:
            ll = line.split()
            logg[i], Teff[i], z[i], microturb[i], a[i], b[i] = list(
                    map(float, ll[:6]))
            band[i], method[i], model[i] = list(map(str, ll[6:]))

        except:
            print(ll)

    ## Prepare conditions to select a given photometric band
    condMT = n.equal(microturb, 2.0)
    assert sum(condMT) > 0, "No entries for microturbulence = 2.0"
    condMOD = model == ATMmodel
    assert sum(condMOD) > 0, "No entries for model {}".format(ATMmodel)
    condMET = (method == 'L')
    assert sum(condMET) > 0, "No entries for method = 'L' and microturbulence = 2.0"
    cond1 = condMT * condMOD * condMET
    assert sum(cond1) > 0
    
    LDCs = {}
    for ii, photband in enumerate(n.sort(photbands)):
        if ii == 0:
            sys.stdout.write('... band: %-12s'%photband)
        else:
            sys.stdout.write('\b'*12+'%-12s'%photband)
        sys.stdout.flush()

        ##
        try:
            cond2 = n.logical_and(cond1, band == trad_to_claret[photband])
        except KeyError:
            warnings.warn('Warning! Photband {} could not be initialised. '
                          'Should be ok if you plan to fit the limb darkening'
                          ' coefficients'.format(photband))
            continue
                

        loggc = n.compress(cond2, logg)
        teffc = n.compress(cond2, Teff)
        zc = n.compress(cond2, z)
        ac = n.compress(cond2, a)
        bc = n.compress(cond2, b)

        ## Create regular grid of values
        loggu = n.unique(loggc); teffu = n.unique(teffc); zu = n.unique(zc)

        length = len(loggu)*len(teffu)*len(zu)
        
        logg_reg = n.zeros((length,), float)
        teff_reg = n.zeros((length,), float)
        z_reg = n.zeros((length,), float)
        a_reg = n.zeros((length,), float)
        b_reg = n.zeros((length,), float)

        i = 0
        impossible_combinations = []
        for loggi in loggu:
            condlogg = n.less(n.abs(loggc - loggi), 1e-3)
            for teffi in teffu:
                condteff = n.less(n.abs(teffc - teffi), 0.1)
                cond_loggteff = n.logical_and(condlogg, condteff)
                for zi in zu:
                    condz = n.less(n.abs(zc - zi), 1e-3)
                    condAll = n.logical_and(cond_loggteff, condz)

                    logg_reg[i] = loggi; teff_reg[i] = teffi; z_reg[i] = zi

                    aa = n.compress(condAll, ac)
                    bb = n.compress(condAll, bc)

                    if len(aa) != 0 and len(bb) != 0:
                        a_reg[i] = aa; b_reg[i] = bb
                    else:
                        ## Complete grid

                        # Get typical step in each variable
                        dlogg = n.mean(n.diff(loggu))
                        dteff = n.mean(n.diff(teffu))
                        dz = n.mean(n.diff(zu))

                        # Compute the normalized distance to all points in file
                        dist2 = ((loggc - loggi)/dlogg)**2.0 + \
                                ((teffc - teffi)/dteff)**2.0 +\
                                ((zc - zi)/dz)**2.0

                        # Choose the closest point in the grid
                        ind = n.argsort(dist2)[0]
                        
                        a_reg[i] = ac[ind]
                        b_reg[i] = bc[ind]

                    i += 1
        
        #invalues = n.array((loggc, Teffc, zc)).transpose()
        invalues = n.array((logg_reg, teff_reg, z_reg)).transpose()

        # Generate interpolated functions
        a_interp = interpolate.LinearNDInterpolator(invalues, a_reg)
        b_interp = interpolate.LinearNDInterpolator(invalues, b_reg)

        LDCs[photband] = [a_interp, b_interp]
    return LDCs


def get_LD(teff, logg, z, photband, verbose = False):

    ## Warnings for out-of-bounds parameters
    print_warning = False
    if teff > 50000.0:
        teff=50000.0
        print_warning = True; wstring = ['Teff', 'maximum', 50000]

    if teff < 3500.00:
        teff=3500.00
        print_warning = True; wstring = ['Teff', 'minimum', 3500]

    if logg > 5.0:
        logg=5.0
        print_warning = True; wstring = ['logg', 'maximum', 5.0]

    if logg < 0.0:
        logg=0.0
        print_warning = True; wstring = ['logg', 'minimum', 0.0]

    if z > 1.0:
        z=1.0
        print_warning = True; wstring = ['z', 'maximum', 1.0]

    if z < -2.5:
        z=-2.5
        print_warning = True; wstring = ['z', 'minimum', -2.5]

    if print_warning and verbose:
        print('Sorry, Claret limb darkening tables out of range')
        print('%s changed to it %s value: %f'%(wstring[0], wstring[1], wstring[2]))


    if photband == 'CoRoT-R' or photband == 'CoRoT-G' or photband == 'CoRoT-B':
        if photband == 'CoRoT-R': ccc=0
        if photband == 'CoRoT-G': ccc=1
        if photband == 'CoRoT-B': ccc=2
        ldcu = LDCs(logg, teff, z, 'SDSS-U')
        ldcg = LDCs(logg, teff, z, 'SDSS-G')
        ldcr = LDCs(logg, teff, z, 'SDSS-R')
        ldci = LDCs(logg, teff, z, 'SDSS-I')
        ldcz = LDCs(logg, teff, z, 'SDSS-Z')
        return [phot.CoRoT_LDC_weights[ccc,0]*ldcu[0]+phot.CoRoT_LDC_weights[ccc,1]*ldcg[0]+phot.CoRoT_LDC_weights[ccc,2]*ldcr[0]+phot.CoRoT_LDC_weights[ccc,3]*ldci[0]+phot.CoRoT_LDC_weights[ccc,4]*ldcz[0],phot.CoRoT_LDC_weights[ccc,0]*ldcu[1]+phot.CoRoT_LDC_weights[ccc,1]*ldcg[1]+phot.CoRoT_LDC_weights[ccc,2]*ldcr[1]+phot.CoRoT_LDC_weights[ccc,3]*ldci[1]+phot.CoRoT_LDC_weights[ccc,4]*ldcz[1]]
    else: 
        return LDCs(logg, teff, z, photband)

    
