"""
Module for the analysis of Markov Chains.
"""

from math import *
import pickle
import numpy as np
import gzip
import scipy
import string
import glob

from .Objects_MCMC import Chain
from .tools import state_constructor
from .priors import prior_constructor, compute_priors
from .. import *
# from .. import resultpath
# from .. import Msun, Rsun, Mjup, Rjup, G, au, Rsun2Rjup, Msun2Mjup, Mearth
# from .. import tools
from ..paths import resultpath
from .. import AstroClasses, ObjectBuilder
from ..limbdarkening import get_LD
from .. import photometry as phot
from .. import isochrones as iso


class VDchain(object):
    def __init__(self, vd, beta, target, runid, filename=None):
        self._value_dict = vd
        self.beta = beta
        self.target = target
        self.runid = runid.rstrip('_').lstrip('_')
        self.filename = filename

    def get_value_dict(self):
        return self._value_dict


def make_solution_file(chain, pastisfile=None, BI=0.2, best=True,
                       outputname=False):
    """
    Build the .solution file of a given chain.

    Parameters
    ----------
    chain: file, string, dict, Chain instance, VDchain, or list instance
        The chain(s) for which the solution file will be created. It can either
        be a file object containing the chain (i.e. a .chain file), a string
        with the full path of that file, a Chain object, a dictionary
        containing the traces of the parameters, or a list of VDchain instances.
        If it is a dictionary, then the pastisfile parameter in mandatory.

    pastisfile: string, optional
        The name of the .pastis configuration file used to construct the
        .solution file. If it is not given, the configuration dictionaries
        are taken as those in the Chain attribute _inputdicts. If the first
        parameter is a dicitonary, then this parameter in not optional.

    Other parameters
    ----------------
    BI: float
        The fraction of the chain to be discarded for the computation of the
        state to be written in the solution file. This is used if the traces
        contain the burn-in period.

    best: boolean
        Decides whether the state with the best highest logL is used or a median
        of all values.

    outputname: boolean
        If true, the name of the solution file is returned.

    See also
    --------
    get_best_values, get_median_values
    """

    if isinstance(chain, file):
        C = pickle.load(chain)

        fout = open(chain.name.replace('.chain', '.solution'), 'w')

    elif isinstance(chain, str):
        # Read chain from file
        f = open(chain, 'r')
        C = pickle.load(f)
        f.close()

        if '.chain' in chain:
            fout = open(chain.replace('.chain', '.solution'), 'w')
        elif '.mcmc' in chain:
            fout = open(chain.replace('.mcmc', '.solution'), 'w')
        else:
            print('Unrecognized file format.')
            return

    elif isinstance(chain, Chain) or isinstance(chain, dict) or isinstance(
            chain, VDchain) or isinstance(chain, list) or isinstance(chain,
                                                                     n.array):
        C = chain
        foutname = pastisfile.replace('.pastis', '.solution')
        foutname = foutname.replace('configfiles', 'resultfiles')
        fout = open(foutname, 'w')

    else:
        print('First parameter must either be a filename, a file object, or'
              ' a MCMC Chain object.')
        return

    # Read input dicts
    if pastisfile is not None:
        f = open(pastisfile, 'r')
        dd = pickle.load(f)
        f.close()

    elif isinstance(C, Chain):
        dd = C._inputdicts

    else:
        print('Error! If the input is a dictionary, or a VDchain, then a'
              ' .pastis file (configuration file) must be given!')
        return

    # Get best values or median values
    if best:
        bestdict = get_best_values(C, BI)
    else:
        bestdict = get_median_values(C, BI)

    input_dict = dd[1].copy()

    for obj in input_dict.keys():
        for param in input_dict[obj].keys():
            if isinstance(input_dict[obj][param], str):
                continue

            try:
                # Copy jump parameter values from bestdict to input_dict
                input_dict[obj][param][0] = bestdict[obj][param][0]
            except:
                # Less fixed parameters as they are.
                pass

    pickle.dump([dd[0], input_dict, dd[2], dd[3]], fout)
    fout.close()

    print('Saved solution file to %s' % fout.name)

    if outputname:
        return fout.name
    else:
        return


def confidence_intervals(C, q=0.6827, hdi=None, nbins=50, burnin=0.0,
                         percentile=True, cumulative=False, difftol=0.3,
                         pmsym='+/-', report='mean', output=False):
    """
    Computes confidence intervals for all parameters in a given chain, and
    print it to the standard output.

    Parameters
    ----------
    C, Chain instance, VDchain instance or dictionary
       The chain for which the confidence intervals are to be computed. An
       actual Chain instance, or a dictionary containing all the traces for
       the chain can be given.

    Other parameters
    ----------------
    q, float
        Fraction of entire distribution to use for the computation of the
        vals (default: 0.6827, i.e. 68.27% confidence interval, 1-sigma).

    hdi, float or list
        The fraction of the mass of the distribution to include in the
        Highest Denstity Interval. It can be a float or an iterable object
        with a series of values. If None, no hdi is computed
    
    nbins, int
        The number of bins used for the computation. If negative, the
        computation is done using the absolute value of nbins for the mode,
        but using the cumulative distribution for the confidence intervals.

    burnin, float
        The fraction to be discarded from the beginning of each trace. Useful
        if the trace contain the burn-in period.

    percentile, boolean
        Defines whether confidence intervals are obtained using the percentile
        function instead of interpolating histogram bins.

    cumulative, boolean
        Defines if instead of returning median and confidence limits, return
        q*100% cumulative value.
        
    difftol, float
        Maximum fractional difference between upper and lower limits above which
        both values are printed. If the difference is below this value, print
        only the largest of both values.

    pmsym, string
        Symbol to use for +/-. It can be used to include latex code

    report, string
        Controls which statistics is reported. It can be 'mean', 'median',
        'mode' or 'map'.

    output, boolean
        Decides whether the value and confidence intervals are returned in a
        dictionary or not.
        
    """
    if isinstance(C, Chain):
        vd = C.get_value_dict()
    elif isinstance(C, VDchain):
        vd = C.get_value_dict()
    elif isinstance(C, dict):
        vd = C
    else:
        print('Warning! Instance not recognized. Assuming it is VDchain '
              'instance.')
        try:
            vd = C.get_value_dict()
        except NameError:
            raise NameError('Failed to load posterior sample. Check input.')

    # Try to compute index of MAP
    try:
        mapindex = n.argmax(vd['posterior'])
    except KeyError:
        mapindex = n.nan
        if report is 'map':
            print('Posterior not given, will return posterior mode.')
            report = 'mode'

    if cumulative:
        # Minimum and maximum prob
        qmin = 0
        qmax = q

    else:
        # Minimum and maximum prob
        qmin = (1 - q) * 0.5
        qmax = (1 + q) * 0.5

    # Initialise output dictionary.
    outputdict = {}

    # Find index after burn in.
    istart = int(n.round(len(vd[vd.keys()[0]]) * burnin))

    for param in n.sort(vd.keys()):

        # Get median and mode, and dispersion
        statdict = {'median': n.median(vd[param][istart:]),
                    'mean': n.mean(vd[param][istart:]),
                    'sigma': n.std(vd[param][istart:]),
                    'min': n.min(vd[param][istart:]),
                    'max': n.max(vd[param][istart:])
                    }

        # Try to add MAP value
        try:
            statdict['map'] = vd[param][mapindex]
        except IndexError:
            statdict['map'] = n.nan

        # Compute histogram for all cases that require it.
        if report is 'mode' or not percentile or hdi is not None:
            # Compute histogram
            m, bins = n.histogram(vd[param][istart:], nbins, normed=True)
            x = bins[:-1] + 0.5 * n.diff(bins)

            modevalue = x[n.argmax(m)]

            statdict['mode'] = modevalue

        ###
        # Find confidence intervals
        ###
        if percentile:
            # Method independent of bin size.
            lower_limit = n.percentile(vd[param][istart:], 100 * qmin)
            upper_limit = n.percentile(vd[param][istart:], 100 * qmax)

        else:
            # Original method
            ci = m.astype(float).cumsum() / n.sum(m)

            if not cumulative:
                imin1 = n.argwhere(n.less(ci, qmin))

                if len(imin1) >= 1:
                    imin1 = float(imin1.max())
                else:
                    imin1 = 0.0

                imin2 = imin1 + 1.0
                lower_limit = scipy.interp(qmin, [ci[imin1], ci[imin2]],
                                           [x[imin1], x[imin2]])
            else:
                lower_limit = 0.0

            imax1 = float(n.argwhere(n.less(ci, qmax)).max())
            imax2 = imax1 + 1.0
            upper_limit = scipy.interp(qmax, [ci[imax1], ci[imax2]],
                                       [x[imax1], x[imax2]])

        try:
            reportvalue = statdict[report]
        except KeyError:
            raise NameError('report statistics not recognised.')

        hc = upper_limit - reportvalue
        lc = reportvalue - lower_limit

        statdict['{:.2f}-th percentile'.format(qmin * 100)] = lower_limit
        statdict['loweruncertainty'] = lc
        statdict['{:.2f}-th percentile'.format(qmax * 100)] = upper_limit
        statdict['upperuncertainty'] = hc

        # ## Trick to show good number of decimal places.
        """
        # Express hc and lc in scientific notation.
        s = '{:e}'.format(min(hc, lc))
        # Find the exponent in scientific notation
        a = s[s.find('e') + 1:]
        """
        # Express hc and lc in scientific notation.
        s = '{:e}'.format(statdict['sigma'])
        # Find the exponent in scientific notation
        a = s[s.find('e') + 1:]  # if negative, the number of decimals is equal to the exponent + 1
        if a[0] == '-':
            ndecimal = int(a[1:]) + 1
        # Else, just one decimal place.
        else:
            ndecimal = 1
        # print('sigma {}: {}, {}'.format(param, statdict['sigma'], ndecimal))

        statdict['ndec'] = int(ndecimal)
        # If upper and lower errors do not differ by more than difftol*100%
        # only report one of them
        maxerr = max(hc, lc)

        if cumulative:

            formatdict = {'param': param, 'value': upper_limit,
                          'ndec': int(ndecimal), 'pmsym': pmsym,
                          'err': maxerr, 'q': qmax * 100}

            formatdict.update(statdict)

            print('{param} ({q}-th percentile): {value:.{ndec}f} '
                  '(mean = {mean:.{ndec}f}, '
                  'median = {median:.{ndec}f}, '
                  'map = {map:.{ndec}f}, '
                  'min = {min:.{ndec}f}, '
                  'max = {max:.{ndec}f}, '
                  'sigma = {sigma:.1e})'.format(**formatdict))

        else:
            if abs(hc - lc) / maxerr < difftol:

                formatdict = {'param': param, 'value': reportvalue,
                              'ndec': int(ndecimal), 'pmsym': pmsym,
                              'err': maxerr}

                formatdict.update(statdict)
                print('{param}: {value:.{ndec}f} {pmsym} {err:.1e} '
                      '(mean = {mean:.{ndec}f}, '
                      'median = {median:.{ndec}f}, '
                      'map = {map:.{ndec}f}, '
                      'min = {min:.{ndec}f}, '
                      'max = {max:.{ndec}f}, '
                      'sigma = {sigma:.1e})'.format(**formatdict))

            else:
                hcr = n.round(hc, ndecimal)
                lcr = n.round(lc, ndecimal)
                formatdict = {'param': param, 'value': reportvalue,
                              'ndec': int(ndecimal), 'pmsym': pmsym,
                              'err+': hcr, 'err-': lcr}

                formatdict.update(statdict)

                print('{param}: {value:.{ndec}f} +{err+:.1e} -{err-:.1e} '
                      '(mean = {mean:.{ndec}f}, '
                      'median = {median:.{ndec}f}, '
                      'map = {map:.{ndec}f}, '
                      'min = {min:.{ndec}f}, '
                      'max = {max:.{ndec}f}, '
                      'sigma = {sigma:.1e})'.format(**formatdict))

        # Compute HDI
        if hdi is not None:
            hdi = n.atleast_1d(hdi)

            for hdii in hdi:
                hdints = compute_hdi(bins, m, q=hdii)

                statdict['{}%-HDI'.format(hdii*100)] = hdints

                hdistr = '{0:.1f}% HDI: '.format(hdii * 1e2)
                for jj, interval in enumerate(hdints):
                    if jj > 0:
                        hdistr += ' U '

                    if n.isnan(interval[0]):
                        continue
                    hdistr += (
                        '[{0:.' + str(ndecimal + 1) + 'f}, {1:.' + str(
                            ndecimal + 1) + 'f}]').format(*interval)
                print(hdistr)

        if output:
            outputdict[param] = statdict

    if output:
        return outputdict
    else:
        return


def compute_derived_params(vd, pastisfile, sampling=1, **kwargs):
    getmags = kwargs.pop('getmags', False)
    bandpasses = kwargs.pop('bandpasses', ['Johnson-V', ])

    # Read dictionaries from pastisfile
    f = open(pastisfile, 'r')
    dd = pickle.load(f)
    f.close()

    objs = dd[1].keys()

    N = len(vd[vd.keys()[0]][::sampling])

    # Define dictionary that will contain derived parameters
    # and auxiliary dictionary vvd
    dvd = {}
    vvd = {}

    # Get names of stars in systems and jump parameters
    jump_params = []
    fixed_params = []

    # Get filters of light curves
    filters = []
    for key in dd[2].keys():
        try:
            filters.append(dd[2][key]['filter'])
        except:
            continue

    sorted_objs = []
    for oo in objs:
        if 'PlanSys' in oo:
            sorted_objs.append(dd[1][oo]['star1'])

        for param in dd[1][oo].keys():
            if type(dd[1][oo][param]) == list:
                # Skip parameters with None
                if dd[1][oo][param][0] is None:
                    continue

                if dd[1][oo][param][1] >= 1:
                    # This is a jump parameter
                    jump_params.append(oo + '_' + param)
                elif dd[1][oo][param][1] == 0:
                    fixed_params.append(oo + '_' + param)

    # Add remaining objects
    for oo in objs:
        if oo not in sorted_objs:
            sorted_objs.append(oo)

    # Add fixed and jump parameters to auxiliary dicitionary dvd
    for pp in fixed_params:
        dvd[pp] = n.zeros(N) + dd[1][pp.split('_')[0]][pp.split('_')[1]][0]

    for pp in jump_params:
        dvd[pp] = vd[pp][::sampling]

    for pp in ['logL', 'posterior']:
        dvd[pp] = vd[pp][::sampling]

    # Search for Host and Target stars
    indhost = map(string.find, sorted_objs, ['Host'] * len(sorted_objs))
    indtarg = map(string.find, sorted_objs, ['Target'] * len(sorted_objs))

    # Target and/or Host stars are present
    if any(n.array(indhost) != -1) or any(n.array(indtarg) != -1):
        # Write initialization condition
        condtarg = 'verticeT' in iso.__dict__
    else:
        # If no Target or Host, target tracks initialization considered done
        condtarg = True

    # Search for Blended stars
    indblend = map(string.find, sorted_objs, ['Blend'] * len(sorted_objs))

    # Blend stars are present
    if any(n.array(indblend) != -1):
        # Write initialization condition
        condblend = 'verticesY' in iso.__dict__
    else:
        condblend = True

    # Initialize only if needed
    if 'AMz' not in phot.__dict__ or not condtarg or not condblend:
        datadict, a = readdata(dd[2])
        initialize(dd[0], datadict, dd[1])

    reload(AstroClasses)
    reload(ObjectBuilder)

    # Compute parameters for stars in systems
    for oo in sorted_objs:
        if 'Target' in oo or 'Host' in oo or 'Blend' in oo:

            print('%s: computing stellar parameters.' % oo)
            # Compute M, R
            if 'Host' in oo:
                host = True
                y = dvd[oo + '_dens']
            elif 'Target' in oo:
                host = False
                y = dvd[oo + '_logg']

            # Include arrays to hold derived parameters in dvd
            dvd[oo + '_mact'] = n.zeros(N)
            dvd[oo + '_R'] = n.zeros(N)
            dvd[oo + '_logL'] = n.zeros(N)
            # Theoretical limb darkening
            for filt in filters:
                dvd[oo + '_ua_theo_' + filt] = n.zeros(N)
                dvd[oo + '_ub_theo_' + filt] = n.zeros(N)

            if 'Target' in oo:
                dvd[oo + '_logage'] = n.zeros(N)
                dvd[oo + '_dens'] = n.zeros(N)

            elif 'Blend' in oo:
                dvd[oo + '_teff'] = n.zeros(N)
                dvd[oo + '_dens'] = n.zeros(N)
                dvd[oo + '_logg'] = n.zeros(N)

            elif 'Host' in oo:
                dvd[oo + '_logage'] = n.zeros(N)
                dvd[oo + '_logg'] = n.zeros(N)

            for i in range(N):

                if (i + 1) % 100 == 0:
                    print('Step %d out of %d' % ((i + 1), N))
                if 'Target' in oo or 'Host' in oo:
                    # Interpolate track
                    try:
                        Mact, logL, logage = iso.get_stellarparams_target(
                            dvd[oo + '_z'][i], y[i],
                            log10(dvd[oo + '_teff'][i]),
                            planethost=host
                        )
                    except:
                        continue

                    # Save partial results
                    dvd[oo + '_logage'][i] = logage

                elif 'Blend' in oo:

                    # Interpolate track
                    logT, logg, logL, Mact = \
                        iso.get_stellarparams_target(
                            dvd[oo + '_z'][i],
                            dvd[oo + '_logage'][i],
                            dvd[oo + '_minit'][i]
                        )

                    Teff = 10 ** logT

                    # Save partial results
                    dvd[oo + '_teff'][i] = Teff
                    dvd[oo + '_logg'][i] = logg

                # Save mass and luminosity
                dvd[oo + '_mact'][i] = Mact
                dvd[oo + '_logL'][i] = logL

                # Compute stellar radius
                L = 10 ** dvd[oo + '_logL'][i]
                try:
                    alphaS = dvd[oo + '_alphaS'][i]
                except KeyError:
                    alphaS = 0.0
                R = n.sqrt(
                    L * (5777.0 / dvd[oo + '_teff'][i]) ** 4.0 / (1 - alphaS))
                dvd[oo + '_R'][i] = R

                # Compute logg for host star
                if 'Host' in oo:
                    mm = Mact * Msun
                    rr = R * Rsun

                    g_SI = G * mm / rr ** 2
                    g_cgs = g_SI * 1e2

                    dvd[oo + '_logg'][i] = n.log10(g_cgs)

                # Compute stellar density for blends and targets
                elif 'Blend' in oo or 'Target' in oo:
                    dvd[oo + '_dens'][i] = Mact / R ** 3

                # Compute theoretical limb darkening
                for filt in filters:
                    dvd[oo + '_ua_theo_' + filt][i], \
                    dvd[oo + '_ub_theo_' + filt][i] = \
                        get_LD(dvd[oo + '_teff'][i], dvd[oo + '_logg'][i],
                               dvd[oo + '_z'][i], filt)

                if getmags:
                    # CONTINUE HERE!!!
                    logg = dvd[oo + '_logg'][i]

        elif 'Fit' in oo or 'Planet' in oo:

            print('%s: checking if it belongs to a system.' % oo)
            insystem = False
            # Check if it belongs to a system
            for ss in objs:
                if 'PlanSys' in ss:
                    # Iterate over all planets in the system
                    for pp in dd[1][ss].keys():
                        if oo == dd[1][ss][pp]:
                            # Fitobs belong to a system!
                            insystem = True
                            print('%s belongs to system %s.' % (oo, ss))
                            # Get key for the star in the system
                            oos = dd[1][ss]['star1']
                            print('Setting host star to %s' % oos)

            if insystem:

                Ms = dvd[oos + '_mact'] * Msun
                Rs = dvd[oos + '_R']

                k_ms = dvd[oo + '_K1'] * 1e3
                k_kms = dvd[oo + '_K1']
                per = dvd[oo + '_P']
                per_s = dvd[oo + '_P'] * 86400.
                ecc = dvd[oo + '_ecc']
                omega = dvd[oo + '_omega'] * pi / 180.0

                #
                # Separate cases where mass ratio is given/fitted  or not.
                #
                if oo + '_q' in dvd.keys():
                    #
                    # Mass ratio q has been fitted or fixed
                    #
                    q = dvd[oo + '_q']

                elif not oo + '_q' in dvd.keys():
                    #
                    # Mass ratio q will be obtained self-consistently from fit.
                    #

                    if oo + '_incl' in dvd.keys() and not oo + '_b' in dvd.keys():
                        # Inclination is being fitted
                        y = dvd[oo + '_incl']
                        use_b = False

                    elif oo + '_b' in dvd.keys():
                        # Impact parameter b is being fitted
                        y = dvd[oo + '_b']
                        use_b = True

                    else:
                        # No inclination nor b, fix incl = 90
                        dvd[oo + '_incl'] = per * 0.0 + 90.0
                        y = dvd[oo + '_incl']
                        use_b = False

                    mact = n.zeros(N)

                    for i in range(N):
                        try:
                            # Tolerance is 1/1000th of an Earth mass.
                            macti = tools.iterative_mass(k_kms[i], per[i],
                                                         dvd[oos + '_mact'][i],
                                                         ecc[i], y[i], Rs=Rs[i],
                                                         omega=
                                                         dvd[oo + '_omega'][i],
                                                         use_b=use_b,
                                                         tol=0.001 * Mearth)

                        except RuntimeError:
                            print(k_kms[i], per[i], dvd[oos + '_mact'][i],
                                  ecc[i], y[i], Rs[i], dvd[oo + '_omega'][i])
                            continue

                        mact[i] = macti

                    q = mact / dvd[oos + '_mact'][i]

                    dvd[oo + '_mact'] = mact
                    dvd[oo + '_Mp'] = mact * Msun2Mjup

                ##
                ## Compute semi-major axis of relative orbit
                ##
                a3 = 0.25 * G / (pi * pi) * (1 + q) * Ms * per_s ** 2
                dvd[oo + '_a'] = a3 ** (1. / 3.) / au

                ##
                ## Compute a/R
                ##
                dvd[oo + '_ar'] = (a3 ** (1. / 3.)) / (Rs * Rsun)

                # Compute true anomaly at transit center
                nu0 = pi / 2.0 - omega

                # Compute star-planet distance at transit center
                rtran = dvd[oo + '_ar'] * (1 - ecc ** 2) / (
                    1 + ecc * n.cos(nu0))

                ##
                ## Compute b if inclination was a jump parameter
                ##
                if oo + '_incl' in dvd.keys() and not oo + '_b' in dvd.keys():

                    ## Compute b
                    dvd[oo + '_b'] = rtran * n.cos(
                        dvd[oo + '_incl'] * pi / 180.0)
                    incl = dvd[oo + '_incl'] * pi / 180.0

                ##
                ## OR compute inclination if the parameter was b
                ##
                elif oo + '_b' in dvd.keys():
                    cosi = dvd[oo + '_b'] / rtran
                    incl = n.arccos(cosi)  # in radians
                    dvd[oo + '_incl'] = incl * 180.0 / pi  # in degrees

                if oo + '_q' in dvd.keys():
                    ##
                    ## Compute companion mass for case where q was given.
                    ##
                    mact = k_ms * (per_s * Ms ** 2 * (1 + q) ** 2 / (
                        2 * pi * G)) ** (1. / 3.) * \
                           n.sqrt(1 - ecc ** 2) / n.sin(incl)

                    dvd[oo + '_mact'] = mact / Msun
                    dvd[oo + '_Mp'] = mact / Mjup

                ##
                ## Compute impact parameter at secondary eclipe
                ##
                nuS = 3 * pi / 2.0 - omega
                rsec = dvd[oo + '_ar'] * (1 - ecc ** 2) / (1 + ecc * n.cos(nuS))

                dvd[oo + '_bsec'] = rsec * n.cos(dvd[oo + '_incl'] * pi / 180.0)

                ##
                ## Compute radius, density, and surface gravity
                ##
                if oo + '_kr' in dvd:
                    kr = dvd[oo + '_kr']

                    dvd[oo + '_Rp'] = Rs * kr * Rsun2Rjup
                    dvd[oo + '_dens'] = dvd[oo + '_Mp'] / dvd[
                                                              oo + '_Rp'] ** 3  # Jupiter units
                    dvd[oo + '_logg'] = n.log10(
                        1e2 * G * dvd[oo + '_Mp'] * Mjup * (
                            dvd[oo + '_Rp'] * Rjup) ** (-2.0))

                    ##
                    ## Compute transit duration
                    ## (uses eq. 7 of Tingley & Sacket 2005)
                    Rp_m = dvd[oo + '_Rp'] * Rjup
                    Rs_m = Rs * Rsun

                    # Probability of secondary eclipse
                    dvd[oo + '_probsec'] = (1 + dvd[oo + '_kr']) / dvd[
                        oo + '_bsec']

                    # Use star-planet distance at transit center in meters
                    # rtran_m = rtran * Rs * Rsun
                    # r0_m = Rsun * a_rsun * ( 1 - ecc**2 ) / ( 1 + ecc * n.cos(nu0) )

                    # Compute geometrical factor Z
                    # Z = n.sqrt( 1 - (rtran_m * n.cos(incl))**2 / (Rs_m + Rp_m)**2 )
                    Z = n.sqrt(1 - dvd[oo + '_b'] ** 2 / (1 + kr) ** 2)

                    # Duration
                    D = 2 * Z * Rs_m * (1 + kr) * n.sqrt(1 - ecc ** 2) * \
                        (1 + ecc * n.cos(nu0)) ** (-1) * \
                        (per_s / (2 * pi * G * Ms * (1 + q))) ** (1. / 3.)

                    dvd[oo + '_duration'] = D / 3600.0  # In hours

                ##
                ## Compute equilibrium temperature
                ## (assuming black body and isotropic planetary emission)
                if oo + '_albedo2' in dvd:
                    albedo = dvd[oo + '_albedo2']
                    teff = dvd[oos + '_teff']
                    a_rsun = dvd[oo + '_a'] * au / Rsun

                    teq = teff * (1 - albedo) ** (1. / 4.) * n.sqrt(
                        0.5 * Rs / a_rsun)
                    dvd[oo + '_teq'] = teq


            else:
                print('%s not in system: not yet implemented.' % oo)

                ##
                ## Compute stellar density
                ##

                continue

                # Time of passage through periastron
    return dvd


def identify_objects(keys):
    objects = []
    for kk in keys:
        obj = kk.split('_')[0]
        if not obj in objects:
            objects.append(obj)
        else:
            continue

    return objects


def get_best_values(C, BI=0.0, bestparam='posterior'):
    """
    Find the parameter values with the best logL of a given chain.

    Parameters
    ----------
    C: Chain, VDchain, list instance or dictionary instance
        The chain or list of chains for which the best values will be obtained.

    BI: float
        The fraction of the chain to be discarded. This is used if the chain
        contains the burn-in period.

    bestparam: str
        Specifies the parameter that plays the role of the merit function.
        It can be "posterior" or "logL".

    Return
    ------
    bestvalues: dict
        A dictionary with the best values for each parameter, both jumping and
        fixed ones.
    """
    bestvalues = {}
    singlechain = True

    if isinstance(C, Chain):
        vd = C.get_value_dict()
        N = C.N

        try:
            if bestparam == 'posterior':
                y = C.get_posterior()
            elif bestpaaram == 'logL':
                y = C.get_logL()
        except:
            print('Warning! {s} not found in chain.'.format(s=bestparam))

    elif isinstance(C, VDchain):
        vd = C.get_value_dict()
        N = len(vd[vd.keys()[0]])

        try:
            y = vd[bestparam]
        except KeyError:
            print('Warning! {s} not found in chain.'.format(s=bestparam))

    elif isinstance(C, dict):
        vd = C.copy()
        N = len(vd[vd.keys()[0]])

        try:
            y = vd[bestparam]
        except KeyError:
            print('Warning! {s} not found in chain.'.format(s=bestparam))

    elif isinstance(C, list) or isinstance(C, n.array):
        singlechain = False

    else:
        print('Input not recognized.')
        return

    if not singlechain:
        # Find for which chain the maximum posterior occurs
        for i, vd in enumerate(C):
            if i == 0:
                maxpost = n.max(vd._value_dict[bestparam])
                maxind = 0
                continue

            if n.max(vd._value_dict[bestparam]) > maxpost:
                maxpost = n.max(vd._value_dict[bestparam])
                maxind = i

        # Get chain containing best value
        vd = C[maxind]._value_dict
        N = len(vd[vd.keys()[0]])
        y = vd[bestparam]

    indi = n.int(n.round(BI * N))
    indmax = n.argmax(y[indi:])

    for key in vd.keys():
        skey = key.split('_')
        if skey[0] == 'logL' or skey[0] == 'posterior': continue
        if skey[0] in bestvalues:
            pass
        else:
            bestvalues[skey[0]] = {}

        # For each parameter, take the value for the highest logL.
        bvi = vd[key][indi:][indmax]
        bestvalues[skey[0]][skey[1]] = [bvi, ]

    return bestvalues


def get_median_values(C, BI=0.0):
    """
    Find the median values of all the parameters of a given chain.

    Parameters
    ----------
    C: Chain, VDchain, or dictionary instance
        The chain for which the median values will be obtained.

    BI: float
        The fraction of the chain to be discarded. This is used if the chain
        contains the burn-in period.

    Return
    ------
    medianvalues: dict
        A dictionary with the median values for each jumping parameter and the
        values of the fixed ones.
    """
    medianvalues = {}

    if isinstance(C, Chain):
        vd = C.get_value_dict()
        N = C.N
        try:
            y = C.get_posterior()
        except:
            print('Warning! Using logL instead of log(prior*likelihood)')
            y = C.get_logL()

    elif isinstance(C, VDchain):
        vd = C.get_value_dict()
        N = len(vd[vd.keys()[0]])
        try:
            y = vd['posterior']
        except KeyError:
            print('Warning! Using logL instead of log(prior*likelihood)')
            y = vd['logL']

    elif isinstance(C, dict):
        vd = C.copy()
        N = len(vd[vd.keys()[0]])
        try:
            y = vd['posterior']
        except KeyError:
            print('Warning! Using logL instead of log(prior*likelihood)')
            y = vd['logL']


    else:
        print('Input not recognized.')
        return

    indi = n.round(BI * N)
    indmax = n.argmax(y[indi:])

    for key in vd.keys():
        skey = key.split('_')
        if skey[0] == 'logL' or skey[0] == 'posterior': continue
        if skey[0] in medianvalues:
            pass
        else:
            medianvalues[skey[0]] = {}

        # For each parameter, take the median of the chain.
        bvi = vd[key][indi:]
        medianvalues[skey[0]][skey[1]] = [n.median(bvi), ]

    return medianvalues


def read_chains(target, runid, beta=None):
    """
    Load chain results of a given run of a target.

    Parameters
    ----------
    
    target: string
        the name of the target whose chains you want to load. The chains will
        be read from the resultfiles directory corresponding to this target.

    runid: string
        a string corresponding to the runID. In practice, all .mcmc files
        in the target results directory whose name contains the runID will be
        read.

    Returns
    -------

    filenames: list
        a list containing the names of the .mcmc files that were read.

    vds: list
        a list of the value dictionaries of each chain.
        See Chain.get_value_dict()

    Other parameters
    ----------------
    beta: float
        the value of the tempering parameters of the chains to read.
        If not given, the value will be obtained from the runid, if possible.
    """
    chaindir = os.path.join(resultpath, target)
    flist = glob.glob(os.path.join(chaindir, target + '*' + runid + '_Beta*'))

    if beta is not None:
        try:
            betai = float(beta)
            betastr = 'Beta' + str(betai)
        except TypeError:
            raise TypeError('beta parameter cannot be converted to float')
    else:
        # Beta is not given, it will be obtained from the filenames
        betastr = ''
        pass

    filenames = []
    vds = []
    for ffile in n.sort(flist):

        # Filter broken files
        if ffile == 'KOI-189_ALL_FitTargetStar_20120722T18:36:18.mcmc':
            continue

        # Select only files that contain the runid and the beta parameter
        if (runid in ffile) and (betastr in ffile) and (
                    ffile.endswith('mcmc') or ffile.endswith('mcmc.gz')):
            print(ffile)
            if ffile.endswith('mcmc'):
                f = open(os.path.join(chaindir, ffile), 'r')
            elif ffile.endswith('mcmc.gz'):
                f = gzip.open(os.path.join(chaindir, ffile), 'r')
            vd = pickle.load(f)
            f.close()

            # IF beta was not given, take it from filename
            if beta == None:
                indi = ffile.find('Beta')
                if indi != -1:
                    indf = ffile.find('_', indi + 4)
                    if indf != -1:
                        betai = float(ffile[indi + 4: indf])
                    else:
                        betai = float(ffile[indi + 4:])
                else:
                    betai = n.nan

            ## Create VDchain instance
            vdC = VDchain(vd, betai, target, runid, filename=f.name)

            vds.append(vdC)
            filenames.append(ffile)
    return filenames, vds


def read_chains2(self, flist, beta=[], target=[], runid=[]):
    """
    Auxiliary function for widget
    """
    '''
    if beta != []:
	try:
	    betai = float(beta)
	    betastr = 'Beta'+str(betai)
	except TypeError:
	    raise TypeError('beta parameter cannot be converted to float')
    else:
	# Beta is not given, it will be obtained from the filenames
	betastr = ''
	    
    '''
    import Tkinter
    import MeterBar

    betai = []
    betastr = []
    for b in beta:
        betai.append(float(b))
        betastr.append('Beta' + str(b))
    filenames = []
    vds = []
    meterbar = Tkinter.Toplevel()
    m = MeterBar.Meter(meterbar, relief='ridge', bd=3)
    m.pack(fill='x')
    m.set(0.0, 'Loading resultfiles...')
    for i, ffilepath in enumerate(n.sort(flist)):

        ffile = os.path.split(ffilepath)[1]

        # Select only files that contain the runid and the beta parameter
        if (betastr[i] in ffile) and (
                    ffile.endswith('mcmc') or ffile.endswith('mcmc.gz')):
            # print(ffile)
            if ffile.endswith('mcmc'):
                f = open(ffilepath, 'r')
            elif ffile.endswith('mcmc.gz'):
                f = gzip.open(ffilepath, 'r')
            vd = pickle.load(f)
            f.close()

            ## Create VDchain instance
            vdC = VDchain(vd, betai[i], target[i], runid[i], filename=f.name)

            vds.append(vdC)
            filenames.append(ffile)
            m.set(float(i + 1) / len(flist))
    m.set(1., 'Loading resultfiles finished')
    meterbar.destroy()

    return filenames, vds


def merge_chains(vds, BI, CL, N=1e4, pickrandom=False,
                 beta=None):
    """
    Construct a merged value dictionary based on a series of value dictionaries.

    The merged dictionary is constructed by randomly selecting N samples
    from each of the given value dictionries.

    Parameter
    ---------
    vds: list
        A list containing the value dictionaries to be merged, or VDchain
	instances.

    BI: float or list
        The fraction of the chains to be discarded before constructing the
	merged dicionary.
	An iterable object containing the fraction to be discarded for each
	individual chain to be merged can also be given. In this case, the
	number of elements in BI must equal that in vds.

    CL: float or list
        Correlation length: the number of samples to thin the chain (i.e., the
        merged chain will have one element every CL).
	An iterable object containing the correlation length for each
	individual chain to be merged can also be given. In this case, the
	number of elements in CL must equal that in vds.

    N: int
        The number of samples to randomly take from each chain (only if
	pickrandom is True)

    beta: float
        If given, and vds is a list of VDchains, then only those chains with the
        correct value of beta are used.

    Returns
    -------
    vdc: dict or VDchain
        A dictionary containing the constructed traces for each parameter in
	vds, or the corresponding VDchain instance.
    """

    if isinstance(vds[0], dict):
        vdd = vds[0].copy()

    else:
        if not isinstance(vds[0], VDchain):
            # Assume its a VDchain from before reloading the module....
            print(
                'Warning! Instance not recognized. Assuming it is VDchain instance.')
        vdd = vds[0].get_value_dict()

    # Create dictionary that will contain random samples
    if pickrandom:
        vdc = dict((kk, n.zeros(N * len(vds))) for kk in vdd.keys())
    else:
        vdc = dict((kk, n.array([])) for kk in vdd.keys())

    ## Create list of starting indexes
    try:
        iter(BI)
    except TypeError:
        # BI is a float; convert to a list
        BI = [BI] * len(vds)

    # Check if BI contains the same number of elements as vds
    if len(BI) != len(vds):
        raise TypeError('BI must be either a float or contain the same number\
of elements as the input list')

    start_index = n.zeros(len(BI), 'int')
    for i in range(len(vds)):
        if isinstance(vds[i], dict):
            vdd = vds[i].copy()
        else:
            if not isinstance(vds[i], VDchain):
                # Assume its a VDchain from before reloading the module. Print
                # warning
                print(
                    'Warning! Instance not recognized. Assuming it is VDchain instance.')
            vdd = vds[i].get_value_dict()

        if BI[i] >= 1:
            raise ValueError('All BIs must be less than 1.')
        start_index[i] = BI[i] * len(vdd[vdd.keys()[0]])
    ##

    ## Create list of correlation lengths
    try:
        iter(CL)
    except TypeError:
        # CL is a float; convert to a list
        CL = [CL] * len(vds)

    inds0 = n.arange(N * len(vds)).astype(int)
    n.random.shuffle(inds0)

    for i, vdi in enumerate(vds):

        if isinstance(vdi, dict):
            vd = vdi.copy()

        else:
            if not isinstance(vdi, VDchain):
                # Assume its a VDchain from before reloading the module. Print
                # warning
                print(
                    'Warning! Instance not recognized. Assuming it is VDchain instance.')

            if beta is not None:
                # Check beta
                if vdi.beta != beta:
                    print(
                        'Warning! Tempering parameter of chain %s is not the requested one' %
                        os.path.split(vdi.filename)[-1])
                    continue

            vd = vdi.get_value_dict()

        if pickrandom:
            inds = n.arange(len(vd[vd.keys()[0]][start_index[i]:]))
            if len(inds) < N:
                raise ValueError(
                    'All elements of the value dicts must have at least N elements.')
            n.random.shuffle(inds)

            for kk in vd.keys():
                vdc[kk][inds0[i * N:(i + 1) * N]] = vd[kk][start_index[i]:][
                                                        inds][:N]

        else:
            for kk in vd.keys():
                vdc[kk] = n.concatenate(
                    (vdc[kk], vd[kk][start_index[i]::CL[i]]))

    if not isinstance(vdi, dict):
        vdc = VDchain(vdc, beta, target=vdi.target, runid=vdi.runid,
                      filename=None)

    return vdc


def corrlength(x, step=1, BI=0.2, BO=1.0, widget=False, verbose=True,
               plot=False, **kwargs):
    """
    Computes the correlation length of a given trace of a Parameter

    The value of the shift for which the correlation reaches 1/e is printed.

    Parameters
    ----------
    x: ndarray
        An array containing the elements of the trace.

    step: int
        The number of elements to shift x at each step of the computation.

    BI: float
        The fraction of the chain to be discarded for the computation of the
        correlation length. This is used if the traces
        contain the burn-in period.

    BO: float
        The fraction up to which the chain is considered. To use in combination
        of BI. E.g.: a BI of 0.2 and a BO of 0.6 imply that the correlation is
        computed over 40 % of the chain.
    
    widget: boolean
        To display or not the printed information

    plot: boolean
        To perform a plot of the correlation length

    Other parameters
    ----------------
    xlabel: str
        Label for x axis of plot (default: 'Step')

    ylabel: str
        Label for y axis of plot (default: 'Correlation')

    title: str
        Title for plot (default: '')
        
    circular: bool
        Assumes the time series have simmetry, so that when shifted, the last
        part goes to the beginning. This saves time.

    stop: bool
        Defines whether the computation stops when the correlation has fallen
        below 1/e or if it continues to the end.
        
    Returns
    -------
    shifts: ndarray
        The shift values for which the computation was done.
	
    corr: ndarray
        The correlation value for each step.

    corrlength: int
        The number of steps after which the correlation has fallen below 1/e.

    Notes
    -----
    To reduce computation of useless values, the algorithm is stopped when the
    last 500 values of the correlation are below 0.2.
    """

    xlabel = kwargs.pop('xlabel', 'Step')
    ylabel = kwargs.pop('ylabel', 'Correlation')
    title = kwargs.pop('title', '')
    circular = kwargs.pop('circular', True)
    stop = kwargs.pop('stop', True)

    if widget == True:
        verbose = False

    indstart = n.round(len(x) * BI)
    indend = n.round(len(x) * BO)
    x = x[indstart: indend]

    ## Reduce mean
    x = x - x.mean()

    if circular:
        ## Compute values of unchanged chain only once
        xmean = x.mean()
        x2mean = (x ** 2).mean()
        xx2 = xmean ** 2.0
        den = x2mean - xmean ** 2.0
    #
    shifts = n.zeros(len(x) / float(step))
    corr = n.zeros(len(x) / float(step))
    for j in range(len(corr)):
        #
        if circular:
            xs = cshift(x, j * step)
            corr[j] = ((x * xs).mean() - xx2) / den
        else:
            ## Compute for each iteration
            if j == 0:
                corr[j] = 1
            else:
                xs = x[:-j * step]
                xmean = x[j * step:].mean()
                x2mean = (x[j * step:] ** 2).mean()
                xx2 = xmean ** 2.0
                den = x2mean - xmean ** 2.0

                corr[j] = ((x[j * step:] * xs).mean() - xx2) / den

        shifts[j] = j * step

        #
        if (j + 1) % 100 == 0 and verbose:
            print('Step {} out of a maximum of {}'.format(j+1, len(corr)))
            os.sys.stdout.flush()

        if j > 500.0 and stop:
            # Stop iteration if last 500 points are below 1/e
            if n.alltrue(n.less(corr[j + 1 - 500: j + 1], 1.0 / e)):
                break

    shifts, corr = shifts[:j + 1], corr[:j + 1]
    try:
        corrlength = n.min(n.compress(corr < 1.0 / e, shifts))
    except ValueError:
        corrlength = len(x)
        print('Error! Correlation length not found.')
    else:
        if verbose: print('Correlation drops to 1/e after %d steps' %
                          corrlength
                          )

    if plot:
        import pylab as p

        fig1 = p.figure()
        ax = fig1.add_subplot(111)
        ax.plot(shifts, corr)
        ax.set_xlabel(xlabel, fontsize=16)
        ax.set_ylabel(ylabel, fontsize=16)
        ax.set_title(title, fontsize=16)
        ax.axhline(1 / e, ls=':', color='0.5')
        ax.axvline(corrlength, ls=':', color='0.5')
        p.draw()
        print(ax)

    return shifts, corr, corrlength


def corrlength2(x, step=1, BI=0.2, widget=False):
    """
    Computes the correlation length of a given trace of a Parameter

    The value of the shift for which the correlation reaches 1/e is printed.

    Parameters
    ----------
    x: ndarray
        An array containing the elements of the trace.

    step: int
        The number of elements to shift x at each step of the computation.

    BI: float
        The fraction of the chain to be discarded for the computation of the
        correlation length. This is used if the traces
        contain the burn-in period.
    
    widget: boolean
        To display or not the printed information

    Returns
    -------
    shifts: ndarray
        The shift values for which the computation was done.
	
    corr: ndarray
        The correlation value for each step.

    corrlength: int
        The number of steps after which the correlation has fallen below 1/e.

    Notes
    -----
    To reduce computation of useless values, the algorithm is stopped when the
    last 500 values of the correlation are below 0.2.
    """
    x = x - x.mean()

    indstart = n.round(len(x) * BI)

    x = x[indstart:]

    xmean = x.mean()
    x2mean = (x ** 2).mean()
    xx2 = xmean ** 2.0
    den = x2mean - xmean ** 2.0
    #
    shifts = []  # n.zeros(len(x)/float(step))
    corr = []  # n.zeros(len(x)/float(step))
    j = 0
    while j < len(x) / float(step):
        # for j in range(len(corr)):
        #
        xs = cshift(x, j * step)

        shifts.append(j * step)
        corr.append(((x * xs).mean() - xx2) / den)
        #
        if (j + 1) % 100 == 0 and not widget:
            print('Step {} out of a maximum of {}'.format(j+1, len(corr)))
            os.sys.stdout.flush()

        if j > 500.0:
            # Stop iteration if last 500 points are below 0.2
            if n.alltrue(n.less(corr[j + 1 - 500: j + 1], 0.2)): break

        j += 1

    shifts = n.array(shifts)
    corr = n.array(corr)

    # shifts, corr = shifts[:j + 1], corr[:j + 1]   # solve the problem with append ?
    corrlength = n.min(n.compress(corr < 1.0 / e, shifts))
    if not widget: print('Correlation drops to 1/e after %d steps' %
                         n.min(n.compress(corr < 1.0 / e, shifts))
                         )

    return shifts, corr, corrlength


def cshift(x, j):
    j = j % len(x)
    return n.concatenate((x[j:], x[:j]))


def corrlength_multichain(vds, step=1, BI=0.2, plot=False, widget=False,
                          verbose=True, plotCL=False, **kwargs):
    """
    Compute correlation lenght for all parameters of a multichain

    Parameters
    ----------
    vds: list
        a list of the value dictionaries of each chain. See Chain.get_value_dict()
    
    step: int
        The number of elements to shift x at each step of the computation.

    BI: float or list
        The fraction of the chain to be discarded for the computation of the
        correlation length.
	An iterable object containing the fraction to be discarded for each
	individual chain to be merged can also be given. In this case, the
	number of elements in BI must equal that in vds.
        
    plot: bool
        plot the chain correlation lenght for each parameter
        
    plotCL: bool
        plot the chain correlation curve  for each parameter of each chain.

    widget, boolean
        option that displays a widget status bar

    Other parameters
    ----------------
    The remaining keyword arguments are passed to corrlength function
    
    Returns
    -------
    corrlen: dict
        A dictionary containing the correlation lenght for all chains for each 
        parameter (key)
 
    """
    if widget:
        verbose = False

    import Tkinter
    import MeterBar
    ## Check if it is a list of VDchain instances or of dictionaries
    if isinstance(vds[0], VDchain):
        vdd = vds[0].get_value_dict()
    elif isinstance(vds[0], dict):
        vdd = vds[0].copy()
    else:
        # Assume its a VDchain from before reloading the module....
        print('Warning! Class VDchain have changed.')
        vdd = vds[0].get_value_dict()

    ## Create list of starting indexes
    try:
        iter(BI)
    except TypeError:
        # BI is a float; convert to a list
        BI = [BI] * len(vds)

    # Check if BI contains the same number of elements as vds
    if len(BI) != len(vds):
        raise TypeError('BI must be either a float or contain the same number\
of elements as the input list')

    # Create output dictionary
    corrlen = {}
    if widget:
        meterbar = Tkinter.Toplevel()
        m = MeterBar.Meter(meterbar, relief='ridge', bd=3)
        m.pack(fill='x')
        m.set(0.0, 'Computing correlation lenght ...')

    for i, parameter in enumerate(vdd.keys()):  # parameter boucle
        if parameter in ['logL', 'posterior']: continue

        if not widget:
            print('Computing corrlength of ' + parameter)

        corrlengthi = []
        for chain in range(len(vds)):  # chain boucle

            if isinstance(vds[chain], VDchain):
                vdd = vds[chain].get_value_dict()
            elif isinstance(vds[chain], dict):
                vdd = vds[chain].copy()
            else:
                # Assume its a VDchain from before reloading the module....
                print('Warning! Class VDchain might have changed.')
                vdd = vds[chain].get_value_dict()

            shifts, corr, corrlengthc = corrlength(vdd[parameter],
                                                   step=step,
                                                   BI=BI[chain],
                                                   widget=widget,
                                                   verbose=verbose,
                                                   plot=plotCL,
                                                   title=parameter,
                                                   **kwargs)

            corrlengthi.append(corrlengthc)
            if widget: m.set(
                float(i) / len(vdd.keys()) + float(chain + 1) / len(
                    vdd.keys()) / len(vds), parameter)

        if widget: m.set(float(i + 1) / len(vdd.keys()), parameter)

        corrlen[(parameter)] = corrlengthi

    if widget:
        m.set(1., 'Correlation lenght computed')
        meterbar.destroy()

    if plot:
        import pylab as p

        for parameter in corrlen.keys():
            p.figure()
            p.plot(corrlen[(parameter)], 'o')
            p.xlabel('chain index')
            p.ylabel(parameter)

    if plotCL:
        import pylab as p

        p.show()

    return corrlen


def corrlength_multichain2(vds, step=1, BI=0.2, plot=False, widget=False):
    """
    Compute correlation lenght for all parameters of a multichain

    Parameters
    ----------
    vds: list
        a list of the value dictionaries of each chain. See Chain.get_value_dict()
    
    step: int
        The number of elements to shift x at each step of the computation.

    BI: float
        The fraction of the chain to be discarded for the computation of the
        correlation length. This is used if the traces
	contain the burn-in period.
        
    plot: bool
        plot the chain correlation lenght for each parameter
        
    widget, boolean
        option that displays a widget status bar

    Returns
    -------
    corrlen: dict
        A dictionary containing the correlation lenght for all chains for each 
        parameter (key)
 
    """
    import Tkinter
    import MeterBar
    ## Check if it is a list of VDchain instances or of dictionaries
    if isinstance(vds[0], VDchain):
        vdd = vds[0].get_value_dict()
    elif isinstance(vds[0], dict):
        vdd = vds[0].copy()
    else:
        print(type(vds[0]))

    ## Create list of starting indexes
    try:
        iter(BI)
    except TypeError:
        # BI is a float; convert to a list
        BI = [BI] * len(vds)

    # Create output dictionary
    corrlen = {}
    if widget:
        meterbar = Tkinter.Toplevel()
        m = MeterBar.Meter(meterbar, relief='ridge', bd=3)
        m.pack(fill='x')
        m.set(0.0, 'Computing correlation lenght ...')

    for i, parameter in enumerate(vdd.keys()):  # parameter boucle
        if parameter == 'logL': continue
        if not widget:
            print('--------------------')
            print('    ' + parameter)
            print('--------------------')
        corrlengthi = []
        for chain in range(len(vds)):  # chain boucle

            if isinstance(vds[chain], VDchain):
                vdd = vds[chain].get_value_dict()
            elif isinstance(vds[chain], dict):
                vdd = vds[chain].copy()
            else:
                if not widget: print(type(vds[chain]))

            shifts, corr, corrlengthc = corrlength2(vdd[parameter],
                                                    step=step, BI=BI[chain],
                                                    widget=widget)
            corrlengthi.append(corrlengthc)
            if widget: m.set(
                float(i) / len(vdd.keys()) + float(chain + 1) / len(
                    vdd.keys()) / len(vds), parameter)

        if widget: m.set(float(i + 1) / len(vdd.keys()), parameter)

        corrlen[(parameter)] = corrlengthi

    if widget:
        m.set(1., 'Correlation lenght computed')
        meterbar.destroy()

    if plot:
        import pylab as p

        for parameter in corrlen.keys():
            p.figure()
            p.plot(corrlen[(parameter)], 'o')
            p.xlabel('chain index')
            p.ylabel(parameter)

    return corrlen


def corrlenchain(corrlen):
    """
    Compute the maximum correlation length among all parameters of a chain

    Parameters
    ----------
    corrlen: dict
        A dictionary containing the correlation lenght for all chains for each 
        parameter (key)

    Returns
    -------
    cl: list
        Maximum correlation length among all parameters of a chain.
        To use as CL keyword in merge_chains.
    """

    # get correlation length values and put into an array
    acorr = n.array(corrlen.values())

    # reshape (all correlation length values for each chain instead of for each parameter)
    acorr2 = acorr.reshape(len(corrlen.keys()), len(corrlen[corrlen.keys()[0]]))

    # get maximum correlation length in each chain
    cl = n.max(acorr2, axis=0)

    return cl


def checkpriors(vdc, pastisfile, **kargs):
    """
    Plot histogram of each parameter of the merged chain together with 
    its prior taken from a .pastis configuration file 

    Parameter
    ---------
    vdc: dict
        The merged chain, output of the merge_chains function.

    pastisfile: string
        The name of the .pastis configuration file 
        
    **kargs
        Parameters for the hist function

    """

    f = open(pastisfile, 'r')
    dd = pickle.load(f)
    f.close()

    priordict = prior_constructor(dd[1], dd[3])

    if isinstance(vdc, VDchain):
        vdd = vdc.get_value_dict()
    elif isinstance(vdc, dict):
        vdd = vdc.copy()
    else:
        print(type(vdc))

    ### PLOTS
    import pylab as p

    for parameter in vdd.keys():
        if parameter != 'logL' and parameter != 'posterior':
            p.figure()
            pdf, bins, patches = p.hist(vdd[parameter], normed=True, **kargs)
            p.xlabel(parameter)
            xmin = n.min(vdd[parameter])
            xmax = n.max(vdd[parameter])
            x = n.arange(xmin, xmax, (xmax - xmin) / 1000.)
            y = priordict[parameter].pdf(x)
            p.plot(x, y, 'r')

    return


def register_chain(vd, beta=None, target=None, runid=None, force=False):
    """
    Saves a given chain to be used by the ModelComparison module. The chain
    is pickled to a file and the beta value and name of this file are recorded
    in a file under /resultfile/target/

    Parameters
    ----------
    vd: dict or VDchain instance
        The chain to be pickled. It must be already thinned, merged, and free
        of burn-in period.

    Other paramters
    ---------------
    If vd is a dictionary, additional parameters are needed to create the file
    that contains the chain.

    - beta, target, runid: str

    force: boolean
        If true, any existing chain file will be overwritten. Else, an error
        is risen.

    """
    if isinstance(vd, dict):
        if n.any([beta == None, target == None, runid == None]):
            raise Error(
                'If input is a dictionary, parameters beta, target, and runid must be specified.')

        vd = VDchain(vd, beta, target, runid)

    elif isinstance(vd, VDchain):
        beta = vd.beta
        target = vd.target
        runid = vd.runid

    # Prepare output file
    fname = os.path.join(resultpath, target,
                         '%s_%s_Beta%.6f_mergedchain.dat' % (
                             target, runid, beta)
                         )
    if os.path.exists(fname) and not force:
        raise Warning('File %s exists!' % fname)

    fout = open(fname, 'w')
    pickle.dump(vd, fout)
    fout.close()

    # Prepare list file
    fnamel = os.path.join(resultpath, target,
                          '%s_%s_mergedchain.lst' % (target, runid)
                          )
    foutl = open(fnamel, 'a')
    foutl.write('%.6f\t%s\n' % (beta, os.path.split(fname)[-1])
                )
    foutl.close()
    return


def get_multichain(vdchain, beta=None):
    """
    Gets multichain dictionary from list of VDchain objects
    
    Parameters
    ----------
    vdchain: VDchain
        list of VDchain objects

    Returns
    -------
    vds: dict
        multichain dictionary
    """

    vds = []

    for chain in n.arange(len(vdchain)):
        if vdchain[chain].beta == beta or beta == None:
            vds.append(vdchain[chain]._value_dict)
        else:
            continue

    return vds


def print_param(vds, kk, BI=0.5):
    """
    Print the median and std of parameter kk for all chains in vds.
    Also print name of chains and median logL
    """

    for i in range(len(vds)):
        vd = vds[i].get_value_dict()
        fname = os.path.split(vds[i].filename)[-1]
        N = len(vd['logL'])

        error = n.std(vd[kk][BI * N:])
        ## Trick to show good number of decimal places
        s = '%e' % error
        a = s[s.find(
            'e') + 1:]  # Contains power of range of axis in scientific notation
        if a[0] == '-':
            ndecimal = int(a[1:]) + 1
        else:
            ndecimal = 1

        if i == 0:
            print('Filename\tlog(L)\t%s\tsigma(%s)' % (kk, kk))

        fmtstr = '%s\t%.1f\t%.' + str(ndecimal) + 'f\t%.' + str(ndecimal) + 'f'
        print(fmtstr % (fname, n.median(vd['logL'][BI * N:]),
                        n.median(vd[kk][BI * N:]), error)
              )
    return


def get_Mtier(AR, Period, solarunits=True):
    """
    Return the M**(1/3)/R value. Valid only for transiting planets.

    Parameters:
    -----------

    AR : list, array or float
        System scale a/R*.

    Period : list, array or float
        Orbital period of the planet, in days.

    solarunits : boolean
        If true returns the value expressed in solar units.

    Returns:
    --------

    Mtier : list, array or float
        return the M**(1./3.)/R value
    """
    from .. import G, Msun, Rsun
    # if len(AR) <> len(Period) :
    #	raise ValueError('Arguments must have the same dimension')
    cte = (4. * n.pi ** 2 / G) ** (1. / 3.)
    Mtier = cte * AR / (Period * 86400.) ** (2. / 3.)
    return Mtier / Msun ** (1. / 3.) * Rsun


def get_ImpactParameter(AR, inc, ecc, omega, component='primary'):
    """
    Return the impact parameter of the transit.

    Parameters:
    -----------

    AR : float, list or array
        System scale a/R*

    inc : float, list or array
        orbital inclination [deg]

    ecc : float, list or array
        orbital eccentricity

    omega : float, list or array
        argument of periastron [deg]

    component : string
        either 'primary' or 'secondary'

    Returns:
    --------

    b : float, list or array
        Impact parameter

    """
    import numpy as n

    if component == 'primary':
        b = AR * n.cos(inc / 180. * n.pi) * (1. - ecc ** 2) / (
            1. + ecc * n.sin(omega / 180. * n.pi))
    elif component == 'secondary':
        b = AR * n.cos(inc / 180. * n.pi) * (1. - ecc ** 2) / (
            1. - ecc * n.sin(omega / 180. * n.pi))
    # else : raise 'InputError' : 'component must be either primary or secondary'
    return b


def get_TransitDuration(P, AR, k, inc, ecc, omega, component='primary'):
    """
    Return the transit duration of the transit.

    Parameters:
    -----------
    P : float, list or array
        Orbital period [d]

    AR : float, list or array
        System scale a/R*

    k : float, list or array
        Rdius ratio Rp/R*

    inc : float, list or array
        orbital inclination [deg]

    ecc : float, list or array
        orbital eccentricity

    omega : float, list or array
        argument of periastron [deg]

    component : string
        either 'primary' or 'secondary'

    Returns:
    --------

    T14 : float, list or array
        Total transit duration [hours]

    """
    import numpy as n

    b = get_ImpactParameter(AR, inc, ecc, omega, component=component)

    if component == 'primary':
        T14 = P / n.pi * n.arcsin(n.sqrt((1. + k ** 2) - b ** 2) / AR / n.sin(
            inc / 180. * n.pi)) * n.sqrt(1. - ecc ** 2) / (
                  1. + ecc * n.sin(omega * n.pi / 180.))
    elif component == 'secondary':
        T14 = P / n.pi * n.arcsin(n.sqrt((1. + k ** 2) - b ** 2) / AR / n.sin(
            inc / 180. * n.pi)) * n.sqrt(1. - ecc ** 2) / (
                  1. - ecc * n.sin(omega * n.pi / 180.))
    # else : raise 'InputError' : 'component must be either primary or secondary'
    return T14 * 24.


def get_Tp(P, T0, ecc, omega):
    """
    Return the epoch of periastron from the other orbital parameters.

    Parameters:
    -----------

    P : float, list or array
        Orbital period, in days

    T0 : float, list or array
        Transit epoch

    ecc : float, list or array
        Orbital eccentricity

    omega : float, list or array
        Argument of periastron, in degrees

    Returns:
    --------

    Tp : float, list or array
        Epoch of periastron
    """
    import numpy as n

    E_0 = n.arctan2(n.sqrt(1. - ecc ** 2) * n.cos(omega * n.pi / 180.),
                    n.sin(omega * n.pi / 180.) + ecc)
    Tp = T0 - P / (2. * n.pi) * (E_0 - ecc * n.sin(E_0))
    return Tp


###
# CONVERGENCE DIAGNOSTICS
###

def gelmanrubin(vds, BI=0.2, BO=1.0, thinning=1, qs=[0.9, 0.95, 0.99]):
    """
    Compute Gelman & Rubin statistics for a multi-chain run.
    
    Code copied from 
    http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_introbayes_sect008.htm#statug.introbayes.bayesgelm
    
    :return dict Rcs: dictionary of Rc statistic for each parameter.
    """
    if isinstance(vds[0], dict):
        vdi = vds[0].copy()
    elif isinstance(vds[0], VDchain):
        vdi = vds[0].get_value_dict()
    else:
        print(vds[0].__class__)
        
    keys = list(vdi.keys())
    start_index = n.int(BI * len(vdi[keys[0]]))
    end_index = n.int(BO * len(vdi[keys[0]]))

    Ws = {}
    Bs = {}
    Vs = {}
    MUs = {}
    DFs = {}
    PSRF = {}
    Rcs = {}

    for j, kk in enumerate(vdi.keys()):

        values = []
        for i, vd in enumerate(vds):

            if isinstance(vd, dict):
                pass
            elif isinstance(vd, VDchain):
                vd = vd.get_value_dict()
            
            try:
                values.append(vd[kk][start_index: end_index: thinning])
            except KeyError:
                print('Parameter {} missing from some chain'.format(kk))
                continue

        values = n.array(values)
        nn = values.shape[1] # Number of steps
        m = values.shape[0] # Number of walkers

        # Compute within-chain variance
        sm2 = n.var(values, axis=1, ddof=1)  # Variance for each chain
        W = n.mean(sm2)  # Mean variance over all chains
        Ws[kk] = W

        # Compute between-chain variance
        thetamean = n.mean(values, axis=1)  # Mean for each chain
        B = nn * n.var(thetamean, ddof=1)  # Variance of means, multiplied by nn
        Bs[kk] = B

        # Estimate mean (mu) using all chains
        mu = n.mean(values)
        MUs[kk] = mu

        # Estimate variance by weighted average of B and W (eq.3 Gelman & Rubin)
        sig2 = (nn - 1.0) / nn * W + B / nn

        # Rc statistic
        Rcs[kk] = np.sqrt(sig2/W)
        
        # The parameter distribution can be approximated by a Student's t
        # distribution (Gelman & Rubin) with scale sqrt(V):
        # print kk, sig, B/(nn*m)
        V = sig2 + B / (nn * m)
        Vs[kk] = V

        # and degrees of freedom df = 2*V**2/var(V) (see eq. 4 G & R)
        varV = ((nn - 1.0) / nn) ** 2.0 * (1.0 / m) * n.var(sm2) + \
               ((m + 1.0) / (nn * m)) ** 2.0 * (2.0 / (m - 1.0)) * B ** 2.0 + \
               2 * (m + 1.0) * (nn - 1.0) / (m * nn ** 2.0) * \
               (1.0 * nn / m) * (n.cov(sm2, thetamean ** 2.0)[1, 0] - \
                                 2 * mu * n.cov(sm2, thetamean)[1, 0])

        df = 2 * V ** 2.0 / varV
        DFs[kk] = df

        psr = n.sqrt((V / W) * df / (df - 2.0))
        PSRF[kk] = psr

        ## Compute degrees of freedom for F distribution for PSRF
        # see sect 3.5 and 3.7 of G & R
        dfn = m - 1
        dfd = 2.0 * W ** 2 / (n.var(sm2) / m)

        ## Compute 90%, 95% and 99% percentiles for this distribution
        qq = scipy.stats.f.ppf(qs, dfn, dfd)

        lims = n.sqrt(
            ((nn - 1.) / nn + (m + 1.) / (nn * m) * qq) * df / (df - 2.0))
        
        if j == 0:
            print('Parameter\tPSRF\t{:d}%\t{:d}%\t{:d}%'
                  ''.format(*[int(q*100) for q in qs]))
        print('%s\t%.3f\t%.3f\t%.3f\t%.3f' % (kk, psr, lims[0],
                                              lims[1], lims[2]))

    return Rcs


def geweke(X, first=0.1, last=0.5, Nsamples=20, BI=0):
    """
    Compute the Geweke diagnostic of covergence.

    X can be the trace of a given parameter, or an entire chain, in which case
    the diagnostic is computed for all jump parameters.

    Nsamples is the number of samples to take from last part of chain.
    """
    if isinstance(X, objMCMC.Chain):

        # Prepare dictionary for Z
        Z = {}

        # Iterate over all parameters
        for par in X._labeldict.keys():
            if X._labeldict[par].jump:
                xx = X.get_trace(par)

                ## Compute Geweke Z

                # Get first part of the chain after BI
                N = len(xx[BI:])
                xxc = xx[BI:]

                fxx = xxc[:int(first * N)]
                lxx = xxc[int(first * N):]

                Z[par] = []
                for jj in range(Nsamples):
                    # Out of the last part of the chain lxx,
                    # get the specified part
                    a = range(len(lxx))
                    random.shuffle(a)
                    lxxc = lxx[a[:int(len(lxx) * last)]]

                    # Compute mean and variance of each part
                    z = (n.mean(fxx) - n.mean(lxxc)) / sqrt(
                        n.var(fxx) + n.var(lxxc))
                    Z[par].append(z)

    else:
        N = len(X[BI:])
        xxc = xx[BI:]
    return Z


# To find the BI of a given chain
def find_BI(vds, samplesize=0.05, endsample=0.1, backwards=True,
            tolerance=2.0, checkvariance=True, correlen=1, param='logL',
            sigmaclip=False, nsigma=5,
            nitersigmaclip=3, verbose=False):
    """
    Find the BurnIn of a given chain or a group of chains by comparing the
    mean and, possibly, also the variance of the trace of parameter param
    throughtout the chain to that at the end of the chain.
    
    Parameters
    ----------
    vds: list, dict, or VDchain instance.
        The list of chains, or individual chain for which the BI will be
        computed. The list can contain dict instances with the traces of the
        parameters, or VDchain instances. All chains must contain param  as a 
        parameter.

    Other parameters
    ----------------
    samplesize: float.
       Size of the sample used to compare to the end sample, at the end of
       the chain. It must be expressed as a fraction of the total length.

    endsample: float.
        Size of the end part of the logL trace used to compare with the rest 
        of the chain

    tolerance: float.
        Number of standard deviations to use as tolerance.

    checkvariance: boolean.
        Determines if a test on the variance of the samples is also performed.

    correlen: float or list
        Correlation length: the number of samples to thin the chain.
        An iterable object containing the correlation length for each
        individual chain for which the BI will be measured can also be given.
        In this case, the number of elements in correlen must equal that in vds.

    param: string.
        Parameter key used to compute burn-in. Default: log(Likelihood)

    sigmaclip: boolean.
        Determines if samples throughout the chain are sigma-clipped before
        computing its variance and mean.

    nsigma: int.
        Number of sigmas to use for sigma clipping algorithm. Default: 5.

    nitersigmaclip: int.
        Number of times the sigma-clipping algorithm is iterated.
    """

    try:
        vds = list(vds)
    except TypeError:
        pass

    # Create list of correlation lengths
    try:
        iter(correlen)
    except TypeError:
        # correlen is a float; convert to a list
        correlen = [correlen] * len(vds)

    # Check if single chain is given.
    if not isinstance(vds, list):
        vds = [vds, ]

    # Define list to contain results
    zlist = []
    zzlist = []
    BIfrac = []
    for j, vd in enumerate(vds):

        if isinstance(vd, dict):
            y = vd[param][::correlen[j]]

        elif isinstance(vd, VDchain):
            y = vd.get_value_dict()[param][::correlen[j]]

        else:
            print('Warning! VDchain class might have changed.')
            y = vd.get_value_dict()[param][::correlen[j]]

        N = len(y)

        # Select end sample
        yf = y[-N * endsample:]

        if sigmaclip:
            # Sigma clip end sample
            yf, ccs = tools.sigma_clipping(yf, nsigma,
                                           niter=nitersigmaclip
                                           )

        mean_yf = n.mean(yf)
        var_yf = n.var(yf, ddof=1)

        # Explore the chain
        Nsample = N * samplesize

        # Define list to contain results
        zlisti = []
        zzlisti = []
        if not backwards:
            ei = 0
            ef = samplesize
            # Forward
            while ef < 1 - endsample:

                # Select sample to compare to end sample
                ys = y[ei * N: ef * N]

                if sigmaclip:
                    # Sigma clips sample
                    ys, ccs = tools.sigma_clipping(ys, nsigma,
                                                   niter=nitersigmaclip
                                                   )

                var_ys = n.var(ys, ddof=1)

                # Compute Geweke statistics (modified: instead of using the
                # variance of sample, which can be extremely large and hinder
                # a correct estimation of BI, we use twice the variance of the
                # end sample.
                z = (n.mean(ys) - mean_yf) / n.sqrt(2 * var_yf)

                # The variance must also be approximately the same
                # If normal, these variables should be Chi2(N-1), 
                # so their variance is 2*(N - 1), where N is the size of the 
                # sample.
                zz = (var_ys - var_yf) / n.sqrt(
                    2 * (len(ys) - 1) + 2 * (len(yf) - 1))

                # Ad hoc condition on maximum of sample
                """
                print n.max(ys), mean_yf, tolerance*n.sqrt(var_yf)
                condAH = n.max(ys) >= (mean_yf - tolerance*n.sqrt(var_yf))
                """
                zlisti.append(z)
                zzlisti.append(n.abs(zz))

                if checkvariance:
                    accept = n.abs(z) < tolerance and n.abs(zz) < tolerance
                else:
                    accept = n.abs(z) < tolerance
                if accept:
                    print('Chain %d: BI set to '
                          '%.2f' % (j, ef + 0.5 * samplesize))
                    BIfrac.append(min(ef + 0.5 * samplesize, 1.0))
                    zlist.append(zlisti)
                    zzlist.append(zzlisti)
                    break

                ei = ef
                ef = ei + samplesize

            if ef >= 1 - endsample:
                # If BI was not found, print warning and set BI = 1
                print('Burn In fraction not found. Has the chain converged?')
                BIfrac.append(1)
                zlist.append(zlisti)
                zzlist.append(zzlisti)

        else:
            ef = 1.0 - endsample
            ei = ef - samplesize
            # Backwards
            while ei >= -1e-8:
                # Select sample to compare to end sample
                ys = y[ei * N: ef * N]

                if sigmaclip:
                    # Sigma clips sample
                    ys, ccs = tools.sigma_clipping(ys, nsigma,
                                                   niter=nitersigmaclip
                                                   )

                var_ys = n.var(ys, ddof=1)

                # Compute Geweke statistics
                z = (n.mean(ys) - mean_yf) / n.sqrt(2 * var_yf)

                # The variance must also be approximately the same
                # If normal, these variables should be Chi2(N-1), 
                # so their variance is 2*(N - 1), where N is the size of the 
                # sample.
                zz = (var_ys - var_yf) / n.sqrt(
                    2 * (len(ys) - 1) + 2 * (len(yf) - 1))

                # Ad hoc condition on maximum of sample
                condAH = n.max(ys) >= (mean_yf - tolerance * n.sqrt(var_yf))

                """
                ## Compute KS test on the two samples
                zz, probKS = scipy.stats.ks_2samp(yf, ys)
                """

                if verbose:
                    # print z, zz, probKS
                    print(ei, ef, mean_yf, var_yf, n.mean(ys), var_ys, condAH)
                zlisti.append(n.abs(z))
                zzlisti.append(n.abs(zz))

                if checkvariance:
                    accept = n.abs(z) > tolerance or n.abs(
                        zz) > tolerance or -condAH
                else:
                    accept = n.abs(z) > tolerance or -condAH

                if accept:
                    print('Chain %d: BI set to'
                          ' %.2f' % (j, ef + 0.5 * samplesize))
                    BIfrac.append(min(ef + 1.5 * samplesize, 1.0))
                    zlist.append(zlisti)
                    zzlist.append(zzlisti)
                    break

                ef = ei
                ei = ef - samplesize

            if ei < -1e-8:
                # If BI was not found, means that BI = 0!
                print('Burn In fraction not found. Setting BI = 0')
                BIfrac.append(0)
                zlist.append(zlisti)
                zzlist.append(zzlisti)

    return zlist, zzlist, BIfrac


def select_best_chains(vds, BI, CL, nmin=100, tolerance=2.0,
                       KStolerance=0.1 / 1e2,
                       param='posterior', param2=None, fnames=None):
    """
    Identify the chains having a significant worse logL than the best one.
    
    Parameters
    ----------
    vds: list, dict, or VDchain instance.
        The list of chains to compare. The list can contain dict instances
        with the traces of the parameters, or VDchain instances.
        All chains must contain param  as a parameter.

    BI: float or list
        The fraction of the chains to be discarded before comparing
	An iterable object containing the fraction to be discarded for each
	individual chain to be can be given. In this case, the
	number of elements in BI must equal that in vds. Otherwise a single
        BI for all chains is used

    CL: float or list
        The correlation length of the chains to be compared. Chains will be
        thinned by this number before comparing.
	An iterable object containing the fraction to be discarded for each
	individual chain to be can be given. In this case, the
	number of elements in CL must equal that in vds. Otherwise a single
        CL for all chains is used
        
    Other parameters
    ----------------
    tolerance: float.
        Number of standard deviations to use as tolerance.

    KStolerance: float.
        Tolerance for KS test.

    param: string.
        Parameter key used to compute burn-in. Default: log(Likelihood)

    fnames: list.
        List of names of chains to compare.
    """
    try:
        vds = list(vds)
    except TypeError:
        pass

    # Check if single chain is given.
    if not isinstance(vds, list):
        vds = [vds, ]

    # Define list to contain results
    ys = []
    mediany = []
    sigmay = []
    pars = []
    ireject = []
    iaccept = []
    ## Create list of starting indexes
    try:
        iter(BI)
    except TypeError:
        print('Warning!')
        # BI is a float; convert to a list
        BI = [BI] * len(vds)

    ## Create list of correlation lengths
    try:
        iter(CL)
    except TypeError:
        # CL is a float; convert to a list
        CL = [CL] * len(vds)

    # Check if BI and CL contain the same number of elements as vds
    for v, namev in zip((BI, CL), ('BI', 'CL')):
        if len(v) != len(vds):
            raise TypeError(
                '%s must be either a float or contain the same number of elements as the input list' % namev)

    ## Prepare lists
    accepted_chains = []
    rejected_chains = []
    accepted_ind = []

    if fnames is not None:
        accepted_fnames = []
        rejected_fnames = []

    ## Iterate over all chains
    for j, vd in enumerate(vds):

        if isinstance(vd, dict):
            y = vd[param]
            if param2 is not None:
                y2 = vd[param2]

        elif isinstance(vd, VDchain):
            y = vd._value_dict[param]
            if param2 is not None:
                y2 = vd._value_dict[param2]
        else:
            print('Warning! VDchain class might have changed.')
            y = vd._value_dict[param]
            if param2 is not None:
                y2 = vd._value_dict[param2]


                ## Compute median and standard deviation
        yc = y[len(y) * BI[j]:: CL[j]]
        if param2 is not None:
            y2c = y2[len(y2) * BI[j]:: CL[j]]

        if len(yc) <= nmin:
            print('Chain %d: Less than %d steps' % (j, nmin))
            ireject.append(j)
            rejected_chains.append(vd)
            if fnames is not None:
                rejected_fnames.append(fnames[j])
            continue
        else:
            iaccept.append(j)

        median_y = n.median(yc)
        sigma_y = n.std(yc, ddof=1)

        ## Add array to list
        ys.append(yc)
        if param2 != None:
            pars.append(y2c)
        ##
        mediany.append(median_y)
        sigmay.append(sigma_y)

    y0 = ys[n.argmax(mediany)]
    my0 = mediany[n.argmax(mediany)]
    sy0 = sigmay[n.argmax(mediany)]

    probsKS = []
    for i in range(len(mediany)):

        if fnames is not None:
            fname = fnames[iaccept[i]]
        else:
            fname = ''

        ## Compute difference between chain i and best
        z = (mediany[i] - my0) / n.sqrt(sy0 ** 2 + sigmay[i] ** 2)

        ## Compute KS test on the two samples
        zz, probKS = scipy.stats.ks_2samp(y0, ys[i])
        probsKS.append(probKS)

        rejectstr = ''
        if abs(z) > tolerance or probKS < KStolerance:
            rejectstr = '*'
            rejected_chains.append(vds[iaccept[i]])
            if fnames is not None:
                rejected_fnames.append(fname)

        else:
            accepted_chains.append(vds[iaccept[i]])
            accepted_ind.append(iaccept[i])
            if fnames is not None:
                accepted_fnames.append(fname)

        if z == 0:
            rejectstr = '!'

        if param2 != None:

            print('%s%s : %.3f\t%.3f\t%.3f\t%.3e\t%.3f\t%d' % (rejectstr,
                                                               fname,
                                                               mediany[i],
                                                               sigmay[i], z,
                                                               probKS,
                                                               n.median(
                                                                   pars[i]),
                                                               len(pars[i])
                                                               )
                  )
        else:
            print('%s%s : %.3f\t%.3f\t%.3f\t%.3e' % (rejectstr,
                                                     fname, mediany[i],
                                                     sigmay[i], z,
                                                     probKS,
                                                     )
                  )

    if fnames is not None:
        return accepted_chains, accepted_fnames, rejected_chains, rejected_fnames, accepted_ind
    else:
        return accepted_chains, rejected_chains, accepted_ind


def get_priors_from_value_dict(vds, pastisfile):
    """
    Get the priors for all steps of a chain represented by vd.
    Adds a key to vd dictionary called 'prior'

    Parameters
    ----------
    vds: list
        List containing dictionary, or VDchain instances of the chains for
        which to compute the priors.

    pastisfile: str, or file instance
        The file containing the configuration for the run for which to compute
        the priors.        
    """

    # Read configuration file
    if isinstance(pastisfile, file):
        f = pastisfile
    else:
        f = open(pastisfile, 'r')
    dd = pickle.load(f)
    f.close()

    # Construct priordict from configuration file
    priordict = prior_constructor(dd[1], dd[3])

    newvds = []
    for vd in vds:
        print(vd)

        # Get value dict
        if isinstance(vd, VDchain):
            vdd = vd.get_value_dict().copy()
        elif isinstance(vd, dict):
            vdd = vd.copy()
        else:
            print(
                'Warning! Instance not recognized. Assuming it is VDchain instance.')
            vdd = vd.get_value_dict().copy()

        N = len(vdd[vdd.keys()[0]])
        vdd['prior'] = n.zeros(N)

        ### Iterate over all elements in chain.
        for i in range(N):

            if i % 1e4 == 0:
                print('Step %d out of %d' % (i, N))
            # Based on input_dict from pastis file, construct a new dictionary with values
            # step i
            ddi = dd[1].copy()

            for kk in vdd.keys():
                try:
                    k1, k2 = kk.split('_')
                except ValueError:
                    continue

                try:
                    ddi[k1][k2][0] = vdd[kk][i]
                except KeyError:
                    raise KeyError(
                        'Parameter %s of object %s not present in configuration dict.' % (
                            k2, k1
                        )
                    )

            # With modified input_dict, compute state
            Xi, labeldict = state_constructor(ddi)

            # Compute priors on new state
            prior, priorprob = compute_priors(priordict, labeldict)

            # Add prior to dictionary.
            vdd['prior'][i] = prior

        newvds.append(vdd)

    return newvds


def adapt_dict(indict, converters):
    """
    Tool that helps converting an MCMC output dictionary done with a previous
    version of PASTIS to the current one.

    Parameters
    ----------

    indict: dict or VDchain instances
        The dictionary to be adapted.

    converters: dict instance
        A dictionary containing the original names to be changes as keys and
        the output names as values.

    """
    if isinstance(indict, Chain) or isinstance(indict, VDchain):
        vd = indict.get_value_dict()

    elif isinstance(indict, dict):
        vd = indict

    else:
        print(
            'Warning! Instance not recognized. Assuming it is VDchain instance.')
        vd = indict.get_value_dict()

    vdout = {}
    for kk in vd.keys():

        kkk = kk.split('_')
        # If key is to be converted
        if kkk[0] in converters.keys():
            vdout[converters[kkk[0]] + '_' + kkk[1]] = vd.pop(kk)

        else:
            vdout[kk] = vd.pop(kk)

    return vdout


def get_model(vd, pastisfile, sampling=1, phase=False):
    """
    Compute the model for all samples of a given chain.

    Parameters
    ----------
    vd: dict or VDchain instance
        Chain for which to compute models

    t: array or list instance
        Times or phases (see keyword param phase) on which to compute the model.

    pastisfile: string or file instance
        Pastis file containing the configuration used to run the chain.

    
    Other parameters
    ----------------
    sampling: int
        Thinning factor (default = 1).

    phase: bool
        Decides if the given array t is a time or phase array (default = False)
    """

    # Check if vd is a dictionary or a VDchain instance
    if isinstance(vd, VDchain):
        vd = vd._value_dict

    elif isinstance(vd, dict):
        pass

    else:
        print('Type of input not recognized, assuming its a VDchain instance.')
        vd = vd._value_dict

    # Read pastis file and configuration dictionaries
    dd = pickle.load(open(pastisfile, 'r'))

    # Initialize


def compute_hdi_new(x, y, q=0.95):
    """
    Compute the 100*q% Highest Density Region, given samples from the posterior
    distribution and the posterior (up to a constant) evaluated in those
    positions.

    :param np.array x: posterior samples. This could be a multidimensional
    array as long as the samples run along the first axis.
    :param np.array y: posterior density evaluated in x.
    :param float q: fraction of mass contained in HDI.
    """

    # Sort posterior density array
    ind = np.argsort(y)

    # Keep only first q*n samples
    n = len(x)

    # How to report the result ? The question is on mode detection
    print(x)
    return x[ind[-int(n*q):]]


def compute_hdi(binedges, pdf, q=0.95):
    """
    Compute the 100*q% Highest Density Interval, given a normalised distribution
    pdf (len N), sampled in bins m (m has N+1 elements)
    """

    lower_bin_edges = binedges[: - 1]
    upper_bin_edges = binedges[1:]
    bincentre = 0.5 * (binedges[1:] + binedges[:-1])
    binsize = n.diff(binedges)

    # Sort elements from pdf
    isort = n.argsort(pdf)[::-1]

    cumulq = 0
    # Start adding bins until the requested fraction of the mass (q) is
    # reached
    for i, binnumber in enumerate(isort):
        cumulq += pdf[binnumber] * binsize[binnumber]

        if cumulq >= q:
            break

    # Keep only bins in 100*q% HDI
    bins_in_HDI = isort[: i + 1]

    # Sort binindex to find for non-contiguous bins
    sorted_binnumber = n.sort(bins_in_HDI)
    jumpind = n.argwhere(n.diff(sorted_binnumber) > 1)

    # Construct intervals
    HDI = []

    if len(jumpind) == 0:
        # HDI.append([binedges[sorted_binnumber[0]], binedges[sorted_binnumber[-1]+1]])
        HDI.append([lower_bin_edges[bins_in_HDI].min(),
                    upper_bin_edges[bins_in_HDI].max()])

    else:
        ji = 0
        for jf in jumpind[:, 0]:
            HDI.append([lower_bin_edges[sorted_binnumber[ji]],
                        upper_bin_edges[sorted_binnumber[jf]]])
            ji = jf + 1

        HDI.append([lower_bin_edges[sorted_binnumber[ji]],
                    upper_bin_edges[sorted_binnumber[-1]]])

    return HDI
