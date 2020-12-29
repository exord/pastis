import numpy as n
from scipy.interpolate import *
from scipy import optimize, stats
from math import *

from .constants import Msun, Rsun, G 

# Define useful functions employed by all packages
def area(x, y):
    """
    Compute area under curve y(x) over entire range of vectors.
    """
    return n.sum((x[1:] - x[:-1])*0.5*(y[1:] + y[:-1]))


def rebin(x, N):
    """
    Bin a given array x every n points.
    Takes mean of points in each bin.
    """
    if type(x) != type(n.array([])):
        x = n.array(x)
        
    Nx = len(x)
    nxb = int(n.ceil(float(Nx)/N))
    xb = n.zeros(nxb)
    for i, j in zip(range(0, Nx, N), n.arange(nxb)):
        xb[j] = (x[i: i+N]).mean()

    return xb

def rebin_texp(over_t, over_y, init_t, texp):
    """
    ReBin a given array over_y oversampled following the time vector over_t back to the initial time vector init_t with integration time texp.
    
    Inputs:
    -------
    
    over_t: array of N elements
    Oversampled time vector
    
    over_y: array of N elements
    Oversampled data vector
    
    init_t: array of M<N elements
    Binned time vector
    
    texp: array or float
    Binning timescale (in seconds)
    
    Output:
    -------
    
    init_y: array
    Binned data vector with init_t time reference
    """
    if len(over_t) != len(over_y):
        raise IOError('arguments over_t (%i) and over_y (%i) should have the same dimension' %(len(over_t), len(over_y)))
    if n.isscalar(texp):
        texp = n.zeros(len(init_t),float) + texp
    elif len(texp) != len(init_t):
        raise IOError('arguments texp (%i) and init_t (%i) should have the same dimension' %(len(texp), len(init_t)))
    
    init_y = n.zeros(len(init_t))
    for i in xrange(len(init_t)):
        init_y[i] = n.mean(over_y[n.where(n.logical_and(over_t >= init_t[i] - texp[i]/2./86400., over_t <= init_t[i] + texp[i]/2./86400.))])
    
    return init_y
    
    
def Z2z(Z):
    return log10(Z/0.019)


def z2Z(z):
    return 0.019*10**(z)


def trilinear_interpolation(input_array, indices):
    """
    Trilinear interpolation of input_array
    """
    
    output = n.empty(indices[0].shape)
    x_indices = indices[0]
    y_indices = indices[1]
    z_indices = indices[2]

    x0 = int(x_indices)
    y0 = int(y_indices)
    z0 = int(z_indices)
    x1 = x0 + 1
    y1 = y0 + 1
    z1 = z0 + 1

    #Check if xyz1 is beyond array boundary:
    if x1 >= input_array.shape[0]: x1 = input_array.shape[0] - 1
    if y1 >= input_array.shape[1]: y1 = input_array.shape[1] - 1
    if z1 >= input_array.shape[2]: z1 = input_array.shape[2] - 1

    x = x_indices - x0
    y = y_indices - y0
    z = z_indices - z0

    output = (input_array[x0,y0,z0]*(1-x)*(1-y)*(1-z) +
              input_array[x1,y0,z0]*x*(1-y)*(1-z) +
              input_array[x0,y1,z0]*(1-x)*y*(1-z) +
              input_array[x0,y0,z1]*(1-x)*(1-y)*z +
              input_array[x1,y0,z1]*x*(1-y)*z +
              input_array[x0,y1,z1]*(1-x)*y*z +
              input_array[x1,y1,z0]*x*y*(1-z) +
              input_array[x1,y1,z1]*x*y*z)

    return output

def bilinear_interpolation(input_array, indices):
    """
    bilinear interpolation of input_array
    """
    
    output = n.empty(indices[0].shape)
    x_indices = indices[0]
    y_indices = indices[1]

    x0 = int(x_indices)
    y0 = int(y_indices)
    x1 = x0 + 1
    y1 = y0 + 1

    #Check if xy1 is beyond array boundary:
    if x1 >= input_array.shape[0]: x1 = input_array.shape[0] - 1
    if y1 >= input_array.shape[1]: y1 = input_array.shape[1] - 1

    x = x_indices - x0
    y = y_indices - y0

    r1=(1-x)*input_array[x0,y0]+x*input_array[x1,y0]
    r2=(1-x)*input_array[x0,y1]+x*input_array[x1,y1]
    output=(1-y)*r1+y*r2

    return output


def eccanomaly(M, ecc, method='Newton', niterationmax=1e4, tol=1e-8):

    E = n.atleast_1d(M)
    Eo = E.copy()
    
    ecc = n.where(ecc > 0.99, 0.99, ecc)

    niteration = 0

    while n.any(n.abs(E - Eo) > tol/len(E)) or niteration==0:   
        
        Eo = E

        ff = E - ecc*n.sin(E) - M
        dff = 1 - ecc*n.cos(E)

        if method == 'Newton':
            # Use Newton method
            E = Eo - ff / dff

        elif method == 'Halley':
            # Use Halley's parabolic method
            d2ff = ecc*n.sin(E)
        
            discr = dff **2 - 2 * ff * d2ff

            E = n.where((discr < 0), Eo - dff / d2ff,
                         Eo - 2*ff / (dff + n.sign(dff) * n.sqrt(discr))
                     )

        # Increase iteration number; if above limit, break with exception.
        niteration += 1
        if niteration >= niterationmax:
            raise RuntimeError('Eccentric anomaly computation not converged.')
        
    return E


def trueanomaly(ea, ecc):
    """
    Compute true anomaly for eccentricity and eccentric anomaly (ea)

    :param np.array ea: eccentric anomaly in radians
    :param np.array ecc: eccentricity
    """
    nu = 2. * n.arctan2(n.sqrt(1. + ecc) * n.sin(ea / 2.),
                        n.sqrt(1. - ecc) * n.cos(ea / 2.)
                        )

    return nu
    
    

def sigma_clipping(y, k, ey = None, niter = 1, useMAD = False):
    """
    Sigma-clips array y at the k-sigma level.
    Returns clipped array and a boolean array with the condition to be
    used to clip other arrays, if necessary.

    Other parameters
    ----------------
    ey: array.
        Vector with errors for array y, used to weight the sigma computation.
        If None, equal weigth is given to all points.

    niter: int.
        Number of times the clipping is iterated.
        
    useMAD: boolean.
        Decides if the Mean Absolute Deviation is used to compute the
        standard deviation, assuming the distribution is Gaussian.
    """

    if ey != None:
        w = 1/ey**2
    else:
        w = n.zeros(len(y)) + 1.0

    # Array that will contain condition
    cond = n.zeros(len(y), bool)

    # Masked arrays
    y = n.ma.array(y, mask = cond)
    w = n.ma.array(w, mask = cond)
    
    for i in range(niter):
        
        # Compute mean value
        mean_y = n.sum(w*y)/n.sum(w)

        if useMAD:
            sigma = n.median( n.abs(y - mean_y) )*1.4826

        else:
            sigma = n.sqrt( n.average( (y - mean_y)**2.0, weights = w ) )

        ## Condition to mask arrays (i.e. those elements above or below k*sigma
        condi = n.greater(n.abs(y - mean_y), k*sigma)

        y.mask = condi
        w.mask = condi

    return y.compressed(), -y.mask
        
def density_prior(teff, sigma_teff, logg, sigma_logg, z, sigma_z, N = 1e3,
                  verbose = False, EM = 'Dartmouth'):
    """
    Compute the expected stellar density given the measurements of the
    atmospheric parameters: Teff, logg, and z.
    Assumes that the measurements are uncorrelated and Gaussian distributed.

    Parameters
    ---------
    teff: float.
        Mean value of Teff

    sigma_teff: float.
        Width of gaussian distribution for Teff.
        
    logg: float.
        Mean value of logg

    sigma_logg: float.
        Width of gaussian distribution for logg.
        
    z: float.
        Mean value of [Fe/H]

    sigma_z: float.
        Width of gaussian distribution for [Fe/H].
    
    N: int.
        Number of samples to draw from distributions (default: 1e3).
        
    verbose: bool
        Print errors
    
    EM: str
        Evolution model: Dartmouth, Geneva, StarEvol, ...

    Returns
    -------
    dens: array.
        Array with the density values corresponding to the drawn values of
        the stellar atmospheric parameters. Its dimension might differ from N,
        since some combinations of Teff, logg and z might be unrealistic.
    """

    from . import AstroClasses as ac

    ## Check if initialization is needed
    try:
        from .isochrones import maxz
    except ImportError:
        from .isochrones import prepare_tracks_target
        from . import EMdict
        prepare_tracks_target(EMdict[EM])
        reload(ac)
        
    Teff = stats.norm.rvs(teff, sigma_teff, size = N)
    logg = stats.norm.rvs(logg, sigma_logg, size = N)
    z = stats.norm.rvs(z, sigma_z, size = N)

    density = n.zeros(N, float)

    for i in xrange(int(N)):
        if i%50 == 0.0 and verbose:
            print(i)

        try:
            s1 = ac.Target(z = z[i], teff = Teff[i], logg = logg[i],
                           dist = 1.0, ebmv = 0.0)
        except Exception as ee:
            if verbose: print(ee.message)
            continue
        else:
            density[i] = s1.mact/s1.R**3

    return n.array(density)


def fithist(y, p0, nbins = 50, disttype = 'Normal', plot = True, **kwargs):
    """
    Fit histogram to a probability distribution of type disttype

    Parameters
    ----------

    y: array
        Array of data to fit

    p0: iterable
        An iterable containing the initial paramters for the fit

    Other parameters
    ----------------
    nbins: int
        Number of bins used to construct the histogram

    disttype: string
        Distribution to fit. Must be a valid prior type:

        Uniform

        Jeffreys

        Normal
        
        Binormal

        AsymmetricNormal

        TruncatedUNormal

        PowerLaw

        Sine

        ...

    plot: boolean
        If True, the best fit model is plotted.
        
    Additional keyword arguments are passed to scipy.optimize.leastsq function.
    """
    import MCMC.priors as priors

    distdict = priors.distdict

    ## Construct histogram
    N, binedges = n.histogram(y, bins = nbins, normed = True)
    x = binedges[:-1] + 0.5*n.diff(binedges)

    x = x[N != 0]
    N = N[N != 0]
    
    # Compute error in bin numbers
    eN = n.sqrt(N)/sqrt(len(y)*(x[1] - x[0]))
    
    def residuals(p, x, y, ey):
        return (y - distdict[disttype][0].pdf(x, *p))/ey

    popt, epopt, infodict, mesg, wf = optimize.leastsq(residuals,
                                                       p0, args = (x, N, eN),
                                                       full_output = True,
                                                       **kwargs)

    chi2 = n.sum(infodict['fvec']**2)
    ndof = nbins - 1
    print('Chi2: %.2f' % chi2)
    print('Number of degrees of freedom: %d' % ndof)
    print('p-value: %.2e' % stats.chi2.sf(chi2, ndof))

    if plot:
        import pylab as p
        f1 = p.figure()
        ax = f1.add_subplot(111)

        ax.hist(y, nbins, histtype = 'step', lw = 2, color = 'k', normed = True)
        ax.errorbar(x, N, eN, fmt = 'ko')
        
        xx = n.linspace(x.min(), x.max(), 1000)
        ax.plot(xx, distdict[disttype][0].pdf(xx, *popt), lw = 2, color = 'r')



        p.draw()
                
    return popt

    
def exportresults(vd, outfile):
    """
    Writes the contents of dictioanry vd to a tab-separated file, outfile.
    
    Parameters
    ----------
    vd, dict.
    
    outfile, string or file
    
    """
    
    if isinstance(outfile, str):
        fout = open(outfile, 'w')
    else:
        fout = outfile

    ## Get keys and values from dictionary
    values = n.array(vd.values())
    keys = vd.keys()

    ## Write header
    headerstr = ''
    for key in keys:
        headerstr = headerstr+'\t'+key
    
    ## Write header to file
    fout.write(headerstr+'\n')

    ## Write values of dictionary
    fmtstr = '%.12f\t'*len(keys)+'\n'

    for line in values.T:
        fout.write(fmtstr%tuple(line))
    
    fout.close()
    
    return


def leastsq(pastisfile, vd0 = {}, initialize = True, plot = False,
            **kwargs):
    """
    Performs a least square fit using a configuration given by the
    pastisfile parameter and initial values included in dictionary vd0. 

    Parameters
    ----------

    pastisfile: str, file, or list
        The configuration file can be a string with the full path of the
        pastis file, or a file instance. For additional flexibility,
        alternatively the configuration dictionaries can be given in a list.
        In this case, the list must contain: info_dict, input_dict, datadict,
        custompriordict
        
    vd0: dictionary containing the starting parameters
    """
    

    import pickle
    
    ### Construct residual function that will be used for least square
    def resLM(p, keys, input_dict, datadict):

        # Construct initial state from info in input_dict
        X, labeldict = PASTIS_MCMC.state_constructor(input_dict)

        ## Put update values in labeldict
        for i, k in enumerate(keys):
            labeldict[k].set_value(p[i])

        try:
            res = PASTIS_MCMC.get_likelihood(X, input_dict, datadict, labeldict,
                                             False, False, LM = True)
        except EBOPparamError as err:
            print(err.message)
            res = 1e100

        return res
    

    if isinstance(pastisfile, str) or isinstance(pastisfile, file):
        if isinstance(pastisfile, str):
            ## Read configuration file
            f = open(pastisfile, 'r')
        else:
            f = pastisfile
            
        dd = pickle.load(f)
        f.close()

        info_dict = dd[0].copy()
        input_dict = dd[1].copy()
        datadict, sp, lc = DataReader(dd[2])
        custompriordict = dd[3]

    elif iter(pastisfile):
        info_dict = pastisfile[0].copy()
        input_dict = pastisfile[1].copy()
        datadict = pastisfile[2].copy()
        custompriordict = pastisfile[3].copy()
        

    if initialize:
        from . import initialize
        initialize(info_dict, datadict, input_dict)
    
    from MCMC import PASTIS_MCMC

    # Construct initial state from info in input_dict
    X, labeldict = PASTIS_MCMC.state_constructor(input_dict)

    # Prepare dictionary with prior functions
    priordict = PASTIS_MCMC.prior_constructor(input_dict, custompriordict)

    # Start fit at random point in prior
    PASTIS_MCMC.pick_random_point(labeldict, priordict)

    # Replace random starting points by user-defined points
    for key in vd0.keys():
        labeldict[key].set_value(vd0[key])

    # Create parameter array
    p0 = []
    k = []
    for xx in X:
        if xx.jump:
            p0.append(xx.get_value())
            k.append(xx.label)
            

    ## Iterate a few times the fitting procedure
    popt, cov, infodict, mesg, ier = \
        optimize.leastsq(resLM, p0, args = (k, input_dict, datadict),
                         full_output = True, **kwargs
                         )
    
    if ier in [1, 2, 3, 4]:
        ## Solution was found!
        print('Solution succesfully reached!')
        for i in range(len(k)):
            labeldict[k[i]].set_value(popt[i])
            print('%s: %f'%(k[i], popt[i]))
        print('====')
        print('Chi2 = %.4f for %d points.'%(n.sum(infodict['fvec']**2),
                                            len(infodict['fvec'])
                                            )
              )

        if plot == True:
            ## Create solution file
            for kk in input_dict.keys():
                for kkk in input_dict[kk].keys():
                    input_dict[kk][kkk][0] = labeldict[kk+'_'+kkk].get_value()
            f = open('temp.dat', 'w')
            pickle.dump((info_dict, input_dict, datadict, custompriordict), f)
            f.close()

            pp.plot_dict(f.name, mergeRV = True)
        
    else:
        print('Not solution found!')
        
    return popt, k, cov, infodict, mesg, ier
        
def get_E0(ecc, omega):
    """
    Compute the eccentric anomaly at the time of transit.
    """
    nu0 = pi/2.0 - omega
    return 2.0*n.arctan2(n.sin(nu0/2.0) * n.sqrt(1 - ecc),
                         n.cos(nu0/2.0) * n.sqrt(1 + ecc))
    #return n.arctan2(n.sqrt(1. - ecc**2)*n.cos(omega), n.sin(omega) + ecc)

def get_Ea(ecc, omega):
    """
    Compute the eccentric anomaly at the time of secondary transit.
    """
    nu_a = -(pi/2.0 + omega)
    return 2.0*n.arctan2(n.tan(nu_a/2.0)*n.sqrt(1 - ecc), n.sqrt(1 + ecc))
    
def get_T0(ecc, omega, P, Tp): 
    """
    Return the transit epoch from other orbital parameters
    """
    E_0 = get_E0(ecc, omega)
    T0 = Tp + P/(2.*n.pi)*(E_0 - ecc*n.sin(E_0))
    return T0

def get_Tp(ecc, omega, P, T0):
    """ 
    Return the epoch of periastron from other orbital parameters
    """
    E_0 = get_E0(ecc, omega)
    Tp = T0 - P/(2.*n.pi)*(E_0 - ecc * n.sin(E_0))
    return Tp


def get_T02(ecc, omega, P, Tp): 
    """
    Return the secondary transit epoch from other orbital parameters
    """
    E_a = get_Ea(ecc, omega)
    Ta = Tp + P/(2.*n.pi)*(E_a - ecc*n.sin(E_a))
    
    return Ta

def derivative(x, y):
    """
    Compute derivative using slope
    """
    d = n.zeros(len(x))
    d[:-1] = (y[1:] - y[:-1])/(x[1:] - x[:-1])
    d[-1] = d[-2]
    
    #substitute infinites by zeros
    d[n.where(d  == n.inf)] = 0

    return d

    
def iterative_mass(K, P, Ms, ecc, y, omega=None, Rs=None, use_b=False,
                   **kwargs):
    """
    Obtain mass of secondary iteratively, without recourse to q = m2/m1,
    or assuming any value of q.

    The trascendental equation is solved with the Newton Raphson method.

    The output is in solar masses.
    
    Parameters
    ----------
    K: semi-amplitude of RV variation (in km/s)

    P: orbital period of the system (in days)

    Ms: mass of the primary star (in solar masses)

    ecc: orbital eccentricity

    y: orbital inclination (in degrees) or impact parameter

    omega: argument of the pericenter; only needed if use_b = True.

    Rs: radius of the primary star (in solar radii); only needed if use_b = True
    
    use_b : boolean; decides if y represents the impact parameter of the inclination

    The remaining parameters are passed to scipy.optimize.newton
    """
    import scipy.optimize

    K_ms = K*1e3
    P_s = P*24*3600
    Ms = Ms * Msun

    if not use_b:
        incl = y*pi/180.0

        # Define function for Newton-Raphson method.
        def ff(m2, K_ms, P_s, Ms, ecc, incl):
            return K_ms*(P_s * (Ms + m2)**2 / (2*pi*G))**(1./3.) * \
                n.sqrt(1 - ecc**2) / n.sin(incl) - m2

        # Define derivative function.
        def dff(m2, K_ms, P_s, Ms, ecc, incl):
            return 2.0/3.0 * K_ms*(P_s * (Ms + m2)**(-1.0) / (2*pi*G))**(1./3.) * n.sqrt(1 - ecc**2) / n.sin(incl) - 1

        # Define starting points (mass with q = 0.0)
        x0 = n.mean(K_ms*(P_s * Ms**2 / (2*pi*G))**(1./3.) *
                    n.sqrt(1 - ecc**2) / n.sin(incl))

        mact = scipy.optimize.newton(ff, x0, fprime=dff,
                                     args=(K_ms, P_s, Ms, ecc, incl),
                                     **kwargs)

    else:
        #####
        # If b is given instead of incl, things are more complicated.
        #####
        
        # First, check if omega and Rs are defined; if so, convert to radians
        # and meters
        if omega is not None and Rs is not None:
            omega = omega * pi/180.0
            Rs = Rs * Rsun
        else:
            raise TypeError('omega and Rs are needed if impact parameter is'
                            ' given as input.')

        # Compute true anomaly at transit
        nu0 = pi/2.0 - omega

        # Compute cos(incl)*(Ms + m2)**(1/3)
        xx = (4*pi**2/G)**(1./3.) * y * Rs * P_s**(-2./3.) * \
             (1 + ecc*cos(nu0)) / (1 - ecc**2.0)

        def auxf(m2, K_ms, P_s, Ms, ecc, xx):
            return K_ms *P_s**(1./3.) * sqrt( 1 - ecc**2) / \
                ( (2*pi*G)**(1./3.) * sqrt( (Ms + m2)**(2./3.) - xx**2) )

        def ff(m2, K_ms, P_s, Ms, ecc, xx):
            return auxf(m2, K_ms, P_s, Ms, ecc, xx) * (Ms + m2) - m2
            
        def dff(m2, K_ms, P_s, Ms, ecc, xx):
            return auxf(m2, K_ms, P_s, Ms, ecc, xx) * \
                (1 - (Ms + m2)**(2./3.) / (3.0 * ((Ms + m2)**(2./3.) - xx**2)
                                        )
                ) - 1.0

        # Define starting points (mass with q = 0.0)
        x0 = K_ms * P_s**(1./3.) * sqrt( 1 - ecc**2) * Ms /\
             ( (2*pi*G)**(1./3.) * sqrt( Ms**(2./3.) - xx**2) )
        
        mact = scipy.optimize.newton(ff, x0, fprime = dff, 
                                     args = (K_ms, P_s, Ms, ecc, xx),
                                     **kwargs
                                 )

    if mact < 0.0:
        return 0.0
    else:
        return mact/Msun
    

def fsincos(x, *p):
    '''
    Sine function with free phase for fitoc function
    '''
    A, B, shift = p
    return A*n.sin(2.*n.pi*x)+B*n.cos(2.*n.pi*x)+shift

def fitoc(datafile, merged, nperiods = 100., pmin = 0.1, pmax = None, outfile = 'simres.pickle', logscale = True):
    '''
    Fit a sine function to the residuals of an RV solution

    datafile: str
        RV datafile

    merged: str
        merged derived chain (dict), should contain: P, T0, K1, v0, drift1, ecc, omega,  mact

    nperiods: int
        number of computed periods

    pmin: float
        minimum period

    pmax: float
        maximum period

    outfile: str
        name of the output pickle file with the results
        the results are the periods and the computed planet mass (in Jupiter masses) for each point of the chain
        To be ploted with  plotfit

    logscale: bool
        if True the periods are logarithmic sampled, if False linear

    '''
    import AstroClasses as ac
    import sys
    import pickle
    from scipy.optimize import curve_fit

    # read merged chain 
    f = open(merged,'r')
    vd = pickle.load(f)
    f.close()

    # read RV data 
    time, rv, srv = loadtxt_iter(datafile, usecols = (0,1,2), unpack=True, skiprows=2)

    # Period array
    if pmax == None: pmax = 0.5*(n.max(time)-n.min(time))
    if logscale:
        periods = n.logspace(n.log10(pmin), n.log10(pmax), nperiods)
    else:
        periods = n.linspace(pmin, pmax, nperiods)

    nchain = len(vd['FitPlanet1_P'])
    nperiods = len(periods)
    data = n.zeros((nchain, nperiods))
    Ntotal = nchain * nperiods

    # Find keys
    for key in vd.keys():
        if '_P' in key: vdp = vd[key]
        if '_T0' in key: vdt0 = vd[key]
        if '_K1' in key: vdK = vd[key]
        if '_v0' in key: vdv0 = vd[key]
        if '_drift1' in key: vddrift = vd[key]
        if '_ecc' in key: vde = vd[key]
        if '_omega' in key: vdomega = vd[key]
        if '_mact' in key: vdmass = vd[key]
        
    try:
        nothing = vddrift[0]
    except:
        vddrift = n.zeros(len(vdp))

    for i in range(nchain):
        # calulo de residuos
        dd = {'P': vdp[i], 'T0': vdt0[i], 'K1': vdK[i], 'v0': vdv0[i], 'drift1': vddrift[i], 'ecc': vde[i], 'omega': vdomega[i]}
        planet = ac.FitPlanet(**dd)
        oc = rv - planet.get_RV(time)

        for j in range(nperiods):

            ni = i*nperiods + j + 1
            nip = 1e2*float(ni)/Ntotal

            if j==0 and i==0:
                sys.stdout.write('Computing... %02d %%'%nip)
            elif n.round(nip)%1 == 0:
                sys.stdout.write('\b'*4+'%02d %%'%nip)
            sys.stdout.flush()

            phase = (time/periods[j]) % 1.0
                
            # fit phase, oc
            coeff, var_matrix = curve_fit(fsincos, phase, oc, p0 = [0.1, 0.1, 0.1],sigma = srv)
            # error de K = n.sqrt(var_matrix[0][0])
            K = n.sqrt(coeff[0]**2+coeff[1]**2)
            data[i,j] = (K*1e3*(vdmass[i]**(2./3.))*(periods[j]**(1./3.)))/203.0 # Mjup 

    f = open(outfile,'w')
    pickle.dump([periods, data], f)
    f.close()



def fitoc(datafile, merged, nperiods = 100., pmin = 0.1, pmax = None, outfile = 'simres.pickle', logscale = True):
    '''
    Fit a sine function to the residuals of an RV solution

    datafile: str
        RV datafile

    merged: str
        merged derived chain (dict), should contain: P, T0, K1, v0, drift1, ecc, omega,  mact

    nperiods: int
        number of computed periods

    pmin: float
        minimum period

    pmax: float
        maximum period

    outfile: str
        name of the output pickle file with the results
        the results are the periods and the computed planet mass (in Jupiter masses) for each point of the chain
        To be ploted with  plotfit

    logscale: bool
        if True the periods are logarithmic sampled, if False linear

    '''
    import AstroClasses as ac
    import sys
    import pickle
    from scipy.optimize import curve_fit

    # read merged chain 
    f = open(merged,'r')
    vd = pickle.load(f)
    f.close()

    # read RV data 
    time, rv, srv = loadtxt_iter(datafile, usecols = (0,1,2), unpack=True, skiprows=2)

    # Period array
    if pmax == None: pmax = 0.5*(n.max(time)-n.min(time))
    if logscale:
        periods = n.logspace(n.log10(pmin), n.log10(pmax), nperiods)
    else:
        periods = n.linspace(pmin, pamx, nperiods)

    nchain = len(vd['FitPlanet1_P'])
    nperiods = len(periods)
    data = n.zeros((nchain, nperiods))
    Ntotal = nchain * nperiods

    # Find keys
    for key in vd.keys():
        if '_P' in key: vdp = vd[key]
        if '_T0' in key: vdt0 = vd[key]
        if '_K1' in key: vdK = vd[key]
        if '_v0' in key: vdv0 = vd[key]
        if '_drift1' in key: vddrift = vd[key]
        if '_ecc' in key: vde = vd[key]
        if '_omega' in key: vdomega = vd[key]
        if '_mact' in key: vdmass = vd[key]
        
    try:
        nothing = vddrift[0]
    except:
        vddrift = n.zeros(len(vdp))

    for i in range(nchain):
        # calulo de residuos
        dd = {'P': vdp[i], 'T0': vdt0[i], 'K1': vdK[i], 'v0': vdv0[i], 'drift1': vddrift[i], 'ecc': vde[i], 'omega': vdomega[i]}
        planet = ac.FitPlanet(**dd)
        oc = rv - planet.get_RV(time)

        for j in range(nperiods):

            ni = i*nperiods + j + 1
            nip = 1e2*float(ni)/Ntotal

            if j==0 and i==0:
                sys.stdout.write('Computing... %02d %%'%nip)
            elif n.round(nip)%1 == 0:
                sys.stdout.write('\b'*4+'%02d %%'%nip)
            sys.stdout.flush()

            phase = (time/periods[j]) % 1.0
                
            # fit phase, oc
            coeff, var_matrix = curve_fit(fsincos, phase, oc, p0 = [0.1, 0.1, 0.1],sigma = srv)
            # error de K = n.sqrt(var_matrix[0][0])
            K = n.sqrt(coeff[0]**2+coeff[1]**2)
            data[i,j] = (K*1e3*(vdmass[i]**(2./3.))*(periods[j]**(1./3.)))/203.0 # Mjup 

    f = open(outfile,'w')
    pickle.dump([periods, data], f)
    f.close()



    
def plotfitoc(pdata = 'simres.pickle', log = False):
    '''
    plot fitoc results

    pdata: str
        pickle file name with the results of fitoc

    log: bool
        if True the x-axis in the plot is in log scale

    '''
    import pickle
    import pylab as p
    from matplotlib import rc
    rc('text', usetex=True)
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    
    f = open(pdata,'r')
    periods, data = pickle.load(f)
    f.close()
    
    nperiods = len(periods)

    Mass = n.zeros((nperiods,3))
    
    for i in range(nperiods):
        perioddata = data[:,i]
        Mass[i,0] = n.percentile(perioddata, 68.27)
        Mass[i,1] = n.percentile(perioddata, 95.45)
        Mass[i,2] = n.percentile(perioddata, 99.73)

        """
        y, x, z = p.hist(perioddata, nbins, normed = True, cumulative = True)
        xm = 0.5*(x[1:]+x[:-1])

        Mass[i,0] = n.interp(0.6827,y,xm)
        Mass[i,1] = n.interp(0.9545,y,xm)
        Mass[i,2] = n.interp(0.9973,y,xm)
        """
   
    f1 = p.figure()
    ax = f1.add_subplot(111)
    ax.fill_between(periods, 0, Mass[:,2],color='r', facecolor='r', alpha=0.1)
    ax.fill_between(periods, 0, Mass[:,1],color='r', facecolor='r', alpha=0.2)
    ax.fill_between(periods, 0, Mass[:,0],color='r', facecolor='r', alpha=0.5)
    lpad = 10
    fontsize = 16
    if log: ax.set_xscale('log')
    ax.set_xlabel('Period [days]', fontsize = fontsize,labelpad=lpad)
    ax.set_ylabel('Upper limit 2nd planet mass [M$_J$]', fontsize = fontsize,labelpad=lpad)
    ax.set_xlim(min(periods),max(periods))

    for xtl in ax.get_xticklabels():
        xtl.set_fontsize(fontsize)

    for xtl in ax.get_yticklabels():
        xtl.set_fontsize(fontsize)

    p.draw()


def fitoc2(datafile, merged, nperiods = 100., nsamples = 200, pmin = 0.1, pmax = None, outfile = 'simres.pickle', logscale = True):
    '''
    Fit a sine function to the residuals of an RV solution

    datafile: str
        RV datafile

    merged: str
        merged derived chain (dict), should contain: P, T0, K1, v0, drift1, ecc, omega,  mact

    nperiods: int
        number of computed periods

    pmin: float
        minimum period

    pmax: float
        maximum period

    outfile: str
        name of the output pickle file with the results
        the results are the periods and the computed planet mass (in Jupiter masses) for each point of the chain
        To be ploted with  plotfit

    logscale: bool
        if True the periods are logarithmic sampled, if False linear

    '''
    import AstroClasses as ac
    import sys
    import pickle
    from scipy.optimize import curve_fit

    # read merged chain 
    f = open(merged,'r')
    vd = pickle.load(f)
    f.close()

    ### Define function phase
    def o(x):
        return x*0.0 + 1.0

    def s(x):
        return n.sin(2.*n.pi*x)
        
    def c(x):
        return n.cos(2.*n.pi*x)

    
    # read RV data 
    time, rv, srv = loadtxt_iter(datafile, usecols = (0,1,2), unpack=True, skiprows=2)

    # Period array
    if pmax == None: pmax = 0.5*(n.max(time)-n.min(time))
    if logscale:
        periods = n.logspace(n.log10(pmin), n.log10(pmax), nperiods)
    else:
        periods = n.linspace(pmin, pamx, nperiods)

    nchain = len(vd['FitPlanet1_P'])
    nperiods = len(periods)
    data = n.zeros((nsamples, nperiods))
    Ntotal = nsamples * nperiods

    # Find keys
    for key in vd.keys():
        if '_P' in key: vdp = vd[key]
        if '_T0' in key: vdt0 = vd[key]
        if '_K1' in key: vdK = vd[key]
        if '_v0' in key: vdv0 = vd[key]
        if '_drift1' in key: vddrift = vd[key]
        if '_ecc' in key: vde = vd[key]
        if '_omega' in key: vdomega = vd[key]
        if '_mact' in key: vdmass = vd[key]
        if key == 'SOPHIE HE-K5_jitter': jitter = vd[key]
        
    try:
        nothing = vddrift[0]
    except:
        vddrift = n.zeros(len(vdp))

    OC = n.zeros((nchain, len(time)))

    ##
    import random
    indChain  = range(nchain)
    random.shuffle(indChain)
    
    for ii, i in enumerate(indChain[:nsamples]):#range(nchain[:100]):
        # calulo de residuos
        dd = {'P': vdp[i], 'T0': vdt0[i], 'K1': vdK[i], 'v0': vdv0[i], 'drift1': vddrift[i], 'ecc': vde[i], 'omega': vdomega[i]}
        planet = ac.FitPlanet(**dd)
        oc = rv - planet.get_RV(time)

        ## Correct for jitter
        erv = n.sqrt(srv**2 + jitter[i]**2)
                     
        for j in range(nperiods):

            ni = ii*nperiods + j + 1
            nip = 1e2*float(ni)/Ntotal

            if j==0 and ii==0:
                sys.stdout.write('Computing... %02d %%'%nip)
            elif n.round(nip)%1 == 0:
                sys.stdout.write('\b'*4+'%02d %%'%nip)
            sys.stdout.flush()

            phase = (time/periods[j]) % 1.0
                
            # fit phase, oc
            coeff, ecoff, cov = leastsqlin(phase, oc, erv, (s, c, o))
            # error de K = n.sqrt(var_matrix[0][0])
            K = n.sqrt(coeff[0]**2+coeff[1]**2)
            data[ii,j] = (K*1e3*(vdmass[i]**(2./3.))*(periods[j]**(1./3.)))/203.0 # Mjup 

    f = open(outfile,'w')
    pickle.dump([periods, data], f)
    f.close()

def leastsqlin(x, y, ey, f, *args, **kwargs):
    """
    Solves a least-squares linear problem.

    The model function is a*f, with f an iterable instance of functions
    """
    from numpy import linalg

    ## Check if f is a single function
    if hasattr(f, '__call__'):
        At = n.array(f(x, *args, **kwargs))
    elif hasattr(f, '__iter__'):
        At = n.zeros([len(f), len(x)])
        for i, func in enumerate(f):
            At[i] = func(x, *args, **kwargs)

    A = At.T
    V = n.diag(1.0/ey**2.0)

    MA = n.dot(At, n.dot(V, A))
    MB = n.dot(At, n.dot(V, y.reshape(len(y), 1)))

    # Need to invert array MA
    MAi = linalg.inv(MA)

    # Compute parameters and their covariance matrix
    par = n.dot(MAi, MB)
    vpar = MAi
    
    return n.array(par)[:,0],n.sqrt(n.diag(vpar)),vpar

    
def loadtxt_iter(filename, delimiter=None, comments = '#',
                 skiprows=0, dtype = float,
                 usecols = None, unpack = False):
    """
    Memory-light version of numpy.loadtxt
    """
    ## Define function that returns generator 
    def iter_func(usecols = usecols):
        with open(filename, 'r') as infile:
            for _ in range(skiprows):
                next(infile)
            for line in infile:
                ## Strip \n and remove comments
                line = line.split(comments)[0].rstrip()
                if line:
                    ## Split at delimiter
                    line = n.array(line.split(delimiter))
                else:
                    continue
                    
                if usecols == None:
                    usecols = range(len(line))
                for item in line[list(usecols)]:
                    yield dtype(item)
        loadtxt_iter.rowlength = len(usecols)

    # Create data from this generator
    data = n.fromiter(iter_func(usecols = usecols), dtype=dtype)
    
    # Reshape data
    data = data.reshape((-1, loadtxt_iter.rowlength))

    if unpack:
        return data.T
    else:
        return data
