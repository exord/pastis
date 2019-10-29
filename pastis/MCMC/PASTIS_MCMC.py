# PASTIS_MCMC
# 23-5-2012 Input dictionaries (input_dict, datadict and customprior_dict are
#           now copied to Chain attribute _inputdicts). Need to include infodict
#  2-4-2012 Due to errors, I comment: Every 10,000 step, save pickle Chain
# 26-3-2012 bug in MCMC (BoundingError instead of objMCMC.BoundingError)
# 23-3-2012 solved bugs with out of bounds and PASTIS_RV call,
#           del  globals()['global_spectrum']
# v1.0 2012-03-18: Classes Parameter, Prior and Chain are now defined
# in module objMCMC, in order to be able to pickle them.

from math import pi, log10, log, exp
import random
import time
import pickle
from scipy import stats
import scipy
from scipy.stats.distributions import rv_frozen
import numpy as n
import traceback

# Intra-package imports
from .tools import state_constructor, state_deconstructor, get_jitter
from .priors import compute_priors, prior_constructor
from . import Objects_MCMC as objMCMC
from. import PCA

# Imports from upper package level
# from .. import inputdicts
from .. import ObjectBuilder
from .. import AstroClasses as ac
from ..tools import rebin, rebin_texp
from ..photometry import corot_colors
from ..models import SED
from ..models import RVgaussianFitError
from ..models.PHOT import PASTIS_PHOT
from ..models.SED import PASTIS_SED
from ..models.RV import PASTIS_RV
from .. import checkblendmag
#
from ..exceptions import EvolTrackError, OutofIsochroneError, EBOPparamError


def mcmc(input_dict, datadict, customprior_dict, N, chain=None,
         beta=1.0, Npca=n.inf, NupdatePCA=n.inf, kPCA=1.0, Nlastupdate=n.inf,
         BIpca=0, randomstart=True, usePCA=True, mindiff=1e-4, AM=False,
         autocorrect=False, **kwargs):

    """
    Main MCMC function

    Parameters
    ----------
    input_dict: dict.
        Input dictionary containing the information about the parameters
        that will be explored by the MCMC algorithm and their priors.

    datadict: dict.
        Dictionary produced by DataReader containing information about the
        observations used for the computation of the likelihood.

    customprior_dict: dict.
        Dictionary containing information about priors concerning more than one
        parameter.

    N: int.
        Number of steps of the chain.

    beta: float.
        Tempering parameter used to modified the target distribution.

    Npca: int.
        Iteration at which the principal component analysis (PCA) used
        to reduce the correlations of parameters is started.

    NupdatePCA: int.
        Number of steps taken between updated of the covariance matrix used
        for the PCA.

    kPCA: float.
        When updating the covariance matrix every NupdatePCA steps, the last
        kPCA*NupdatePCA iterations used for the computation.

    Nlastupdate: int.
        Iteration number in which updating of the covariance matrix stops.
        
    BIpca: int.
        Number of iterations to reject for the first computation of the
        covariance matrix.

    randomstart: bool.
        Determines whether the chain is started at a random point in paramter
        space drawn from the joint prior.

    usePCA: bool.
        Determines if PCA is used. Added for compatibility with widget.

    mindiff: float.
        Minimum relative difference of covariance and principal components
        analysis for which an update is performed. DESACTIVATED!

    AM: bool.
        Determines if an adaptive step size is used for each parameter or
        if jumps are taken from a multinormal distribution with adaptive
        covariance matrix, as in Haario (2001).

    outputfile: str.
        The name of the outputfile where preliminary steps of the chain will
        be stocked.

    Nsave: int.
        Number of steps between preliminary dumping of results.
    
       
    """

    outputfile = kwargs.pop('outputfile', None)
    Nsave = kwargs.pop('Nsave', 1e4)
    
    # INITIALISATION
    
    # Construct initial state from info in input_dict
    X, labeldict = state_constructor(input_dict)

    # Prepare dictionary with prior functions
    priordict = prior_constructor(input_dict, customprior_dict)

    # If requested, start chain at random point in prior
    if randomstart:
        pick_random_point(labeldict, priordict)

    # Compute prior for initial state
    priorx, priorprobx = compute_priors(priordict, labeldict)

    # Compute likelihood for initial state
    haspassed = False
    while not haspassed:
        try:
            Lx, logLx, likedictx = get_likelihood(X, input_dict, datadict,
                                                  labeldict, False, autocorrect)

        except RVgaussianFitError as rvfite:
            print(rvfite)
            print('Input parameters: {0}, {1}, {2}'.format(rvfite.contrast,
                                                           rvfite.rv0,
                                                           rvfite.sigma))
            pick_random_point(labeldict, priordict)
                
        except (EvolTrackError, OutofIsochroneError):
            print('Something went wrong. Initial objects cannot be constructed '
                  'due to limitations in the stellar evolution tracks.')
            pick_random_point(labeldict, priordict)

        except RuntimeError:
            print('Something went wrong. Initial objects cannot be constructed '
                  'due to a Runtime error (possibly Qhull!).')
            pick_random_point(labeldict, priordict)

        except EBOPparamError as eboperr:
            print('Parameter outside limits for JKTEBOP; '
                  'error message: {}'.format(eboperr))
            pick_random_point(labeldict, priordict)

        except ValueError as verr:
            print('Value Error. Trying new starting point. '
                  'Message: {}'.format(verr))
            traceback.print_exc()
            pick_random_point(labeldict, priordict)

        else:
            haspassed = True

    for xx in X:
        print(xx.label, xx.get_value())

    # Initialise chain if no chain is given
    if chain is None:
        markov_chain = objMCMC.Chain(N, X, Npca, usePCA)
        markov_chain.TrackError = []
        # markov_chain._inputdicts = [inputdicts[0], input_dict, datadict,
        #                 customprior_dict]
        markov_chain._inputdict = [input_dict, datadict, customprior_dict]
    
    elif isinstance(chain, objMCMC.Chain):
        markov_chain = chain
        X = markov_chain.get_current_state()
    else:
        raise TypeError('chain must be None or Chain object')

    # Include posterior and likelihood in chain
    markov_chain._posterior[0] = logLx / log(10) + log10(priorx)
    markov_chain._likelihood[0] = [Lx, logLx, likedictx]
    
    ####
    
    # Iterate desired number of steps
    ti = time.time()

    N = int(N)
    for i in range(1, N):

        # Every 500 steps, print
        if (i+1) % 500 == 0.0:
            print('Step {0:d} out of {1:d}'.format(i+1, N))

        # Every 2500 steps, print state
        if (i+1) % 2500 == 0.0:
            for pp in markov_chain.get_current_state():
                if pp.jump:
                    print('{0}: {1:.6f}'.format(pp.label, pp.get_value()))
            print('{0}: {1:.6f}'.format('logL', logLx))

        # Every Nsave steps, pickle chain value dict
        if outputfile is not None:
            if (i+1) % Nsave == 0.0:
                print('Saving preliminary results.')
                fout = open(outputfile, 'wb')
                vdpre = markov_chain.get_value_dict()
                for kk in vdpre.keys():
                    vdpre[kk] = vdpre[kk][:i]
                pickle.dump(vdpre, fout)
                fout.close()

        ########################
        # Principal components analysis
        #
        # Start after Npca iterations; update coefficient matrix every
        # NupdatePCA (None if no update is wished)
        ########################
        if i == Npca and usePCA:
            print('STARTING PRINCIPAL COMPONENT ANALYSIS')
            
            # Compute coefficient matrix M
            S, meanV = PCA.get_covariance_matrix(markov_chain, bi=BIpca,
                                                 iteration=i)
            M, w = PCA.run_pca(S, iteration=i)
            markov_chain.S = S
            markov_chain.M = M
            markov_chain.Ms = [M]
            markov_chain.Ss = [S]
            markov_chain.distancesM = []
            markov_chain.distancesS = []

            # Estimate proposal scale for principal components
            # sv = PCA.estimate_propscale_pca( C, M, S )
            sv = n.sqrt(w.copy())

            # Create list Parameter objects for Principal Components
            x = markov_chain.get_current_state_jumping_values()
            v = n.dot(M, x - meanV)

            U = []
            for ii in range(len(markov_chain.jumpind)):
                ui = objMCMC.Parameter(v[ii], None, label='U{}'.format(ii + 1),
                                       proposal_scale=sv[ii])
                U.append(ui)

            # Add list of PCs to Chain as _currentstatePCA attribute
            markov_chain._currentstatePCA = U

        if ((i+1 - Npca) % NupdatePCA == 0.0 and Npca < i < Nlastupdate) and \
                usePCA:

            print('Updating covariance matrix for PCA')

            # Compute coefficient matrix M
            if i > kPCA*NupdatePCA:
                BIupdatePCA = int(i - kPCA*NupdatePCA)
            else:
                BIupdatePCA = BIpca

            S, meanV = PCA.get_covariance_matrix(markov_chain, bi=BIupdatePCA,
                                                 iteration=i)
            M, w = PCA.run_pca(S, iteration=i)

            # Compute "distance" of coefficient matrices
            distM = n.linalg.norm(markov_chain.M - M, 'fro') / \
                n.linalg.norm(markov_chain.M, 'fro')
            distS = n.linalg.norm(markov_chain.S - S, 'fro') / \
                n.linalg.norm(markov_chain.S, 'fro')
            markov_chain.distancesM.append((i, distM))
            markov_chain.distancesS.append((i, distS))
            markov_chain.M = M
            markov_chain.S = S
            markov_chain.Ms.append(M)

            # If dist is smaller than mindiff, stop updating covariance
            # if dist <= mindiff:

            # Estimate proposal scale for principal components
            # sv = PCA.estimate_propscale_pca( C, M, S )
            sv = n.sqrt(w.copy())

            # Create list Parameter objects for Principal Components
            x = markov_chain.get_current_state_jumping_values()
            v = n.dot(M, x - meanV)

            U = []
            for ii in range(len(markov_chain.jumpind)):
                ui = objMCMC.Parameter(v[ii], None, label='U%d' % (ii + 1),
                                       proposal_scale=sv[ii])
                U.append(ui)

            """
            for ii in range(len(C.jumpind)):
                U[ii].set_value(v[ii])
            """

            # Add list of PCs to Chain as _currentstatePCA attribute
            markov_chain._currentstatePCA = U

        """
        if ((i+1 - Npca)%1000.0 == 0.0 and i > Npca):
            # Compute new covariance matrix S
            print('Updating jump size for PCA')
            S, meanV = PCA.get_covariance_matrix(C, iteration = i)
            sv = PCA.estimate_propscale_pca( C, M, S )
        """ 

        # Make copy of current state, that will be modified next.
        Y = X
                
        # Make proposal
        if i < Npca or not usePCA:
            if not AM:
                Y, jumppar = make_proposal(X, markov_chain.jumpind)

            else:
                Y = make_proposal2(markov_chain, i)
                jumppar = None

        else:
            Y, jumppar = make_proposal_pca(markov_chain, meanV, onebyone=True)

        # Small tweak for ciclical variables
        for cc in labeldict.keys():
            if 'omega' in cc or 'pomega' in cc or 'Omega' in cc or 'L0' in cc \
                    or 'M0' in cc or 'lambda0' in cc or 'spinorbit' in cc \
                    or 'rotangle' in cc:
                # Change the value of the angle variables
                labeldict[cc]._value = (labeldict[cc].get_value()) % 360.0

            elif 'incl' in cc or 'phi0' in cc:
                # Change value of the inclination
                if labeldict[cc].get_value() > 90.0:
                    labeldict[cc]._value = 180.0 - labeldict[cc].get_value()
                if labeldict[cc].get_value() < -90.0:
                    labeldict[cc]._value = -180.0 - labeldict[cc].get_value()
                    
        # Compute prior on proposal state
        priory, priorproby = compute_priors(priordict, labeldict)

        # If proposal state is not possible, reject and continue
        if priory == 0.0:
            
            """
            ### PRINT INFO (deactivated)
            outparams = []
            outstr = ''
            for kk in priorproby.keys():
                if priorproby[kk] == 0.0:
                    outparams.append(kk)
                    outstr = outstr+' '+kk
             
            print('Iteration %d: Prior 0.0.'%i)
            print('Blame: %s '%outstr)
            """
            markov_chain.reject_proposal(i, priorx, Lx, logLx, likedictx)

            # Check if need to update proposal scale.
            update_proposal_scale(jumppar)

            continue

        # Check if need to update proposal scale.
        update_proposal_scale(jumppar)

        plot_rv = False
        # Compute likelihood for proposal step
        try:
            Ly, logLy, likedicty = get_likelihood(Y, input_dict, datadict,
                                                  labeldict, plot_rv,
                                                  autocorrect)

        except RVgaussianFitError:
            # print(rvfite)
            # print('Input parameters: %f, %f, %f'%(rvfite.contrast, rvfite.rv0,
            #                                        rvfite.sigma)
            #      )
            markov_chain.TrackError.append(i)
            markov_chain.reject_proposal(i, priorx, Lx, logLx, likedictx)
            continue

        except (EvolTrackError, OutofIsochroneError):
            # print('One of the requested stars cannot be constructed because it
            # is outside the limits of the evolution tracks.')
            # Instead of printing, add number of step to TrackError list
            markov_chain.TrackError.append(i)
            markov_chain.reject_proposal(i, priorx, Lx, logLx, likedictx)
            continue

        except RuntimeError:
            print('Runtime Error. Probably and error in Qhull. Rejecting step.')
            markov_chain.reject_proposal(i, priorx, Lx, logLx, likedictx)
            continue

        except EBOPparamError as eboperr:
            print('Parameter outside limits for JKTEBOP; '
                  'error message: {}'.format(eboperr))
            markov_chain.reject_proposal(i, priorx, Lx, logLx, likedictx)
            continue
        
        except ValueError as verr:
            print('Value Error. Rejecting step. '
                  'Message: {}'.format(verr))
            markov_chain.reject_proposal(i, priorx, Lx, logLx, likedictx)
            continue
            
        # Any other error, just reject the step.
        except:
            traceback.print_exc()
            markov_chain.reject_proposal(i, priorx, Lx, logLx, likedictx)
            continue
            
        # Use computed likelihood and priors to produce Metropolis ratio
        try:
            r = priory/priorx*exp(beta*(logLy - logLx))

        # See this in DETAIL!!
        except OverflowError:
            try:
                # If Metropolis ratio overflows, compute rather log(r)
                # and check that the ratio is greater than 1
                logr = log(priory) - log(priorx) + beta*(logLy - logLx)
                if logr >= 0.0:
                    accept = True
            except ValueError:
                # If there's a problem stop execution !!
                raise ValueError('Fatal Error: can not compute Metropolis '
                                 'ratio.\n Priors < 0.0??')

        else:
            # If r is > 1.0, accept step right away.
            if r >= 1.0:
                accept = True
            # If not decide whether to accept step or not
            elif scipy.random.random() <= r:
                accept = True
            else:
                accept = False

        if accept:
            # Update parameter in case of one-by-one proposals
            if jumppar is not None and not AM:
                jumppar.Naccepted += 1
                jumppar.Naccepted_since_update += 1

            markov_chain.accept_proposal(i, priory, Ly, logLy, likedicty)
            Lx, logLx, likedictx = Ly, logLy, likedicty.copy()
            priorx = priory
            X = markov_chain.get_current_state()

        else:
            # print jumppar.get_value()
            markov_chain.reject_proposal(i, priorx, Lx, logLx, likedictx)

    print('Total number of accepted steps: {:d}'.format(markov_chain.Naccept))
    print('Running time (without initialization): '
          '{0:.3f} hours'.format((time.time() - ti)/3600))
    return markov_chain


####
# Function to randomly choose starting point in priors
####
def pick_random_point(labeldict, priordict):
    """
    Set the values of the parameters in labeldict randomly,
    by sampling the prior distributions in priordict
    """
    print('Starting chain at random point in priors.')
    # Iterate over all priors
    for kk in priordict.keys():
        
        if isinstance(priordict[kk], rv_frozen) and labeldict[kk].jump == 1:
            try:
                val = priordict[kk].rvs()
            except ValueError:
                raise ValueError('Prior error for parameter: '+kk)
        else:
            continue
            # Implement custom prior!

        """
        elif isinstance(priordict[kk], list()):
        plist = priordict[kk]
        val = plist[0].get_random_point()
        """

        labeldict[kk].set_value(val)
    return


###
# FUNCTIONS FOR PROPOSAL MAKING
###
def make_proposal(state, jumpind):
    """
    Construct proposal state of chain from given state.

    Uses prescription from Ford (ApJ 642 505, 2006)
    (i.e. change only one parameter at the time)

    state: current chain state (an array of Parameter instances)
    
    jumpind: list with indices of X containing jump parameters
    """
    # Decide randomly which parameter to jump
    parind = random.choice(jumpind)
    par = state[parind]
                           
    # Construct proposal state
    if par.propfunc == 'gaussian':
        # print par.propscale
        
        if par.propscale <= 0.0:
            print('Warning! Proposal scale is zero for parameter {}'
                  .format(par.label))
            newvalue = par.get_value()
        else:
            newvalue = par.get_value() + scipy.random.normal(0.0, par.propscale)

    elif hasattr(par.propfunc, '__call__'):
        newvalue = par.propfunc(par.get_value())

    else:
        raise NameError('Proposal distribution not known')

    # Assign new value to parameter
    par.set_value(newvalue)

    # Update parameter instance
    par.Nstep += 1
    par.Nstep_since_update += 1

    for otherpar in state:
        if otherpar != par:
            otherpar.set_value(otherpar.get_value())

    return state, par


def make_proposal_pca(chain, mean_values, onebyone=False):
    """
    Transform current state current_state to PCA space using matrix M, then
    perform jump in PCA space for all parameters at the same time and transform
    back to normal space.
    """
    current_state = chain.get_current_state()
    current_state_pca = chain.get_current_state(PCA=True)
    
    # Change one parameter at the time
    if onebyone:
        # Decide randomly which parameter to jump
        parind = random.choice(range(len(chain.jumpind)))
        # jumppar = current_state_pca[0]
        jumppar = current_state_pca[parind]
                           
        # Construct proposal state
        if jumppar.propfunc == 'gaussian':
            if jumppar.propscale <= 0.0:
                print('Warning! Proposal scale is zero for parameter {}'
                      .format(jumppar.label))
                newvalue = jumppar.get_value()
            else:
                newvalue = jumppar.get_value() + \
                    scipy.random.normal(0.0, jumppar.propscale)
            
        elif hasattr(jumppar.propfunc, '__call__'):
            newvalue = jumppar.propfunc(jumppar.get_value())
        
        else:
            raise NameError('Proposal distribution not known')

        # Assign new value to parameter
        jumppar.set_value(newvalue)

        # Update parameter instance
        jumppar.Nstep += 1
        jumppar.Nstep_since_update += 1

        # For all other parameters, keep same value
        for otherpar in current_state_pca:
            if otherpar != jumppar:
                otherpar.set_value(otherpar.get_value())
            
    # Change all parameters at the same time
    else:
        jumppar = None

        # Produce proposal state for all jumping parameters
        for par in current_state_pca:
            vprime = par.get_value() + scipy.random.normal(0.0, par.propscale)
            par.set_value(vprime)
    ##

    # Get modified values of PCA
    vprime = chain.get_current_pca_values()
    # Revert back to normal space
    # XX = n.dot((markov_chain.M).T, vprime.T)
    # YY = n.dot((markov_chain.M).T, vprime.T) + meanV
    # print XX.shape, YY.shape
    y = n.dot(chain.M.T, vprime) + mean_values
    
    # Save values in state
    for i, ind in enumerate(chain.jumpind):
        current_state[ind].set_value(y[i])

    return current_state, jumppar


def make_proposal2(markov_chain, i):
    """
    Construct proposal state of markov_chain from given state current_state.

    Uses prescription from Ford (ApJ 642 505, 2006)
    (i.e. change only one parameter at the time)

    current_state: current markov_chain state (an array of Parameter instances)
    
    jumpind: list with indices of current_state containing jump parameters
    """
    current_state = markov_chain.get_current_state()
    xx = markov_chain.get_current_state_jumping_values()
    
    sd = 2.4**2/len(xx)
    eps = 1e-4

    # Update mean of parameters
    if i > 1:
        markov_chain._meanX[i - 1] = markov_chain._meanX[i - 2]*(i - 1)/i + xx/i

    # Update covariance matrix
    if i <= 1:
        S = markov_chain._S0

    elif i < 5000:
        S = sd * n.cov(markov_chain.values[:i, markov_chain.jumpind],
                       rowvar=0) + sd * eps * n.identity(len(xx))

    else:
        xx2 = markov_chain._meanX[i - 2].reshape((len(xx), 1))
        xx1 = markov_chain._meanX[i - 1].reshape((len(xx), 1))

        pxx2 = xx2 * xx2.T
        pxx1 = xx1 * xx1.T

        xxx = xx.reshape((len(xx), 1))
        pxx = xxx * xxx.T

        S = (i - 2.)/(i - 1.) * markov_chain._S + \
            sd/(i-1.)*((i - 1.)*pxx2 - i*pxx1 + pxx + eps*n.identity(len(xx)))
    
    markov_chain._S = S.copy()

    # Produce multivariate proposal
    dX = scipy.random.multivariate_normal(xx.reshape((len(xx),)), S)

    for ii, j in enumerate(markov_chain.jumpind):
        # Assign new value to parameter
        current_state[j].set_value(dX[ii])

    return current_state


def update_proposal_scale(par, target_rate=0.25):
    """
    Check if proposal scale has to be updated, following
    Ford (ApJ 642 505, 2006), for a given parameter

    If so, update it.

    par: Parameter instance
    """
    
    acceptance_rate = par.Naccepted_since_update/float(par.Nstep_since_update)
    number_steps_since_update = par.Nstep_since_update
    update_interval = par.update_interval
    betamu = par.propscale

    # Compute variance of estimated acceptance rate acceptance_rate.
    # This is: psi0*(1 - psi0)/number_steps_since_update
    # If the squared difference (acceptance_rate - psi0)**2 is larger
    # than a few times this value, then the scale must be updated

    update_stat = update_interval * target_rate * (
        1 - target_rate) / number_steps_since_update

    # Update step size
    if (acceptance_rate - target_rate) ** 2 > update_stat and \
        abs(acceptance_rate - target_rate) > 0.1 * target_rate and \
            par.Nstep_since_update >= 100:
        
        # VERSION ROLO #
        if acceptance_rate > 2*target_rate:
            k = (acceptance_rate/target_rate)**2.0
        elif target_rate < acceptance_rate <= 2*target_rate:
            k = acceptance_rate/target_rate
        elif 0.5*target_rate < acceptance_rate <= target_rate:
            k = (acceptance_rate/target_rate)**0.5
        elif 0.2*target_rate < acceptance_rate <= 0.5*target_rate:
            k = (acceptance_rate/target_rate)**1.5
        elif 0.1*target_rate < acceptance_rate <= 0.2*target_rate:
            k = (acceptance_rate/target_rate)**2.0
        else:
            k = 0.01  # Never reduce k by more than 100

        """
        ### VERSION FORD ###
        if acceptance_rate > 0.5*psi0:
            k = (acceptance_rate/psi0)**0.5
        elif acceptance_rate > 0.2*psi0 and acceptance_rate <= 0.5*psi0:
            k = (acceptance_rate/psi0)**1.5
        elif acceptance_rate > 0.1*psi0 and acceptance_rate <= 0.2*psi0:
            k = (acceptance_rate/psi0)**2.0
        else:
            k = 0.01 # Never reduce k by more than 100
        """

        # Check if scale has being increased and then decreased or viceversa.
        # If so, increase mean update interval.
        if par.lastk is not None:
            if (par.lastk - 1)*(k - 1) < 0.0:
                par.update_interval += 1.0
        
        # Update proposal scale
        if par.propscale != 0.0:
            par.propscale *= k
        elif k > 1.0:
            par.propscale = par.get_value()*0.1

        # Set number of steps and accepted steps since last update to zero
        par.Nstep_since_update = 0.0
        par.Naccepted_since_update = 0.0       

        # Update value of lastk used
        par.lastk = k

        # Update dictionary with proposal scales
        par.update_jump_history()

        # Since changing scale "breaks" the chain, record step number
        # change_step.append(i)

        if betamu*k < 1e-5 or betamu < 1e-5:
            print('Scale adjusted for parameter {0};\t\t acceptence rate: '
                  '{1:.2f}%; old scale {2:.2e}; new scale {3:.2e}; '
                  'Value = {4:.3f}'.format(par.label, acceptance_rate*100,
                                           betamu, betamu*k, par.get_value())
                  )
        else:
            print('Scale adjusted for parameter {0};\t\t acceptence rate: '
                  '{1:.2f}%; old scale {2:.6f}; new scale {3:.6f}; '
                  'Value = {4:.3f}'.format(par.label, acceptance_rate*100,
                                           betamu, betamu*k, par.get_value())
                  )
    return


###
# FUNCTIONS FOR LIKELIHOOD
###
def get_likelihood(state, input_dict, datadict, labeldict, autocorrect,
                   LM=False):
    """
    Compute the likelihood of chain state.
    Uses variables from configuration file.
    
    The keyword 'distribution' is the distribution function of the data.
    The possible options are:

    'distrib_params' is a list containing the parameters for the distribution
    function, if any.
    """

    # Construct dictionary for Object Builder
    input_dict = state_deconstructor(state, input_dict)

    # Build objects
    objects = ObjectBuilder.ObjectBuilder(input_dict)
    
    ###
    # HARD CODED magnitude LIMITATION TO BLENDS !! tremenda chapuza ...
    ###
    # print('checkblendmag = '+str(checkblendmag))
    if checkblendmag:
        checkmag = False
    
        for obj in objects:

            # For PiB and TRIPLE
            if isinstance(obj, ac.Triple):
                # Get primary magnitude
                mag_target = obj.object1.get_mag('Johnson-V')

                if isinstance(obj.object2, ac.PlanSys):
                    # Get secondary magnitude
                    mag_binary = obj.object2.star.get_mag('Johnson-V')
                    checkmag = True
                elif isinstance(obj.object2, ac.IsoBinary):
                    # Get secondary magnitude
                    mag_binary = obj.object2.star1.get_mag('Johnson-V')
                    checkmag = True
                
            # For BEB and BTP
            elif isinstance(obj, ac.Target):
                # Get target magnitude
                mag_target = obj.get_mag('Johnson-V')
                checkmag = True
            
            elif isinstance(obj, ac.IsoBinary):
                # Get binary magnitude
                mag_binary = obj.star1.get_mag('Johnson-V')
                checkmag = True
            
            elif isinstance(obj, ac.PlanSys):
                # Get binary magnitude
                mag_binary = obj.star.get_mag('Johnson-V')
                checkmag = True

            else:
                raise NameError('Object class not recognised.')

        # print('checkmag = '+str(checkmag))
        if checkmag:
            # If binary is not fainter than 1 magnitude, raise Exception.
            try:
                if mag_binary - mag_target < 1.0:
                    raise ValueError('Binary is too bright!')
            except (UnboundLocalError, NameError):
                pass
    ###
    ###
    ###

    # Dictionary that will contain all computed likelihoods
    likeout = {}

    # To know if CoRoT colors have already been computed
    has_computed_colors = False

    # Prepare residuals to return if using Least-Squares fit
    res = n.zeros(0)

    # Prepare list with the names of the RV observables
    rvdiags = ['RV', 'CTRS', 'FWHM', 'BIS', 'Vspan', 'Wspan',
               'BiGauss', 'Vasy']
    
    # Iterate for all datasets
    for key in datadict.keys():
        if datadict[key]['type'] == 'RV':
            
            # Get information of instrument
            spectro = datadict[key]['instrument']
            mask = datadict[key]['mask']
            
            # Get value of instrument offsets
            offsetdict = {}

            for rvobs in rvdiags:
                if (key + '_' + rvobs + 'offset') in labeldict:
                    offsetdict[rvobs] = {}
                    offsetdict[rvobs]['offset'] = \
                        labeldict[key+'_'+rvobs+'offset'].get_value()

            # Get value of jitter for each observable
            jitterdict = {}

            for rvobs in rvdiags:
                if (key + '_' + rvobs + 'offset') in labeldict:
                    jitterdict[rvobs] = get_jitter(datadict[key], key,
                                                   labeldict,
                                                   observable=rvobs)

            # Prepare dictionary with info for PASTIS_RV
            rvdict = {'spectro': spectro, 'mask': mask}
            rvdict.update(offsetdict)

            # Check if need to compute oversampled RVcurve
            if 'overtime' in datadict[key] and len(datadict[key]['overtime']) \
                    != len(datadict[key]['data']['time']):
                # Compute theoretical OVERSAMPLED RVcurve
                oversampled_rv_output = PASTIS_RV(datadict[key]['overtime'],
                                                  rvdict,
                                                  *objects)

                # Rebin model flux to original sampling rate
                rvoutput = {}
                if 'texp' in datadict[key]:
                    for rvkey in oversampled_rv_output.keys():
                        rvoutput[rvkey] = \
                            rebin_texp(datadict[key]['overtime'],
                                       oversampled_rv_output[rvkey],
                                       datadict[key]['data']['time'],
                                       datadict[key]['texp'])

                elif 'texp' in datadict[key]['data']:
                    for rvkey in oversampled_rv_output.keys():
                        rvoutput[rvkey] = \
                            rebin_texp(datadict[key]['overtime'],
                                       oversampled_rv_output[rvkey],
                                       datadict[key]['data']['time'],
                                       datadict[key]['data']['texp'])

                else:
                    for rvkey in oversampled_rv_output.keys():
                        rvoutput[rvkey] = rebin(oversampled_rv_output[rvkey],
                                                datadict[key]['sampling'])
            else:
                # Compute theoretical RV curves and likelihood
                rvoutput = PASTIS_RV(datadict[key]['data']['time'],
                                     rvdict,
                                     *objects
                                     )

            # Compute likelihood for each RV observable
            for ii, jj in zip(rvdiags,
                              ['vrad', 'cont', 'fwhm', 'bis',
                               'vspan', 'wspan', 'bigauss', 'vasy']):
                
                if ii not in rvdict:
                    continue

                L, logL = likelihood(datadict[key]['data'][jj],
                                     datadict[key]['data']['s'+jj],
                                     rvoutput[ii],
                                     jitter=jitterdict[ii])
                
                likeout[key+'_'+ii] = [L, logL]

                if LM:
                    resi = (datadict[key]['data'][jj] - kk) / \
                        datadict[key]['data']['s'+jj]
                    res = n.concatenate((res, resi))

        if datadict[key]['type'] == 'PHOT':
            try:
                filtre = datadict[key]['filter']
            except KeyError:
                filtre = key

            # Get value of contamination and Foot
            cont = labeldict[key+'_contamination'].get_value()
            foot = labeldict[key+'_foot'].get_value()
            
            # Get value of ephemeris offset, if it exists
            try:
                dt0 = labeldict[key+'_dt0'].get_value()
            except KeyError:
                dt0 = 0.0

            # Get value of jitter for this dataset
            jitter = get_jitter(datadict[key], key, labeldict)

            # If filter is one of CoRoT colors
            if (filtre == 'CoRoT-R' or filtre == 'CoRoT-G' or
                    filtre == 'CoRoT-B') and not has_computed_colors:
                for key2 in datadict.keys():
                    if 'filter' in datadict[key2]:
                        if datadict[key2]['filter'] == 'CoRoT-R':
                            meanR = datadict[key2]['MeanFlux']
                            contR = labeldict[key2+'_contamination'].get_value()

                        if datadict[key2]['filter'] == 'CoRoT-G':
                            meanG = datadict[key2]['MeanFlux']
                            contG = labeldict[key2+'_contamination'].get_value()

                        if datadict[key2]['filter'] == 'CoRoT-B':
                            meanB = datadict[key2]['MeanFlux']
                            contB = labeldict[key2+'_contamination'].get_value()

                # Add limbdarkening weights and filters for CoRoT colors
                if 'global_spectrum' not in SED.__dict__:
                    SED.compute_global_spectrum(*objects)
                corot_colors(meanR, meanG, meanB, contR, contG, contB)
                has_computed_colors = True
                
            # Check if need to compute oversampled lightcurve
            if 'overtime'in datadict[key] and len(datadict[key]['overtime']) \
                    != len(datadict[key]['data']['time']):
                # Compute theoretical OVERSAMPLED lightcurve
                osft = PASTIS_PHOT(datadict[key]['overtime'], filtre,
                                   datadict[key]['is_phase'], cont, foot, dt0,
                                   *objects)

                # Rebin model flux to original sampling rate
                if 'texp' in datadict[key]:
                    ft = rebin_texp(datadict[key]['overtime'], osft,
                                    datadict[key]['data']['time'],
                                    datadict[key]['texp'])
                elif 'texp' in datadict[key]['data']:
                    ft = rebin_texp(datadict[key]['overtime'], osft,
                                    datadict[key]['data']['time'],
                                    datadict[key]['data']['texp'])
                else:
                    ft = rebin(osft, datadict[key]['sampling'])
            else:
                # Compute theoretical light curves and likelihood
                ft = PASTIS_PHOT(datadict[key]['data']['time'], filtre,
                                 datadict[key]['is_phase'], cont, foot, dt0,
                                 *objects)

            # TESTING AUTOCORRECTION
            if autocorrect:

                tf, ff, efc = PASTIS_CORRECT(datadict[key]['data']['time'],
                                             datadict[key]['data']['flux'],
                                             datadict[key]['data']['sflux'],
                                             *objects
                                             )
                if tf is None:
                    # No transit
                    L = 0.0
                    logL = -n.inf

                else:
                    ft = PASTIS_PHOT(tf, filtre, datadict[key]['is_phase'],
                                     cont, foot, dt0, *objects)

                    print(len(ff), len(efc), len(ft))
                    L, logL = likelihood(ff, efc, ft)
            ###
            else:
                L, logL = likelihood(datadict[key]['data']['flux'],
                                     datadict[key]['data']['sflux'], ft,
                                     jitter=jitter)

                if LM:
                    error = n.sqrt(datadict[key]['data']['sflux']**2 +
                                   jitter**2)
                    # error = datadict[key]['data']['sflux']
                    resi = (datadict[key]['data']['flux'] - ft)/error
                    res = n.concatenate((res, resi))

            likeout[key] = [L, logL]

        if datadict[key]['type'] == 'SED':

            # Compute theoretical magnitudes and likelihood
            magst = PASTIS_SED(datadict[key]['data']['band'], *objects)

            # Get value of jitter for this dataset
            jitter = get_jitter(datadict[key], key, labeldict)

            L, logL = likelihood(datadict[key]['data']['mags'],
                                 datadict[key]['data']['smags'], magst,
                                 jitter=jitter)

            if LM:
                resi = (datadict[key]['data']['mags'] - magst) / \
                    datadict[key]['data']['smags']
                res = n.concatenate((res, resi))

            likeout[key] = [L, logL]
    ##
            
    # Combine likelihoods
    L = n.array(list(likeout.values()))[:, 0]
    logL = n.array(list(likeout.values()))[:, 1]

    Lcomb = n.prod(L)
    logLcomb = n.sum(logL)

    # CHANGED WHEN GOING TO PACKAGE FORMULATION ####
    if 'global_spectrum' in SED.__dict__:
        del SED.global_spectrum

    if LM:
        return res
    else:
        return Lcomb, logLcomb, likeout


def likelihood(y, ey, yt, distribution='normal', jitter=0.0, distrib_params=()):

    # Add jitter term
    ey = n.sqrt(ey*ey + jitter*jitter)
    
    # Multiply error by jitter term
    # ey = jitter*ey
    
    res_norm = (y - yt)/ey

    if distribution == 'normal':
        dist_func = stats.norm.pdf  # Be careful with the normalization check!
        logp = -0.5*(res_norm*res_norm + n.log(2.0*pi*ey*ey))
        
    elif distribution == 'poisson':
        dist_func = stats.poisson.pmf

    else:
        raise NameError('Distribution function not found.')

    if len(distrib_params) > 0:
        p = dist_func(res_norm, *distrib_params)
    else:
        p = dist_func(res_norm)
    
    # If the data are independent, the likelihood is the product of all
    # probabilities
    """
    if n.any(n.equal(p, 0)):
        print('Warning! Likelihood out of resolution.')
        return 0.0, -1e308
    """
    
    L = n.prod(p)
    logL = n.sum(logp)

    return L, logL
