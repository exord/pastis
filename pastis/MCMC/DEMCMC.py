import numpy as n
import scipy
import random
import types
import imp
import time
from math import log10, e
from scipy import stats, interpolate
from scipy.stats.distributions import rv_frozen

# Intra-package imports
from tools import state_constructor, state_deconstructor
from priors import compute_priors, prior_constructor
import Objects_MCMC as objMCMC
import PCA

# Imports from upper package level
#from .. import inputdicts
from .. import ObjectBuilder
from ..AstroClasses import *
from ..tools import rebin
from ..photometry import corot_colors
from ..models import SED
from ..models import RVgaussianFitError
from ..models.PHOT import PASTIS_PHOT
from ..models.SED import PASTIS_SED
from ..models.RV import PASTIS_RV

from .PASTIS_MCMC import *

def DEMCMC(input_dict, datadict, customprior_dict, N, Nchains,
           randomstart = True, **kwargs):
    """
    Implementation of DE-MC code
    (Caj. J. F. Ter Braak, Stat Comput (2006) 16:239-249)
    """

    autocorrect = kwargs.pop('autocorrect', False)
    gamfactor = kwargs.pop('gamma', 0.4)
    beta = kwargs.pop('beta', 1.0)
    randfactor = kwargs.pop('randomfactor', 0.0)
    comment = kwargs.pop('comment', '')
    if comment != '':
        comment = comment+'_'

    resume = kwargs.pop('resume', None)
    save =  kwargs.pop('save', True)
    
    ### INITIALISATION
    
    # Construct initial state from info in input_dict
    X, labeldict = state_constructor(input_dict)

    # Prepare dictionary with prior functions
    priordict = prior_constructor(input_dict, customprior_dict)

    jumpind = []
    for i, xx in enumerate(labeldict.keys()):
        if labeldict[xx].jump:
            jumpind.append(i)


    ## DEFINE GRAND STATE XX
    XX = n.zeros((N, Nchains, len(jumpind) + 3))

    priorx = n.zeros(Nchains)
    logLx =  n.zeros(Nchains)

    accepted = n.zeros(Nchains)

    ti = time.time()
    for j in range(Nchains):

        # If requested, start chains at random point in prior
        if randomstart:
            pick_random_point(labeldict, priordict)

        for kk, ii in enumerate(jumpind):

            if resume == None:
                XX[0, j, kk] = labeldict[labeldict.keys()[ii]].get_value()

            else:
                XX[0, j, kk] = resume[-1, j, kk]
                labeldict[labeldict.keys()[ii]].set_value(XX[0, j, kk])


        priorx[j], priorprobx = compute_priors(priordict, labeldict)

        # Compute likelihood for initial state
        haspassed = False

        while not haspassed:
            try:
                Lx, logLx[j], likedictx = get_likelihood(X,
                                                         input_dict,
                                                         datadict,
                                                         labeldict,
                                                         False, autocorrect)

            except RVgaussianFitError as rvfite:
                print(rvfite)
                print('Input parameters: %f, %f, %f'%(rvfite.contrast, 
                                                      rvfite.rv0,
                                                      rvfite.sigma)
                      )

            except (EvolTrackError, OutofIsochroneError):
                print('Something went wrong. Initial objects cannot be constructed due to limitations in the stellar evolution tracks.')
                pick_random_point(labeldict, priordict)

            except RuntimeError:
                print('Something went wrong. Initial objects cannot be constructed due to a Runtime error (possibly Qhull!).')
                pick_random_point(labeldict, priordict)

            except EBOPparamError as eboperr:
                print('Parameter outside limits for JKTEBOP; error message: %s'\
                          %eboperr.message)
                pick_random_point(labeldict, priordict)

            else:
                haspassed = True

    XX[0, :, -3] = priorx
    XX[0, :, -2] = logLx
    XX[0, :, -1] = logLx*log10(e) + n.log10(priorx)


    N = int(N)
    for i in xrange(1, N):

        # Every 500 steps, print
        if (i+1)%500 == 0.0:
            print('Step %d out of %d'%(i+1, N))
                


        ####
        # Make proposal for each chain
        ####
        for j in range(Nchains):
            
            # Every 500 steps, print
            if (i+1)%500 == 0.0 or ii == 1:
                print('Chain %d: posterior = %f'%(j, XX[i - 1, j, -1]))
            
                      
            reject = False

            """
            ### Sort the chains according to their merit
            ind = n.argsort(XX[i - 1, :, -1])[::-1]
            # Choose two chains from the best half at random
            indbest = ind[: Nchains/2.]
            random.shuffle(indbest)
            c1, c2 = indbest[:2]
            """
           
            ### Pick the two chains that will construct the jump vector
            chains = n.concatenate((n.arange(j), n.arange(j + 1, Nchains)))
            random.shuffle(chains)
            c1, c2 = chains[:2]

            ## Construct proposal vector
            xp = XX[i - 1, j, :-3] + gamfactor*(XX[i - 1, c1, :-3] - XX[i - 1, c2, :-3]) + n.random.randn(len(XX[i - 1, j, :-3]))*XX[i - 1, j, :-3]*randfactor

            for kk, ii in enumerate(jumpind):
                labeldict[labeldict.keys()[ii]].set_value(xp[kk])

                
            # Compute prior on proposal state
            priory, priorproby = compute_priors(priordict, labeldict)

            # If proposal state is not possible, reject and continue
            if priory == 0.0:
                XX[i, j, :] = XX[i - 1, j, :]
                continue

                
            ### Compute likelihood on proposal step
            try:
                Ly, logLy, likedicty = get_likelihood(X, input_dict, datadict,
                                                      labeldict, False,
                                                      autocorrect)

            except RVgaussianFitError as rvfite:
                print(rvfite)
                print('Input parameters: %f, %f, %f'%(self.contrast, self.rv0,
                                                      self.sigma)
                      )
                reject = True

            except (EvolTrackError, OutofIsochroneError):
            #print('One of the requested stars cannot be constructed because it is outside the limits of the evolution tracks.')
	    ## Instead of printing, add number of step to TrackError list
                reject = True

            except RuntimeError:
                print('Runtime Error. Probably and error in Qhull. Rejecting step.')
                reject = True
	    
            except EBOPparamError as eboperr:
                
                #print('Parameter outside limits for JKTEBOP; error message: %s'\
                #          %eboperr.message)
                reject = True
            
            else:
                reject = False
                

            if reject:
                XX[i, j, :] = XX[i - 1, j, :]
                continue

            try:
                r = priory/XX[i - 1, j, -3]*exp(beta*(logLy - XX[i - 1, j, -2]))

            except OverflowError:
                try:
                    # If Metropolis ratio overflows, compute rather log(r)
                    # and check that the ratio is greater than 1
                    logr = log(priory) - log(XX[i - 1, j, -3]) \
                        + beta*(logLy - XX[i - 1, j, -2])

                    if logr >= 0.0:
                        accept = True

                except ValueError:
                    # If there's a problem stop execution !!
                    print(priory, XX[i - 1, j, -3])
                    print('Chain: %d; Step: %d'%(j, i))
                    return XX, priory, logLy
                    raise ValueError('Fatal Error: can not compute Metropolis ratio.\n Priors < 0.0??')


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
                XX[i, j, :-3] = xp
                XX[i, j, -3] = priory
                XX[i, j, -2] = logLy
                XX[i, j, -1] = logLy*log10(e) + n.log10(priory)
                accepted[j] = accepted[j] + 1
            else:
                XX[i, j, :] = XX[i - 1, j, :]

                

    print('Running time: %.4f hours'%((time.time() - ti)/3600.0))

    """
    ## Create dictionary
    vd = {}

    for i, kk in enumerate(jumpind):
        vd[n.array(labeldict.keys())[jumpind]] = XX[:, :, i]
    vd['prior'] = XX[:, :, -3]
    vd['likelihood'] = XX[:, :, -2]
    """
    
    if save:
        import datetime
        dt = datetime.datetime.isoformat(datetime.datetime.now())
        f = open('/data/PASTIS/resultfiles/testDEMC/%stestDEMC_nchains%d_gamma%.4f_%s.dat'%(comment, Nchains, gamfactor, dt), 'w')
        import pickle
        pickle.dump([XX[::10,:, :],
                     accepted, n.array(labeldict.keys())[jumpind]], f)
        f.close()
    
    return XX, accepted, n.array(labeldict.keys())[jumpind]
