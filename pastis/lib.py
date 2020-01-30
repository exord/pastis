"""
Module containing useful functions to link PASTIS MCMC posterior samples with
the bayev package.
"""
import os
import pickle
import importlib
import numpy as np

from .paths import resultpath, configpath
from .exceptions import EBOPparamError, EvolTrackError
from . import MCMC
from .DataTools import readdata


def read_pastis_file(target, simul, pastisfile=None):
    """Read configuration dictionary."""
    if pastisfile is None:
        # Get input_dict
        configname = os.path.join(configpath, target,
                                  target + '_' + simul + '.pastis')
    else:
        configname = pastisfile

    try:
        f = open(configname)
    except IOError:
        raise IOError('Configuration file {} not found!'.format(configname))

    dd = pickle.load(f)
    f.close()

    return dd


def get_datadict(target, simul, pastisfile=None):

    config_dicts = read_pastis_file(target, simul, pastisfile)
    return readdata(config_dicts[2])[0]


def get_priordict(target, simul, pastisfile=None):

    config_dicts = read_pastis_file(target, simul, pastisfile)
    return MCMC.priors.prior_constructor(config_dicts[1], {})


def get_posterior_samples(target, simul, mergefile=None,
                          suffix='_Beta1.000000_mergedchain.dat'):

    if mergefile is None:
        mergepath = os.path.join(resultpath, target,
                                 target + '_' + simul + suffix)
    else:
        mergepath = mergefile

    f = open(mergepath, 'r')
    vdm = pickle.load(f)
    f.close()

    return vdm


def pastis_init(target, simul, posteriorfile=None, datadict=None,
                pastisfile=None):

    # Read configuration dictionaries.
    configdicts = read_pastis_file(target, simul, pastisfile)

    infodict, input_dict = configdicts[0], configdicts[1].copy()

    # Construct prior instances
    priordict = get_priordict(target, simul, pastisfile=pastisfile)

    # Read data dictionary.
    if datadict is None:
        datadict = get_datadict(target, simul, pastisfile=pastisfile)

    # Obtain PASTIS version the merged chain was constructed with.
    vdm = get_posterior_samples(target, simul, mergefile=posteriorfile)
    modulename = vdm.__module__.split('.')[0]

    # Import the correct PASTIS version used to construct a given posterior
    # sample
    pastis = importlib.import_module(modulename)

    # To deal with potential drifts, we need initialize to fix TrefRV.
    pastis.initialize(infodict, datadict, input_dict)

    # import PASTIS_rhk.MCMC as MCMC
    # MCMC.PASTIS_MCMC.get_likelihood

    importlib.import_module('.MCMC.PASTIS_MCMC', package=pastis.__name__)
    importlib.import_module('.AstroClasses', package=pastis.__name__)
    importlib.import_module('.ObjectBuilder', package=pastis.__name__)
    importlib.import_module('.models.RV', package=pastis.__name__)

    reload(pastis.AstroClasses)
    reload(pastis.ObjectBuilder)
    reload(pastis.models.RV)
    reload(pastis.MCMC.PASTIS_MCMC)

    return priordict, datadict, input_dict


def pastis_loglike(samples, paramdict, jumppardict, input_dict, datadict):
    """
    A wrapper to run the PASTIS.MCMC.get_likelihood function.

    Computes the loglikelihood on a series of points given in samples using
    PASTIS.MCMC.get_likelihood.

    :param np.array samples: parameter samples on which to compute log
    likelihood. Array dimensions must be (np x k), where *np* is the number of
    samples and *k* is the number of model parameters.

    :param dict jumppardict: dictionary containing the Parameter 
    instances for all jump parameters. Keys must be the parameter names. 

    :return:
    """
    # Prepare output arrays
    loglike = np.zeros(samples.shape[0])

    for i, s in enumerate(samples):

        for parameter_index, param_name in enumerate(jumppardict):

            if paramdict[param_name].jump == 0:
                continue
            
            # Modify input_dict
            paramdict[param_name].set_value(s[parameter_index])

        chain_state = list(paramdict.values())
        
        try:
            # Compute likelihood for this state
            ll, loglike[i], likeout = MCMC.PASTIS_MCMC.get_likelihood(
                chain_state, input_dict, datadict, paramdict, False,
                False)

        except (ValueError, RuntimeError, EvolTrackError, EBOPparamError):
            print('Error in likelihood computation, setting lnlike to -np.inf')
            loglike[i] = -np.inf
            pass

    return loglike


def pastis_logprior(samples, paramdict, pdfdict):
    """
    A wrapper to run the PASTIS.MCMC.get_likelihood function.

    Computes the loglikelihood on a series of points given in samples using
    PASTIS.MCMC.get_likelihood.

    :param np.array samples: parameter samples on which to compute log
    likelihood. Array dimensions must be (np x k), where *np* is the number of
    samples and *k* is the number of model parameters.

    :param list params: parameter names.

    :return:
    """
    # Prepare output arrays
    logprior = np.zeros(samples.shape[0])

    for i, s in enumerate(samples):

        for parameter_index, param_name in enumerate(paramdict):
            
            # Modify input_dict
            paramdict[param_name].set_value(s[parameter_index])

        # Construct chain state
        #chain_state, labeldict = MCMC.tools.state_constructor(input_dict)

        # Compute prior distribution for this state
        prior_probability = MCMC.priors.compute_priors(pdfdict, 
                                                       paramdict)[0]
        logprior[i] = np.log(prior_probability)

    return logprior


# def lnpost(samples, params, input_dict, datadict, priordict=None):

def lnpost(samples, input_dict, paramdict, jumpparamdict, 
           datadict, priordict):
    
    x = samples.copy()
    x = np.reshape(x, (1, len(x)))
    
    lnprior = pastis_logprior(x, jumpparamdict, priordict)
    # Do not compute likelihood for points outside prior.
    if lnprior[0] == -np.inf:
        return -np.inf
    lnlike = pastis_loglike(x, paramdict, jumpparamdict, input_dict, 
                            datadict)

    return (lnprior + lnlike)[0]
