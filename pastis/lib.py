"""
Module containing useful functions to link PASTIS MCMC posterior samples with
the bayev package.
"""
import os
import pickle
import importlib
import numpy as np

from . import resultpath, configpath, EBOPparamError, EvolTrackError
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


def pastis_loglike(samples, params, input_dict, datadict):
    """
    A wrapper to run the PASTIS.MCMC.get_likelihood function.

    Computes the loglikelihood on a series of points given in samples using
    PASTIS.MCMC.get_likelihood.

    :param np.array samples: parameter samples on which to compute log
    likelihood. Array dimensions must be (np x k), where *np* is the number of
    samples and *k* is the number of model parameters.

    :param list params: parameter names. Must be in the PASTIS format: \
    objectname_parametername.

    :return:
    """
    # Prepare output arrays
    loglike = np.zeros(samples.shape[0])

    for s in range(samples.shape[0]):

        for parameter_index, full_param_name in enumerate(params):

            # Modify input_dict
            obj_name, param_name = full_param_name.split('_')
            input_dict[obj_name][param_name][0] = samples[s, parameter_index]

        # Construct chain state
        chain_state, labeldict = MCMC.tools.state_constructor(input_dict)

        try:
            # Compute likelihood for this state
            ll, loglike[s], likeout = MCMC.PASTIS_MCMC.get_likelihood(
                chain_state, input_dict, datadict, labeldict, False,
                False)

        except (ValueError, RuntimeError, EvolTrackError, EBOPparamError):
            print('Error in likelihood computation, setting lnlike to -np.inf')
            loglike[s] = -np.inf
            pass

    return loglike


def pastis_logprior(samples, params, input_dict, pdfdict=None):
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

    if pdfdict is None:
        # Build prior instances.
        pdfdict = MCMC.priors.prior_constructor(input_dict, {})

    for s in range(samples.shape[0]):

        for parameter_index, full_param_name in enumerate(params):

            # Modify input_dict
            obj_name, param_name = full_param_name.split('_')[:2]
            input_dict[obj_name][param_name][0] = samples[s, parameter_index]

        # Construct chain state
        chain_state, labeldict = MCMC.tools.state_constructor(input_dict)

        # Compute prior distribution for this state
        prior_probability = MCMC.priors.compute_priors(pdfdict, labeldict)[0]
        logprior[s] = np.log(prior_probability)

    return logprior


def lnpost(samples, params, input_dict, datadict, priordict=None):

    x = samples.copy()
    x = np.reshape(x, (1, len(x)))

    lnprior = pastis_logprior(x, params, input_dict, priordict)
    # Do not compute likelihood for points outside prior.
    if lnprior[0] == -np.inf:
        return -np.inf
    lnlike = pastis_loglike(x, params, input_dict, datadict)

    return (lnprior + lnlike)[0]
