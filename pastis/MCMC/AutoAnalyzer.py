"""
Module for automatic batch analysis of chains
"""

import os
import shutil
import analysis
import datetime
import numpy as np
from math import exp
from .. import ModelComparison
from .. import resultpath, configpath


def analyze_chains(target, hyp1, hyp2=None, **kwargs):
    """
    Analyse all chains of target. Compares hypotheses hyp1 and hyp2.

    Parameters
    ----------
    target: string.
        Name of target. This parameter determines the folder where the .mcmc
        files of the hypotheses are located.

    Other parameters
    ----------------
    solutionfile: boolean.
        Decides if a solution file is produced for each beta.

    suffix: str.
        Suffix to append to hypotheses names in order to find pastis files.
        (only used if solution file is to be created).
        
    BIparam: string.
        Name of parameter used to compute BI of chains. Default: posterior

    BIpad: float.
        Fraction of chain to discard after BI to insure convergence of all
        parameters. Default: 0.25
        
    CLstep: int.
        Step size used for the computation of the correlation length.
        If negative, this value will be used for all chains, and no computation
        will be done.
        Default: 10
        
    CLbi: float.
        BI used for the CL computation. Usually a large value to increase
        speed.
        Default: 0.8

    Nmin: int.
        Minimum number of independent points per chain to use for merging.
        Default: 100
        
    MCmethod: string.
        Method to pass to compare_models function. Default: 'beta'

    MCintmethod: string.
        Method to pass to compare_models function for integration.
        Default: 'simple'

    beta: float.
        If given, only chains with this beta will be analyzed.

    """

    solutionfile = kwargs.pop('solutionfile', True)
    suffix = kwargs.pop('suffix', None)
    biparam = kwargs.pop('BIparam', 'posterior')
    bipad = kwargs.pop('BIpad', 0.25)
    clstep = kwargs.pop('CLstep', 10)
    clbi = kwargs.pop('CLbi', 0.8)
    nmin = kwargs.pop('Nmin', 100)
    mcmethod = kwargs.pop('MCmethod', 'beta')
    mcintmethod = kwargs.pop('MCintmethod', 'simple')
    singlebeta = kwargs.pop('beta', None)

    # Check if needed files exist
    if solutionfile:
        if hyp2 is not None:
            hyps = (hyp1, hyp2)
        else:
            hyps = (hyp1,)

        pfiles = []
        for hyp in hyps:
            if suffix is not None:
                pasfile = target + '_' + hyp + '_' + suffix + '.pastis'
            else:
                pasfile = target + '_' + hyp + '.pastis'

            pastisfile = os.path.join(configpath, target, pasfile)

            if not os.path.exists(pastisfile):
                pastisfile = pastisfile.replace(target + '_', '')

            if not os.path.exists(pastisfile):
                raise Exception('None {} and {} configuration files found!'
                                ''.format(pastisfile,
                                          pastisfile.replace(hyp1,
                                                             target + '_' + hyp)
                                          )
                                )

            pfiles.append(pastisfile)
        pastisfile = pfiles[0]
        if hyp2 is not None:
            pastisfile2 = pfiles[1]

    # Open and initilize log file
    if hyp2 is not None:
        logname = target + '_' + hyp1 + '_' + hyp2 + '.log'
    else:
        logname = target + '_' + hyp1 + '.log'

    now = datetime.datetime.isoformat(datetime.datetime.now())
    logfilepath = os.path.join(resultpath, target, logname)
    if os.path.exists(logfilepath):
        logfile = open(logfilepath, 'a+')
        logfile.write('=' * 80 + '\n')
        logfile.write('=' * 80 + '\n')
    else:
        logfile = open(logfilepath, 'w')
        logfile.write('Log file for {0}; hypothesis {1}, '
                      '{2}\n'.format(target, hyp1, hyp2))
    logfile.write('Analysis start time: {}\n'.format(now))

    if hyp2 is None:
        logfile.write('Only one hypothesis given; model comparison will not '
                      'be performed.\n')
        comparemodels = False
    else:
        logfile.write('Comparing models:\n{0}\n{1}\n'.format(hyp1, hyp2))
        comparemodels = True
    ##

    ###
    # ## Read chains
    ###
    fnames1, vds1 = analysis.read_chains(target, hyp1, beta=singlebeta)

    # Log
    logfile.write('Chains read for hypothesis1\n')
    for ff in fnames1:
        logfile.write(ff + '\n')

    if comparemodels:
        fnames2, vds2 = analysis.read_chains(target, hyp2, beta=singlebeta)

        # Log
        logfile.write('Chains read for hypothesis2\n')
        for ff in fnames1:
            logfile.write(ff + '\n')
    ###

    #
    # Get list of betas
    #
    betas1 = []

    for vd in vds1:
        betas1.append(vd.beta)
    betas1 = np.array(betas1)

    # Log
    logfile.write('Beta values for hypothesis1: ')
    for bb in np.unique(betas1):
        logfile.write('{:.4e} '.format(bb))
    logfile.write('\n')

    if comparemodels:
        betas2 = []
        for vd in vds2:
            betas2.append(vd.beta)
        betas2 = np.array(betas2)

        # Log
        logfile.write('Beta values for hypothesis2: ')
        for bb in np.unique(betas2):
            logfile.write('{:.4e}  '.format(bb))
        logfile.write('\n')

    logfile.write('=' * 80 + '\n')
    ###

    ###
    # Iterate over all betas
    ###
    for bb in np.unique(betas1):

        print('Beta = %.2e' % bb)
        # Analyse chains for each hypothesis and register result
        logfile.write('Hypothesis {0}\n'.format(hyp1))

        vdm = single_analysis(vds1, fnames1, betas1, bb, biparam, bipad, clstep,
                              clbi, nmin, logfile)
        ###
        # Make solution file
        ###
        if solutionfile:
            solfile = analysis.make_solution_file(vdm, pastisfile=pastisfile,
                                                  outputname=True)
            solfdest = solfile.replace('.solution', '_beta{0:.3e}.'
                                                    'solution'.format(bb))
            shutil.move(solfile, solfdest)

            # log
            logfile.write('Solution file saved to {0}\n'.format(solfdest))

        #
        # Register chain
        #
        analysis.register_chain(vdm, beta=str(bb), target=target, runid=hyp1,
                                force=True)
        # log
        logfile.write('Merged chain registered.\n')
        logfile.write('=' * 80 + '\n')
        print('=' * 80)

    #
    # Repeat for hypothesis 2
    #
    if comparemodels:
        for bb in np.unique(betas2):
            print('Beta = %.2e' % bb)
            vdm = single_analysis(vds2, fnames2, betas2, bb, biparam, bipad,
                                  clstep, clbi, nmin, logfile)

            ###
            # Make solution file
            ###
            if solutionfile:
                solfile = analysis.make_solution_file(vdm,
                                                      pastisfile=pastisfile2,
                                                      outputname=True)

                solfdest = solfile.replace('.solution', '_beta{0:.3e}.'
                                                        'solution'.format(bb))
                shutil.move(solfile, solfdest)

                # log
                logfile.write('Solution file saved to %s\n' % solfdest)

            ###
            # Register chain
            ###
            analysis.register_chain(vdm, beta=str(bb), target=target,
                                    runid=hyp2, force=True)
            # log
            logfile.write('Merged chain registered.\n')
            logfile.write('=' * 80 + '\n')
            print('=' * 80)

        ####
        # MODEL COMPARISON
        ####
        bbs1, bbs2, ml1, ml2, evid = \
            ModelComparison.compare_models(target, hyp1, hyp2,
                                           intmethod=mcintmethod,
                                           method=mcmethod,
                                           plot=False,
                                           minbeta=1e-10)

        # log
        ids = [hyp1, hyp2]

        if mcmethod == 'beta':
            try:
                evidence_ratio = exp(max(evid) - min(evid))
            except OverflowError:
                logfile.write(
                    '## Warning! Overflow Error, evidence ratio is the '
                    'exponent of e.')
                evidence_ratio = (max(evid) - min(evid))

        else:
            evidence_ratio = max(evid) / min(evid)

        logfile.write(
            'Evidence ratio M[{0}]/M[{1}] = {2:.3e}\n'.format(
                ids[np.argmax(evid)],
                ids[np.argmin(evid)],
                evidence_ratio
            ))
    ##
    now = datetime.datetime.isoformat(datetime.datetime.now())
    logfile.write('Analysis finish time: %s\n' % now)

    logfile.close()
    return


def single_analysis(vds, fnames, betas, beta, burninparam, burninpad,
                    corrlength_step, corrlength_burnin, nmin, logfile):
    ###
    # Analysis for each beta
    ###
    #vds = np.atleast_1d(vds)
    #fnames = np.atleast_1d(fnames)
    #fnames = np.atleast_1d(fnames)

    vds = np.compress(betas == beta, vds)
    fnames = np.compress(betas == beta, fnames)

    ###
    # Reject chains with problems
    ###
    condnan = []
    for i in range(len(vds)):
        ll = vds[i].get_value_dict()[burninparam]
        condnan.append(np.any(np.isnan(ll)))
    condnan = np.array(condnan)

    fnames_mal = np.compress(condnan, fnames)

    vds = np.compress(-condnan, vds)
    fnames = np.compress(-condnan, fnames)

    if len(fnames_mal) > 0:
        # Log
        logfile.write('Chains rejected due to NaNs:\n')
        for ff in fnames_mal:
            logfile.write(ff + '\n')
        logfile.write('--\n')

    # Log
    logfile.write('{0} chains for beta {1:.4e}\n'.format(len(vds), beta))

    # Find BI (variance is not checked to avoid problems with low beta chains.
    z, zz, bi = analysis.find_BI(vds, param=burninparam, verbose=False,
                                 checkvariance=True)
    bi = np.array(bi)
    BI = bi + (1 - bi) * burninpad

    # log
    logfile.write('Burn-In estimated using {0}.\n'.format(burninparam))

    if corrlength_step > 0:
        # Compute CL
        corrdict = analysis.corrlength_multichain(vds, step=corrlength_step,
                                                  BI=corrlength_burnin,
                                                  plot=False,
                                                  verbose=False)
        # Log
        logfile.write('Correlation length computed with parameters: step {0}; '
                      'BI {1:.3f}\n'.format(corrlength_step, corrlength_burnin))

        # Analysis of CLs
        kk = np.array(corrdict.keys())
        vv = np.array(corrdict.values())

        for i in range(len(vds)):
            # Get the parameter with largest value for each chain
            maxparam = kk[np.argmax(vv[:, i])]
            maxvalue = np.max(vv[:, i])

            # Check if CL needs to be rerun, if the chain used is too short
            if maxvalue >= len(vds[i].get_value_dict()[maxparam]) * (1 - BI[i]):
                rerun = True
            else:
                rerun = False

            # log CL and BI
            logfile.write('Chain {0}\n'.format(fnames[i]))
            logfile.write('burn in = {0:.2f}.\n'.format(BI[i]))
            logfile.write(
                'Maximum CL for parameter {0}: {1}\n'.format(maxparam,
                                                             maxvalue))

            if rerun:
                logfile.write('## WARNING! CL for chain {0} should be computed '
                              'using a smaller BI value!\n'.format(fnames[i]))
                print('WARNING! CL for chain {0} should be computed using a '
                      'smaller BI value!'.format(fnames[i]))

        CL = analysis.corrlenchain(corrdict)

    else:
        CL = np.full(len(vds), abs(corrlength_step))
    ##

    ###
    # Select best chains
    ###
    # if beta > 1e-4:
    try:
        vdsa, fnamesa, vdsr, fnamesr, iaccept = \
            analysis.select_best_chains(vds, CL=CL, BI=BI, nmin=nmin,
                                        fnames=fnames, param=burninparam,
                                        )

    except:
        logfile.write('## ERROR! Problem selecting best chains!\n')
        logfile.close()
        raise Exception

    # Log
    logfile.write('{0} chains are rejected.\n'.format(len(vds) - len(vdsa)))

    if len(vds) - len(vdsa) > 0:
        logfile.write('Rejected chains:\n')
        for fn in fnamesr:
            logfile.write(fn + '\n')
        logfile.write('--\n')

    if len(vdsa) < int(0.5 * len(vds)) and beta > 1e-4:
        print('WARNING! More than half the chains have been rejected.')
        # log
        logfile.write(
            '## WARNING! More than half the chains have been rejected.\n')

    if len(vdsa) < 2 and beta > 1e-4:
        print('WARNING! Only best chain chosen!')
        # log
        logfile.write('## WARNING! Only best chain was chosen!\n')

    ###
    # Merge chains
    ###
    vdm = analysis.merge_chains(vdsa, BI=BI[iaccept], CL=CL[iaccept],
                                beta=beta)

    # log
    anyparam = np.random.choice(vdm.get_value_dict().keys())
    logfile.write('Merged chain for beta {0:.4e} has '
                  '{1} elements.\n'.format(beta,
                                           len(vdm.get_value_dict()[anyparam])))

    return vdm
