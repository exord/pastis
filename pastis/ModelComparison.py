import os
import numpy as n
import pickle

import scipy.interpolate
import scipy.integrate

from math import exp, log, log10, e
from . import resultpath
from . import PASTISroot
resultdir = PASTISroot+'/resultfiles/'

# To correctly handle TeX
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

import pylab as p



def loadchains(listfile, target):
    """
    Load all the chains in a given listfile (constructed with
    MCMC.analysis.register_chain)
    """
    listfile = os.path.split(listfile)[-1]
    fin = open(os.path.join(resultpath, target, listfile), 'r')
    chainlist = fin.readlines()
    fin.close()

    betas = []
    vds = []

    for ll in chainlist:

        if ll.startswith('#'): continue
        
	betas.append(float(ll.split()[0]))

	chainfile = os.path.join(resultpath, target, ll.split()[1])
	ff = open(chainfile, 'r')
	vds.append(pickle.load(ff))
	ff.close()

    return n.array(betas), vds

	
def compute_logL(vds, median = False):
    """
    Compute the mean logL for a list of chains, usually chains with
    different tempering parameter beta. The chains have to be thinned and
    the burn-in period has to be already clipped. 

    Parameters
    ----------
    vds: iterable
        An iterable object with the VDchain instances of different chains.
    """
    meanLogL = []
    for i in range(len(vds)):
        logL = vds[i].get_value_dict()['logL']
        if median:
            meanLogL.append(n.median(logL))
        else:
            meanLogL.append(n.mean(logL))

    return n.array(meanLogL)
			

def get_evidence(C1, C2, betas1, betas2 = None, method = 'beta',
                 intmethod = 'simple', **kwargs):
    """
    Get the odd ratio for two list of chains having different values of the
    tempering parameter beta.
    Each list can be evaluated at different betas.
    """

    TPMlag = kwargs.pop('TPMlag', 5e3)
    h = TPMlag
    TPMlambda = kwargs.pop('TPMlambda', 1e-2)
    use_median = kwargs.pop('use_median', False)
    BI = kwargs.pop('BI', 0.0)
    
    betas1 = n.array(betas1)
    if betas2 == None:
	betas2 = betas1.copy()
    else:
	betas2 = n.array(betas2)
	
    if method == 'beta':
        ml1 = compute_logL(C1, median = use_median)
        ml2 = compute_logL(C2, median = use_median)

        ml1 = ml1[n.argsort(betas1)]
        ml2 = ml2[n.argsort(betas2)]
        betas1 = n.sort(betas1)
        betas2 = n.sort(betas2)

        #Integrate for all betas
        logE1 = integrate_logL(betas1, ml1, intmethod)
        logE2 = integrate_logL(betas2, ml2, intmethod)

        return betas1, betas2, ml1, ml2, logE1, logE2
    
    elif method == 'HM' or method == 'TPM' or method == 'M':

        methoddict = {'M' : 'Median',
                      'HM' : 'Harmonic Mean',
                      'TPM' : 'Truncated Posterior-Mixture'
                      }
                      
        
        if (not 1 in betas1) or (not 1 in betas2):
            print('Error! Need chains with beta = 1 for %s estimation.'%\
                      (methoddict[method])
                  )
            return

        vd1 = n.compress(betas1 == 1, C1)[0]
        vd2 = n.compress(betas2 == 1, C2)[0]

        l1 = vd1.get_value_dict()['logL']
        l2 = vd2.get_value_dict()['logL']

        l1 = l1[len(l1)*BI:]
        l2 = l2[len(l2)*BI:]
        
        # Mean estimation
        if method == 'M':
            logE1 = n.sum(l1)/len(l1)
            logE2 = n.sum(l2)/len(l2)

        # Harmonic mean estimation
        elif method == 'HM':
            n1 = len(l1); n2 = len(l2)
            logE1 = log(n1) + l1[0] - log(1 + n.sum(n.exp(-(l1[1:] - l1[0]))))
            logE2 = log(n2) + l2[0] - log(1 + n.sum(n.exp(-(l2[1:] - l2[0]))))

        # Truncated posterior-mixture
        elif method == 'TPM':
            try:
                lp1 = vd1.get_value_dict()['posterior']
                lp2 = vd2.get_value_dict()['posterior']


            except KeyError:
                raise Exception('Error! \'posterior\' key not available. Impossible to use TPM method.')

            lp1 = lp1[len(lp1)*BI:]
            lp2 = lp2[len(lp2)*BI:]

            ## Change base of lp1 and lp2 to e
            lp1 = lp1/log10(e)
            lp2 = lp2/log10(e)

            ## Get log(prior)
            lprior1 = lp1 - l1
            lprior2 = lp2 - l2
            
            l1lambda = log(1 - TPMlambda)
            llambda = log(TPMlambda)


            logE1 = computeTPM(lp1, lprior1, h, TPMlambda)
            logE2 = computeTPM(lp2, lprior2, h, TPMlambda)

            """
            ### CHECK THE USE OF SLICING AT h!
            rho1 = l1lambda + lp1[h:] + \
                n.log(1 + n.exp(llambda - l1lambda + lp1[:-h] - lp1[h:])
                      )

            rho2 = l1lambda + lp2[h:] + \
                n.log(1 + n.exp(llambda - l1lambda + lp2[:-h] - lp2[h:])
                      )
 
            ## Compute estimate
            logE1 = lp1[h] - lprior1[h] + \
                log(1 + n.sum(n.exp(lp1[h:] - lp1[h] + rho1[h] - rho1)
                              )
                    ) - \
                log(1 + n.sum(n.exp(lprior1[h:] - lprior1[h] + rho1[h] - rho1)
                              )
                    )

            logE2 = lp2[h] - lprior2[h] + \
                log(1 + n.sum(n.exp(lp2[h:] - lp2[h] + rho2[h] - rho2)
                              )
                    ) - \
                log(1 + n.sum(n.exp(lprior2[h:] - lprior2[h] + rho2[h] - rho2)
                              )
                    )
            """

        return logE1, logE2
        

def computeTPM(lp, lprior, TPMlag, TPMlambda):
    

    ## Define auxiliary variables
    h = TPMlag
    log_lambda = log( TPMlambda / (1.0 - TPMlambda) )
    
    gamma_i = lp[h] + log( 1 + exp( log_lambda + lp[0] - lp[h] ) ) - lp[h + 1:] - n.log( 1 + n.exp( log_lambda + lp[1: -h] - lp[h + 1:]) )

    logE = lp[h] - lprior[h] + log( 1 + n.sum( n.exp (lp[h + 1:] - lp[h] + gamma_i) ) ) \
        - log( 1 + n.sum( n.exp (lprior[h + 1 :] - lprior[h] + gamma_i ) ) )

    return logE


def integrate_logL(x, y, method = 'simple'):
    if method == 'simple':
        ymin = y.min()
        y = y + 2*abs(ymin)
	return n.sum((x[1:] - x[:-1])*(y[1:] + y[:-1])*0.5) - 2*abs(ymin)
    

    elif method == 'quad':
	# Interpolate function
	ff = scipy.interpolate.interp1d(x, y, kind = 'quadratic',
					bounds_error = False, fill_value = 0.0)
	return scipy.integrate.quad(ff, 0, 1)[0]

    elif method == 'gquad':
	# Interpolate function
	ff = scipy.interpolate.interp1d(x, y, kind = 'quadratic',
					bounds_error = False, fill_value = 0.0)
	return scipy.integrate.fixed_quad(ff, 0, 1)[0]

    else:
	raise Exception('Invalid integration method')
	

def compare_models(target, id1, id2, intmethod = 'simple', plot = True,
                   method = 'beta', use_median = False, minbeta = 0.0,
                   offset = False, **kwargs):
    """
    Compare models id1 and id2 for target.
    This requieres the files .lst to be present in the resultfiles directory
    of target for each id.
    """
    npar1 = kwargs.pop('npar1', None)
    npar2 = kwargs.pop('npar2', None)
    ndata = kwargs.pop('ndata', None)
    
    verbose = kwargs.pop('verbose', True)
    
    f1 = '%s_%s_mergedchain.lst'%(target, id1)
    f2 = '%s_%s_mergedchain.lst'%(target, id2)
    betas1, vds1 = loadchains(f1, target)
    betas2, vds2 = loadchains(f2, target)

    vds1 =  n.compress(betas1 >= minbeta, vds1)
    vds2 =  n.compress(betas2 >= minbeta, vds2)
    betas1 =  n.compress(betas1 >= minbeta, betas1)
    betas2 =  n.compress(betas2 >= minbeta, betas2)
    
    
    ids = [id1, id2]
    
    if method == 'beta':

        betas1, betas2, ml1, ml2, logE1, logE2 = get_evidence(vds1, vds2,
                                                              betas1, betas2,
                                                              method = method,
                                                              intmethod = intmethod,
                                                              use_median = use_median,
                                                              **kwargs
                                                              
                                                              )
        
    elif method == 'BIC':

        if ndata == None:
            print('Provide ndata to run BIC computation')
            return
            
        if npar1 == None:
            npar1 = len(vds1[0]._value_dict) - 2

        BIC1 = get_BIC(vds1[0], npar = npar1, ndata = ndata, **kwargs)

        if npar2 == None:
            npar2 = len(vds2[0]._value_dict) - 2

        BIC2 = get_BIC(vds2[0], npar = npar2, ndata = ndata, **kwargs)


        print('BIC_{0}: {1}; BIC_{2}: {3}'.format(id1, BIC1, id2, BIC2))
        print('-DeltaBIC/2.0 = -(BIC_{1} - BIC_{0}) / 2 ~ ln(E_{1}/E_{0}) = {2}'.format(id1, id2, (BIC2 - BIC1)*-0.5))
        return

    else:
        logE1, logE2 = get_evidence(vds1, vds2, betas1, betas2, method = method,
                                    **kwargs)

        
              

    # Compute evidence ratio
    evid = [logE1, logE2]
    try:
        evidence_ratio = exp(max(evid) - min(evid))
    except OverflowError:
        print('Warning! Overflow Error, evidence is the exponent of e.')
        evidence_ratio = (max(evid) - min(evid))

    # Print
    evidencestr = 'Evidence ratio M[%s]/M[%s] = %.3e'%(ids[n.argmax(evid)],
                                                       ids[n.argmin(evid)],
                                                       evidence_ratio
                                                       )
    evidencestr = evidencestr.replace('_', ' ')

    if verbose: print(evidencestr)

    """
    if plot and method == 'beta':
	f1 = p.figure()
	ax = f1.add_subplot(111)

	ax.plot(betas1, ml1, 'o-', label = id1)
	ax.plot(betas2, ml2, 'o-', label = id2)
	ax.set_ylabel('<ln(L)>', fontsize = 16)
	ax.set_xlabel('Tempering parameter beta', fontsize = 16)
	ax.legend(loc = 0)
	p.show()
    """
    
    if plot and method == 'beta':
	f1 = p.figure()
	ax = f1.add_subplot(111)

        if offset:
            oset = abs( 2*min(ml1 - ml2) )
            y = ml1 - ml2 + oset
            ax.loglog(betas1, y, 'o-k', mfc = 'None', ms = 8)
        else:
            y = ml1 - ml2
            oset = 0.0
            ax.semilogx(betas1, y, 'o-k', mfc = 'None', ms = 8)
            
	ax.set_ylabel('$<\ln(\mathcal{L}_1)> - <\ln(\mathcal{L}_2)>$', fontsize = 16)
	ax.set_xlabel('Tempering parameter beta', fontsize = 16)

        estr = 'Bayes\' Factor: %.2e'%evidence_ratio
        ax.text(0.5, 0.7, estr, transform = ax.transAxes, ha = 'left',
                size = 16)

        # Shadow area under curve and horizontal line
        ## With offset
        ## Without offset
        ax.fill_between(betas1, y, y2 = oset, color = '0.55')
        ax.axhline(offset, ls = ':', color = 'k')
        
        
	#p.show()

    if method == 'beta':
        return betas1, betas2, ml1, ml2, evid, f1

    else:
        return evid[0], evid[1]
        #return (evid[0] - evid[1])


def get_BIC(vd, npar, ndata, **kwargs):

    if isinstance(vd, dict):
        logLmax = vd['logL'].max()
    else:
        logLmax = vd._value_dict['logL'].max()

    return -2 * logLmax + npar * log(ndata)
