# 31-07-2012 Corrected bug in get_best_values() method of Chain that was making
#            it choose the wrong values.
# 22-05-2012 Implementation of Binormal and PowerLaw prior types. Also
#            included "sin" type, but should make sure that x is in radians. 
# 18-04-2012 Modification of numerous methods and inclusion of new ones
#            to incorporate Principal Component Analysis
# 26-03-2012 Inclination limits changed. Inclination of triple 
#            are not afected by the posible limitation of JKTEBOP.

# CLASSES used for MCMC
##
import sys
import numpy as n
import scipy
from math import log10, log


# PARAMETER
class Parameter(object):
    def __init__(self, value, object, jump=True, label=None, family=[],
                 proposal_function='gaussian', proposal_scale=None):

        self._value = value
        self._formervalue = value
        self._object = object  # Object it modifies

        self.jump = jump
        self.label = label  # A name to identify object
        self.family = family # Potential list of parameter dependeces

        # If jump parameter, at initialization, set some mandatory (and
        # dynamic!) attributes
        if jump:
            # Number of steps since this the jump scale of this paremeter
            # was last updated.
            self.Nstep_since_update = 0.0

            # Number of steps accepted since last update of proposal scale
            self.Naccepted_since_update = 0.0

            # Number of accepted steps
            self.Naccepted = 0.0

            # Number of proposals changes made for this parameter
            self.Nstep = 0.0

            # Set proposal function
            # It can be the string 'gaussian' or a function object
            self.propfunc = proposal_function

            # Set proposal scale
            if proposal_scale is None:
                self.propscale = abs(self.get_value())*0.1
            else:
                self.propscale = proposal_scale

            # Set frequency of scale update (see Ford (2006))
            self.update_interval = 2.0
            
            # Set value of last scale change (initially, None)
            self.lastk = None

            # Prepare list to hold jump history
            self.jumphistory = []

    def update_jump_history(self):
        self.jumphistory.append((self.Nstep, self.propscale))
            
    def get_object(self):
        return self._object

    def get_value(self):
        return self._value

    def get_formervalue(self):
        return self._formervalue

    def set_value(self, value):
        try:
            value = float(value)
            self._formervalue = self.get_value()
            self._value = value
        except TypeError:
            raise Exception('Value cannot be converted to float.')
        except ValueError:
            print('String {} cannot be converted to float. It cannot be '
                  'the argument of set_value()'.format(value))

    def revert_to_old_value(self):
        self.set_value(self.get_formervalue())
        
    def copy(self):
        """
        Produce a copy of the paremeter
        """
        pp = Parameter(self.get_value(), self.get_object())

        for attr in dir(self):
            # Avoid copying methods
            if type(getattr(self, attr)) == types.MethodType:
                continue
            
            # Avoid everything starting with '__'
            # (JUST TO BE SAFE, I DON'T UNDERSTAND THIS YET)
            if attr.startswith('__'):
                continue
            setattr(pp, attr, getattr(self, attr))

        return pp


# CHAIN
class Chain(object):
    """
    A Class of Markov Chains.

    Methods:
    MAKE DOC FOR METHODS
    """

    def __init__(self, N, theta0, Npca = n.inf, usePCA = True):

        self._currentstate = theta0 #Initial state
        self.N = N  #Number of steps in chain
        self.Naccept = 0 # Number of accepted steps

        self.usePCA = usePCA
        self.Npca = Npca # Number of steps after which to initiate PCA

        
        self.values = n.zeros((self.N, len(theta0)), dtype = float)
        self.append(0)	
        self._likelihood = n.zeros((self.N, 3), dtype = object)
        self._posterior = n.zeros(self.N)
        
        # Create labeldictionary
        self._labeldict = dict((theta0[i].label, theta0[i])
                               for i in range(len(theta0)))
        
        # Initialise jumpscale dictionary
        self._jumpscaledict = {}
        for par in theta0:
            if par.jump: par.update_jump_history()


        # Get jump parameters indices
        jumpind = []
        for jj, par in enumerate(theta0):
            if par.jump: jumpind.append(jj)
        self.jumpind = jumpind

    	    ### NEW: FOR MA
        self._meanX = n.zeros((self.N, len(jumpind)), dtype = float)
        self._meanX[0] = self.get_current_state_jumping_values()
        self._S0 = n.diag(self.get_current_state_jumping_values()*0.1)

        xx = self.get_current_state_jumping_values()
        for i in range(len(xx)):
            for j in range(len(xx)):
                if j > i:
                    self._S0[i,j] = self._S0[j,i] = 0.1
		    
	   ###
	
        # For PCA
        self._currentstatePCA = None
        if Npca < N:
            njump = len(jumpind)
            # Initialice PCAvalues array
            self.PCAvalues = n.zeros((N - Npca, njump), dtype = float)
            self.M = n.identity(njump) # Change-of-base matrix (identity)

            
    def append(self, i):

        # Get values of X and add them to values list
        vals = self.get_current_state_values()
        vals = vals.reshape((1, len(vals)))
        self.values[i] = vals
        
        # For PCA
        if i >= self.Npca and self.usePCA:
            # Also update PCAvalues array
            vals = self.get_current_pca_values()
            vals = vals.reshape((1, len(vals)))
            self.PCAvalues[i - self.Npca] = vals
        
    def accept_proposal(self, step, prior, L, logL, Ldict):
        self.Naccept += 1
        self.append(step)
        self._likelihood[step] = [L, logL, Ldict]
        self._posterior[step] = logL / log(10) + log10(prior)
        return

    def reject_proposal(self, step, prior, L, logL, Ldict):
        # Recover former value of all parameters
        for par in self._currentstate:
            par.revert_to_old_value()

        # Recover former value of all Principal components
        if step > self.Npca and self.usePCA:
            for par in self._currentstatePCA:
                par.revert_to_old_value()

        # Append this state (the former one) to the chain
        self.append(step)
        self._likelihood[step] = [L, logL, Ldict]
        self._posterior[step] = logL / log(10) + log10(prior)
        
        return
    
    def get_current_state(self, PCA = False):
        if PCA:
            return self._currentstatePCA
        else:
            return self._currentstate

    def get_state(self, i):
        """
        Return state i of chain, except for proposal scale
        """
        S = []
        X = self._currentstate
        vals = self.values[i, :]

        for j, par in enumerate(X):
            spar = Parameter(vals[j], '', jump = par.jump, label = par.label)
            S.append(spar)
        return S

    def get_current_state_values(self):
        values = []
        for par in self._currentstate:
            values.append(par.get_value())
        #return n.resize(n.array(values), (1, len(values)))
        return n.array(values)

    def get_current_state_jumping_values(self):
        values = []
        for par in self._currentstate:
            if par.jump: values.append(par.get_value())
        #return n.resize(n.array(values), (1, len(values)))
        return n.array(values)
    
    def get_current_pca_values(self):
        values = []
        for par in self._currentstatePCA:
            values.append(par.get_value())
        return n.array(values)
            
        
    def get_value_dict(self):
        X = self.get_current_state()
        valuedict = {}
        for i in range(len(X)):
            if X[i].jump:
                valuedict[X[i].label] = self.values[:, i]

        valuedict['logL'] = self.get_logL()
        valuedict['posterior'] = self.get_posterior()
        
        return valuedict

    def get_median_values(self, BI):
        """
        Return dictionary with median values of trace
        """
        medianvalues = {}
        vd = self.get_value_dict()
        for key in self._labeldict.keys():
            skey = key.split('_')
            if skey[0] in medianvalues:
                pass
            else:
                medianvalues[skey[0]] = {}
                
            if key in vd.keys():
                bvi = n.median(vd[key][round(BI*self.N):])
                medianvalues[skey[0]][skey[1]] = [bvi, ]
            else:
                bvi = self._labeldict[key].get_value()
                medianvalues[skey[0]][skey[1]] = [bvi, ]

        return medianvalues
            
    def get_best_values(self, BI = 0.0):
        """
        Return dictionary with median values of trace
        """
        bestvalues = {}
        vd = self.get_value_dict()
        logL = self.get_logL()
        indi = n.round(BI*self.N)
        for key in self._labeldict.keys():
            skey = key.split('_')
            if skey[0] in bestvalues:
                pass
            else:
                bestvalues[skey[0]] = {}
                
            if key in vd.keys():
                indmax = n.argmax(logL[indi:])
                bvi = vd[key][indi:][indmax]
                bestvalues[skey[0]][skey[1]] = [bvi, ]
            else:
                bvi = self._labeldict[key].get_value()
                bestvalues[skey[0]][skey[1]] = [bvi, ]

        return bestvalues

    def confidence_intervals(self, q = 0.6827, nbins = 50, BI = 0.2):
        """
        Computes  median and sigmas for the parameters in chain, using the
        q*100% confidence limits, and nbins.
        """
        vd = self.get_value_dict()
        istart = n.round(BI*self.N)
        
        # Minimum and maximum prob
        qmin = (1 - q)*0.5
        qmax = (1 + q)*0.5

        for param in vd.keys():
            # Compute histogram, and cumulative distribution
            m, bins = n.histogram(vd[param][istart:], nbins, cumulative = True)

            x = bins[:-1] + 0.5*n.diff(bins)
            ci = m.astype(float).cumsum()/n.sum(m)

            imin1 = float(n.argwhere(n.less(qmin)).max())
            imin2 = imin1 + 1.0
            funmin = scipy.interp([ci[imin1], ci[imin2]], [x[imin1], x[imin2]])
            lower_limit = funmin(qmin)

            imax1 = float(n.argwhere(n.less(qmax)).max())
            imax2 = imax1 + 1.0
            funmax = scipy.interp([ci[imax1], ci[imax2]], [x[imax1], x[imax2]])
            upper_limit = funmin(qmax)

            medianvalue = n.median(vd[param][istart:])
            
            hc = upper_limit - medianvalue
            lc = medianvalue - lower_limit
            # print('%s: %.6f + %.1e - %.1e'%(param, medianvalue, hc[i], lc[i]))
        return 'CACA'
        
    def get_trace(self, paramname):
        vd = self.get_value_dict()
        return vd[paramname]

    def get_jumpscale_dict(self):
        X = self.get_current_state()
        jumpdict = {}
        for par in X:
            if par.jump:
                jumpdict[par.label] = par.jumphistory
        return jumpdict
        
    def get_jumpscale_history(self, paramname):
        """
        Return the values that the jump scale of a given parameter has
        taken, and the chain step in which it was taken.
        """
        jumpdict = self.get_jumpscale_dict()
        return jumpdict[paramname]

    def get_L(self):
        L = []
        for Llist in self._likelihood:
            L.append(Llist[0])
        return n.array(L)
    
    def get_logL(self):
        logL = []
        for Llist in self._likelihood:
            logL.append(Llist[1])
        return n.array(logL)
        
    def get_Ldict(self):
        Ldict = {}
        for i, Llist in enumerate(self._likelihood):

            for key in Llist[2]:
                if i == 0:
                    Ldict[key] = [Llist[2][key],]
                else:
                    Ldict[key].append(Llist[2][key])

        for key in Ldict:
            Ldict[key] = n.array(Ldict[key])
        
        return Ldict          

    def get_posterior(self):
        return self._posterior
    
    def savetofile(self, file):
        import pickle
        f = open(file, 'w')
        pickle.dump(self, f)
        f.close()
        return

    def savetoasciifile(self, file):
        """
        Save values of chains to an ascii file
        """
        kkeys = self._labeldict.keys()
        vd = self.get_value_dict()
        
        f = open(file, 'w')
        for i, key in enumerate(kkeys):
            if i != (len(kkeys) - 1): f.write(key+'\t')
            else: f.write(key+'\n')

        for i, key in enumerate(kkeys):
            if i != (len(kkeys) - 1): f.write('-'*len(key)+'\t')
            else: f.write('-'*len(key)+'\n')

        for i in range(len(self.values)):
            for jj, j in enumerate(kkeys):
                if jj != (len(kkeys) -1): f.write('%.5f\t'%vd[j][i])
                else: f.write('%.5f\n'%vd[j][i])
        f.close()
       
        
    def plot_trace(self, paramname, **kwargs):
        """
        Plot the chain of a given parameter

        - kwargs are passed to plot method, except for:

        * xlabel & ylabel: labels used for x and y axis. This
        defauts to 'Chain step' and parmname, respectively.
        """
        xlabel = kwargs.pop('xlabel', 'Chain step')
        ylabel = kwargs.pop('ylabel', paramname)


        try:
            p = sys.modules['pylab']
        except KeyError:
            import pylab as p

        fig = p.figure()
        ax = fig.add_subplot(111)

        ax.plot(self.get_trace(paramname), **kwargs)

        ax.set_xlabel(xlabel, fontsize = 16)
        ax.set_ylabel(ylabel, fontsize = 16)
        fig.show()

    def plot_hist(self, paramname, **kwargs):
        """
        Plot histogram of a given parameter

        - kwargs are passed to hist method, except for:

        * xlabel & ylabel: labels used for x and y axis. This
        defauts to 'Chain step' and N or PDF, respectively.

        * title: title of figure. Defaults to paramname
        """
        xlabel = kwargs.pop('xlabel', 'Chain step')
        ylabel = kwargs.pop('ylabel', 'N')
        title = kwargs.pop('title', paramname)

        try:
            p = sys.modules['pylab']
        except KeyError:
            import pylab as p

        fig = p.figure()
        ax = fig.add_subplot(111)

        bb = ax.hist(self.get_trace(paramname), **kwargs)
        ax.set_xlabel(xlabel, fontsize = 16)

        if n.sum(bb[0]*(bb[1][1:] - bb[1][:-1])) == 1.0 and ylabel == 'N':
            ylabel = 'PDF'
            
        ax.set_ylabel(ylabel, fontsize = 16)
        ax.set_title(title)
        fig.show()

    def plot_correlation(self, paramname1, paramname2, **kwargs):
        """
        Plot the chain of a given parameter

        - kwargs are passed to plot method, except for:

        * xlabel & ylabel: labels used for x and y axis. This
        defauts to paramname1 and parmname2, respectively.

        Also, linestyle is fixed to '', so that point markers are used.
        """
        xlabel = kwargs.pop('xlabel', paramname1)
        ylabel = kwargs.pop('ylabel', paramname2)

        try:
            p = sys.modules['pylab']
        except KeyError:
            import pylab as p

        fig = p.figure()
        ax = fig.add_subplot(111)

        ax.plot(self.get_trace(paramname1), self.get_trace(paramname2),
                ls='', **kwargs)

        ax.set_xlabel(xlabel, fontsize = 16)
        ax.set_ylabel(ylabel, fontsize = 16)
        fig.show()

    def plot_pyramid(self):
        """
        Construct pyramid plot for a chain.
        """
        raise NotImplemented
        return
