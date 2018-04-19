from ..MCMC.analysis import VDchain
import numpy as n
import pylab as p

params = {
          'text.usetex': False,
          }
p.rcParams.update(params)

def plot_multitraces(vds, tracekey, sampling = 100, BI = 0.5, **kwargs):
    """
    Plots traces of a parameter for a series of chains.

    Parameters
    ----------
    vds, list
        A list of value dictionaries containing the parameter to be plotted or
	of VDchain instances.
	All dictionaries or chains in the list must contain the parameter
        to plot.

    tracekey, string
        The name of the paremeter to plotted, as specified by the keys of
	the value dictionnaries in vds.

    Other parameters
    ----------------
    sampling, int
        The thinning factor. Only one point every "sampling" ones will be
	plotted.

    BI, float 
        The fraction of the chain to be discarded from the beginning. This is
        used if the traces contain the burn-in period.

    offset, float
        If different from 0.0, the traces will be offset vertically. The offset
	is computed as the median standard deviation over all traces times the
	number specified in this parameter.

    show, boolean
        Decides whether the pylab.show() command is executed
	
    Notes
    -----
    The remaining keyword parameters are passed to the plot function.
    """
    show = kwargs.pop('show', True)
    offset = kwargs.pop('offset', 0.0)
    f1 = p.figure()
    ax = f1.add_subplot(111)

    print p.isinteractive()


    #compute offset
    scatters = n.zeros(len(vds))
    for i in range(len(vds)):

	if isinstance(vds[i], VDchain):
	    vdi = vds[i].get_value_dict()
	elif isinstance(vds[i], dict):
	    vdi = vds[i].copy()
	else:
	    # Assume its a VDchain from before reloading the module....
	    print('Warning! Instance not recognized. Assuming it is VDchain instance.')
	    vdi = vds[i].get_value_dict()


	if i == 0:
	    start_index = n.round(BI*len(vdi[vdi.keys()[0]]))

	scatters[i] = n.std(vdi[tracekey][start_index:])

    oset = offset*n.median(scatters)

    for i, vd in enumerate(vds):

	if isinstance(vd, VDchain):
	    vdi = vd.get_value_dict()
	elif isinstance(vd, dict):
	    vdi = vd.copy()
        else:
            vdi = vd.get_value_dict()

        ax.plot(vdi[tracekey][start_index::sampling] - i*oset, label = str(i),
                **kwargs)
        ax.set_ylabel(tracekey, fontsize = 16)
        ax.set_xlabel('Chain step %% %d'%sampling, fontsize = 16)
        ax.yaxis.set_major_formatter(p.ScalarFormatter(useOffset = False))
        ax.legend(loc = 0, numpoints = 1, ncol = len(vds))

    if show:
	p.show()
    return ax


def plot_correlations(vds, key1, key2, sampling = 100, BI = 0.2, **kwargs):
    """
    Plot the correlation diagram between to parameters, for a series of
    chains.

    Parameters
    ----------
    vds, list
        The value dictionaries containing the parameter to be plotted.
	All dictionaries in the list must contain both parameters to plot.

    key1, key2, strings
        The name of the paremeter to plotted, as specified by the keys of
	the value dictionnaries in vds.
    
    Other parameters
    ----------------
    sampling, int
        The thinning factor. Only one point every "sampling" ones will be
	plotted.

    BI, float 
        The fraction of the chain to be discarded from the beginning. This is
        used if the traces contain the burn-in period.

    Notes
    -----
    The remaining keyword parameters are passed to the plot function.
    """
 
    f1 = p.figure()
    ax = f1.add_subplot(111)

    if isinstance(vds[0], VDchain):
	vdi = vds[0].get_value_dict()
    elif isinstance(vds[0], dict):
	vdi = vds[0].copy()

    start_index = n.round(BI*len(vdi[vdi.keys()[0]]))
    for i, vd in enumerate(vds):

	if isinstance(vd, VDchain):
	    vdi = vd.get_value_dict()
	elif isinstance(vd, dict):
	    vdi = vd.copy()

        ax.plot(vdi[key1][start_index::sampling],
                vdi[key2][start_index::sampling], label = str(i), **kwargs)
        ax.set_xlabel(key1); ax.set_ylabel(key2)
        ax.legend(loc = 0, numpoints = 1)

    p.show()
    return



def parallel_chains(C, BI = 0.2, sampling = 100, **kwargs):
    """
    Plot all traces of a given chain or value_dict C
    """
    
    plotL = False
    try:
        vd = C.get_value_dict()
        indi = n.round(C.N*BI)
        logL = C.get_logL()[indi:]
        plotL = True
    except AttributeError:
        vd = C.copy()
        indi = n.round(len(vd[vd.keys()[0]])*BI)
        try:
            logL = vd['logL']
            plotL = True
        except KeyError:
            pass
        
    nplots = len(vd) + 1.0

    f1 = p.figure()
    f1.subplots_adjust(hspace =0.0)

    print nplots
    for j, key in enumerate(n.sort(vd.keys())):
        ax = f1.add_subplot(nplots, 1, j + 1)
        ax.plot(vd[key][indi::sampling], color = 'k', **kwargs)
        ax.set_ylabel(key, rotation = 'horizontal')

        # Set major tick positions
        yl = ax.get_ylim(); dyl = yl[1]-yl[0]
        yts = n.arange(yl[0], yl[1], dyl/3.0)
        ax.yaxis.set_major_locator(p.FixedLocator(yts[1:]))

        # Trick to get nice labels
        s = '%e'%dyl
        a = s[s.find('e')+1:] # Contains power of range of axis in scientific notation
        if a[0] == '-':
            ndecimal = int(a[1:]) + 1
        else:
            ndecimal = 1
                
        ax.yaxis.set_major_formatter(p.FixedFormatter(n.round(yts[1:],
                                                              ndecimal)
                                                      )
                                     )
        for ytl in ax.get_yticklabels():
            ytl.set_fontsize(11)
                    
        ax.xaxis.set_major_formatter(p.NullFormatter())
        """
        #if j != len(vd.keys()) -1:
	#for ytl in ax.get_yticklabels():
	ax.yaxis.set_major_formatter(p.NullFormatter())
        """
        
    if plotL :
        ax = f1.add_subplot(nplots, 1, nplots)
        ax.plot(logL[indi::sampling], color = 'k')
        ax.set_ylabel('ln L')
    

    p.draw()
    return

"""
def compare_hist(vds1, vds2, keys = 'all', BI = 0.5):
    start_index = n.round(BI*len(vds1[0][vds1[0].keys()[0]]))
    if keys == 'all':
        keys = vds1[0].keys()

    for kk in keys:
        if kk in vd:
            f1 = p.figure()
            ax1 = f1.add_subplot(111)

            #
            ax1.hist(vd[kk], 50, histtype = 'step', color = 'b', lw = 2,
                     label = 'OldVersion', normed = True)
            
            if len(vds1) >= 1:
                for vdd in vds1:
                    ax1.hist(vdd[kk][start_index:], 50, histtype = 'step',
                             color = 'k', lw = 1, normed = True)
                ax1.patches[-1].set_label('NormalPrior')

            if len(vds2) >= 1:
                for vdd in vds2:
                    ax1.hist(vdd[kk][start_index:], 50, histtype = 'step',
                             color = 'r', lw = 1, normed = True)
                ax1.patches[-1].set_label('WidePrior')

            ax1.legend(loc = 0)
            ax1.set_xlabel(kk, fontsize = 16)

    return
"""
