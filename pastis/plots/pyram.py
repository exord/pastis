import numpy as n
from scipy import stats
import pylab as p

from ..photometry import Allpbands


def pyramid(C, BI=0.0, Nlast=None, sample=1, label=None, filename=None):

    fig = p.figure(figsize = (10, 10))
    
    vd = C.get_value_dict()

    nBI = n.round(len(vd.values()[0])*BI, 0)
    if Nlast is not None:
        nf = Nlast
    else:
        nf = len(vd.values()[0])+1
    
    Nvar = len(vd.keys())

    fig.subplots_adjust(bottom= 0.1, top=0.95, hspace = 0.0, wspace = 0.0)

    axs = []
    for i, key1 in enumerate(n.sort(vd.keys())[1:]):
        for j, key2 in enumerate(n.sort(vd.keys())[:-1]):
            if j > i: continue
            ax = fig.add_subplot(Nvar-1, Nvar-1, i*(Nvar - 1) + j + 1)
            axs.append(ax)
            ax.plot(vd[key2][nBI: nf: sample], vd[key1][nBI: nf: sample], 'k,')

            # Set position of major ticks
            xl = ax.get_xlim(); dxl = xl[1]-xl[0]
            xts = n.arange(xl[0], xl[1], dxl/4.0)
            ax.xaxis.set_major_locator(p.FixedLocator(xts[1:]))

            yl = ax.get_ylim(); dyl = yl[1]-yl[0]
            yts = n.arange(yl[0], yl[1], dyl/4.0)
            ax.yaxis.set_major_locator(p.FixedLocator(yts[1:]))

            if i != len(vd.keys())-2:
                ax.xaxis.set_major_formatter(p.NullFormatter())
            else:
                ## Trick to get nice ticklabels
                s = '%e'%dxl
                a = s[s.find('e')+1:] # Power of range of axis in scientific notation
                if a[0] == '-':
                    ndecimal = int(a[1:]) + 1
                else:
                    ndecimal = 1

                ax.xaxis.set_major_formatter(p.FixedFormatter(n.round(xts[1:],
                                                                      ndecimal)
                                                              )
                                             )
                
                for xtl in ax.get_xticklabels():
                    xtl.set_rotation(90)
                    xtl.set_fontsize(11)

            # Modify key
            key2 = key2[key2.find('_'):]
            ax.set_xlabel(key2, fontsize = 10)

            if j != 0:
                ax.yaxis.set_major_formatter(p.NullFormatter())
            else:
                # Trick to get nice ticklabels
                s = '%e'%dyl
                a = s[s.find('e')+1:] # Power of range of axis in scientific notation
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
                    
                ax.set_ylabel(key1, fontsize = 10)

    labeltext = label+'\nNiteration = %d;\nBurnIn = %d%%'%(nf -1, BI*100)
    fig.text(0.6, 0.75, labeltext, fontsize = 16, ha = 'left', va = 'center')
    if filename == None : 
        p.show()
    else : 
        p.savefig(filename, bbox_inches = 'tight')
    
    return 


def pyramid_cont(C, bi=None, nlast=None, sample=1, filename=None,
                 **kwargs):
    """
    Plot the correlation diagram (a.k.a. the pyramid) for all parameters
    in a given chain.

    Parameters
    ----------
    C: dict, VDchain instance, or Chain instance
        The chain for which the pyramid plot will be done. The actual Chain
    instance or a dictionnary containing the values of the paremeters can
    be given.

    bi: float or None
        The burnin (in number of steps) to discard from the plot. This is used
        if the traces contain the burn-in period.

    nlast: int or None, optional
        Last element of the trace to be used for the plot. If not given, all
        the chain (potentially after the burn-in period) will be used.

    sample: int
        The thinning factor for the chains. Only one element over "sample" ones
        will be taken for the pyramid construction.


    Additional kwargs
    -----------------
    addhists, boolean
        if True, histograms for each parameter are added on the top of each
        column of the pyramid (an to the right of the bottom line).

    nbins, int
        number of bins in each direction used to compute confidence regions
        (default: 10)
    
    fontsize, int
        numbers font size

    labelfontsize, int
        labels font size

    objectlabel, boolean
        If false, the label of the object are deleted from the parameter label
        in the plot

    colorlist, tuple
        a tuple with the colors used for the contour plots. The color codes
        must be understood by the pylab module.

    confidence_levels, tuple
        a tuple with the confidence levels to plots
        (default: 0.6827, 0.9545, 0.9973).
        An arbitrary number of levels is supported.

    mapdict, dict
        an optional dictionary with the values to makr for each parameter.

    mapfmt, dict
        an optional dictionary with keywords parameters passed to axvline.
    
    label, string
        The label to append to the plot.

    overplot, boolean
        If true, the pyramid is done in the current figure.

    frequentist, boolean.
        If True contours are chosen according to the frequentist approach,
        as detailed in Sect. 9.7.6 of Frodesen.

    
    """
    
    ##
    # Get parameters from kwargs
    addhists = kwargs.pop('addhists', True)

    nbins = kwargs.pop('nbins', 10)
    fontsize = kwargs.pop('fontsize', 11)
    labelfontsize = kwargs.pop('labelfontsize', 10)    
    objectlabel = kwargs.pop('objectlabel', True)    
    colorlist = kwargs.pop('colorlist', ('0.85', '0.65', '0.45'))
    cmap = kwargs.pop('cmap', None)
    interpolation = kwargs.pop('interpolation', 'none')
    confidence_levels = kwargs.pop('confidence_levels',
                                   (0.6827, 0.9545, 0.9973))
    mapdict = kwargs.pop('mapdict', None)
    mapfmt = kwargs.pop('mapfmt', {})
    label = kwargs.pop('label', '')
    overplot = kwargs.pop('overplot', False)
    frequentist = kwargs.pop('frequentist', False)
    transpose = kwargs.pop('transpose', False)

    # Prepare map line fmt dictionary
    mapdefault = {'color': 'r',
                  'ls': ':',
                  'lw': 2,
                  'alpha': 1}
    mapdefault.update(mapfmt)
    mapfmt = mapdefault.copy()

    if overplot:
        fig = p.gcf()
    else:
        fig = p.figure(figsize = (8.5, 8.5))
    
    
    if isinstance(C, dict):
        vd = C.copy()
    else:
        vd = C.get_value_dict()

    # Check argument types
    assert (type(bi) is int) or (bi is None), "bi must be an integer or None"
    assert (type(nlast) is int) or (nlast is None), ("nlast must be an integer"
                                                      " or None")
    
    Nvar = len(vd.keys())

    fig.subplots_adjust(bottom= 0.15, left = 0.15, top=0.9, right=0.9,
            hspace = 0.0, wspace = 0.0)

    axs = []

    if transpose:
        keys1 = n.sort(vd.keys())[:-1]
        keys2 = n.sort(vd.keys())[1:]
    else:
        keys2 = n.sort(list(vd.keys()))[:-1]
        keys1 = n.sort(list(vd.keys()))[1:]
        
    for i, key1 in enumerate(keys1):
        for j, key2 in enumerate(keys2):
            if j > i: continue
            ax = fig.add_subplot(Nvar - 1, Nvar - 1, i*(Nvar - 1) + j + 1)
            axs.append(ax)

            #ax.plot(vd[key2][nBI: nf: sample], vd[key1][nBI: nf: sample], 'k,')
            H, xedges, yedges = n.histogram2d(vd[key2][bi: nlast: sample],
                                              vd[key1][bi: nlast: sample],
                                              bins = (nbins, nbins),
                                              normed = False)

            if frequentist:
                # Find Frequentist contour limits
                lims = find_frequentist_limits(Nvar, confidence_levels)
            else:
                # Find Bayesian contour limits
                lims = find_limits(H, confidence_levels)

            lims = lims[::-1]; lims.append(n.sum(H))

            ax.contourf(0.5*(xedges[:-1] + xedges[1:]),
                        0.5*(yedges[:-1] + yedges[1:]),
                        H.T, lims, cmap=None, colors=colorlist[::-1])

            """
            # Imshow version (WORKS!)
            ax.imshow(H.T, extent = (xedges[0], xedges[-1],
                                    yedges[0], yedges[-1]),
                        interpolation=interpolation, aspect='auto',
                        cmap=cmap, origin='lower')
            """
            
            """
        ax.contourf(xedges[1:], yedges[1:],
                        H.T, lims, colors = colorlist[::-1])
        """
            # Mark value
            if mapdict is not None and key2 in mapdict:
                ax.axvline(mapdict[key2], **mapfmt)
            
            # Set position of major ticks
            xl = ax.get_xlim(); dxl = xl[1]-xl[0]
            xts = n.arange(xl[0], xl[1], dxl/4.0)
            ax.xaxis.set_major_locator(p.FixedLocator(xts[1:]))

            yl = ax.get_ylim(); dyl = yl[1]-yl[0]
            yts = n.arange(yl[0], yl[1], dyl/4.0)
            ax.yaxis.set_major_locator(p.FixedLocator(yts[1:]))

            if i != len(vd.keys())-2:
                ax.xaxis.set_major_formatter(p.NullFormatter())
                ax.xaxis.set_ticks_position('none')
            else:
                ## Trick to get nice ticklabels
                s = '%e'%dxl
                a = s[s.find('e')+1:] # Power of range of axis in scientific notation
                if a[0] == '-':
                    ndecimal = int(a[1:]) + 1
                else:
                    ndecimal = 1

                ax.xaxis.set_major_formatter(p.FixedFormatter(n.round(xts[1:],
                                                                      ndecimal)
                                                              )
                                             )
                ax.xaxis.set_ticks_position('bottom')

                
                for xtl in ax.get_xticklabels():
                    xtl.set_rotation(90)
                    xtl.set_fontsize(fontsize)

                
                if objectlabel:
                    key2m = key2.replace('_', '\n')

                    key2m = key2m.replace('fitobsPlanet1\n','')
                    key2m = key2m.replace('CoRoT-W1','CoRoT')
                    key2m = key2m.replace('HARPS-G2','HARPS')
                    
                else:
                    objk = key2[:key2.find('_') - 1]
                    if not objk in Allpbands:
                        print(objk)
                        # Modify key
                        key2m = key2[key2.find('_')+1:]#.title()
                    else:
                        key2m = key2.replace('_', '\n')

                ax.set_xlabel(key2m, fontsize = labelfontsize)

            if j != 0:
                ax.yaxis.set_major_formatter(p.NullFormatter())
                ax.yaxis.set_ticks_position('none')
            else:
                # Trick to get nice ticklabels
                s = '%e'%dyl
                a = s[s.find('e')+1:] # Power of range of axis in scientific notation
                if a[0] == '-':
                    ndecimal = int(a[1:]) + 1
                else:
                    ndecimal = 1
                
                ax.yaxis.set_major_formatter(p.FixedFormatter(n.round(yts[1:], ndecimal)))
                ax.yaxis.set_ticks_position('left')
                for ytl in ax.get_yticklabels():
                    ytl.set_fontsize(fontsize)
                    
                if objectlabel:
                    key1m = key1.replace('_', '\n')

                    key1m = key1m.replace('fitobsPlanet1\n','')
                    key1m = key1m.replace('CoRoT-W1','CoRoT')
                    key1m = key1m.replace('HARPS-G2','HARPS')

                else:
                    if not key1[:key1.find('_') - 1] in Allpbands:
                        # Modify key
                        key1m = key1[key1.find('_')+1:]#.title()
                    else:
                        key1m = key1.replace('_', '\n')

                ax.set_ylabel(key1m, fontsize = labelfontsize)


            ### Add plot of histogram above top panels
            if addhists:
                if i == j:
                    axpos = ax.get_position().get_points()
                    width = axpos[1, 0] - axpos[0, 0]
                    height = axpos[1, 1] - axpos[0, 1]
                    axh = fig.add_axes((axpos[0, 0], axpos[1, 1], width,
                    height/2)
                    )

                    axh.hist(vd[key2][bi: nlast: sample], nbins,
                                     histtype = 'step', color = 'k')
                            
                    axh.set_xlim(xl)
                    axh.xaxis.set_ticks_position('bottom')
                    axh.yaxis.set_ticks_position('none')
                    axh.xaxis.set_major_formatter(p.NullFormatter())
                    axh.yaxis.set_major_formatter(p.NullFormatter())
                    axh.xaxis.set_major_locator(p.FixedLocator(xts[1:]))
    
                    if mapdict is not None and key2 in mapdict:
                        axh.axvline(mapdict[key2], **mapfmt)

                if i == len(vd.keys())-2: #Last row
                    axh2 = fig.add_axes((axpos[1, 0], axpos[0, 1], width/2,
                                         height))

                    axh2.hist(vd[key1][bi: nlast: sample], nbins, 
                              histtype='step', color='k', 
                              orientation='horizontal')
                    axh2.set_ylim(yl)
                    axh2.xaxis.set_ticks_position('none')
                    axh2.yaxis.set_ticks_position('left')
                    axh2.xaxis.set_major_formatter(p.NullFormatter())
                    axh2.yaxis.set_major_formatter(p.NullFormatter())
                    axh2.yaxis.set_major_locator(p.FixedLocator(yts[1:]))

                    if mapdict is not None and key1 in mapdict:
                        axh2.axhline(mapdict[key1], **mapfmt)

    #labeltext = label+'\nNiteration = %d;\nBurnIn = %d%%'%(nf -1, BI*100)
    labeltext = label
    fig.text(0.6, 0.75, labeltext, fontsize = 16, ha = 'left', va = 'center')
    if filename == None : 
        p.show()
    else : 
        p.savefig(filename, bbox_inches = 'tight')
    return fig

def find_limits(H, conflevels):
    """
    Find the limits for a pyramid plot using a Bayesian approach.
    """

    N = n.sum(H)

    # Sort number of elements per 2D bin.
    ind = n.argsort(H.flatten())[::-1]
    Hs  = n.sort(H.flatten())[::-1]

    bin1 = []; numpoints1 = 0
    bin2 = []; numpoints2 = 0
    bin3 = []; numpoints3 = 0

    ###
    ii = 0
    numpoints = 0
    bins = []
    for j, conflevel in enumerate(conflevels):
        bins.append([])
        while numpoints <= conflevel*N:
            # Add bins until number of points adds up to conflevel*N
            numpoints += Hs[ii]
            bins[j].append(ind[ii])
            ii += 1

    lims = []
    for bb in bins:
        lims.append(H.flatten()[bb].min())
        
    return lims
    
    
def find_frequentist_limits(k, conflevels):
    """
    Find the limits for a pyramid plot using a Frequentist approach.
    """
    #
    # Given k degrees of freedom, find levels at which likelihood
    # function must be cut
    lims = []
    for cl in conflevels:
        lims.append(0.5*stats.chi2.ppf(cl, k))
    return lims



