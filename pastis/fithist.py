def fithist(y, par0, nbins = 50, disttype = 'Normal', plot = True, **kwargs):
    """
    Fit histogram to a probability distribution of type disttype

    Parameters
    ----------

    y: array
        Array of data to fit

    par0: iterable
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

    xlabel = kwargs.pop('xlabel', '')

    distdict = priors.distdict

    # Construct histogram
    ni, binedges = n.histogram(y, bins=nbins)
    x = binedges[:-1] + 0.5*n.diff(binedges)

    # Estimated number of events
    distdict[disttype][0].pdf(x, *par0) * n.diff(binedges)[0] * len(y)

    x = x[ni != 0]
    ni = ni[ni != 0]
    # Total number of observed events
    ntotal = n.sum(ni)
    
    # Compute error in bin numbers
    #eN = np.sqrt( N / (len(y) * (x[1] - x[0])) + N**2 / len(y) )

    def residuals(par, x, ni, bin_edges):

        # Estimated number of events
        ni0 = distdict[disttype][0].pdf(x, *par) * n.diff(bin_edges)[0] * ntotal

        return (ni - ni0) / n.sqrt(ni0)

    popt, epopt, infodict, mesg, wf = optimize.leastsq(residuals,
                                                       par0,
                                                       args=(x, ni, binedges),
                                                       full_output = True,
                                                       **kwargs)

    chi2 = n.sum(infodict['fvec']**2)
    ndof = len(ni) - 1
    print('Chi2: %.2f' % chi2)
    print('Number of degrees of freedom: %d' % ndof)
    print('p-value: %.2e' % stats.chi2.sf(chi2, ndof))

    if plot:
        import pylab as p
        f1 = p.figure()
        ax = f1.add_subplot(111)

        ni0 = distdict[disttype][0].pdf(x, *popt) * \
              n.diff(binedges)[0] * ntotal

        ax.hist(y, nbins, histtype='step', lw=2, color='k', normed=False)
        ax.errorbar(x, ni, n.sqrt(ni0), fmt='ko')
        
        xx = n.linspace(x.min(), x.max(), 1000)
        ax.plot(xx, distdict[disttype][0].pdf(xx, *popt) *
                n.diff(binedges)[0] * ntotal, lw=2, color='r')

        ax.set_xlabel(xlabel)

        p.draw()
                
    return popt, stats.chi2.sf(chi2, ndof)
