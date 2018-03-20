#!/usr/bin/env python

#Plot tools for pastis

## Setting matplotlib to use Tex
from matplotlib import rc
rc('text', usetex=True)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


import numpy as n
import pickle
import pylab as p
import colorsys
import itertools

from ..DataTools import readdata
from .. import photometry as phot
from ..models.PHOT import PASTIS_PHOT
from ..models.SED import PASTIS_SED, compute_global_spectrum
from ..models.RV import PASTIS_RV
from ..tools import rebin
## Load ObjectBuilder
from ..AstroClasses import *
from ..ObjectBuilder import ObjectBuilder
from ..MCMC.tools import state_constructor

aacolumnwidth = 255.76535 # column width of aa template in pts
aatextwidth = 523.5307 # text width of aa template in pts

inches_per_point = 1.0/72.27 # according to http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

golden_mean = (n.sqrt(5)-1.0)/2.0 # To obtain aesthetic figures

"""
params = {'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'lines.markersize': 5,
          'interactive': True
          }
p.rcParams.update(params)
"""

#p.rcdefaults()

def plot_dict(dictionary, **kwargs):
    """
    plot a dictionary (.pastis) created with the widget: data, model 
    Example: plot_dict('/data/PASTIS/configfiles/test.pastis')

    Parameters
    ----------
    dictionary: str
       A full path and filename of the solution file to plot. This file is
       can be produced with MCMC.analysis.make_solution_file and contains a
       pickled list with all configuration dictionaries.


    Additional keyword parameters
    -----------------------------

    - Parameters controlling the size of figures and axes

    axsize (l, b w, h): list
        A 4-element list used to control the size of the main axes in the
        figures.

    ressize (l, b, w, h'): list
        Similar for the residual lower panel.

    figsize (w, h): list
        Size of figures in inches.

    ticksize: int
        Size of ticklabels.
        
    nplots: int
        Number of plots per column in A&A template (default = 1).
        If 0, use figsize argument to set size of figures.

    nplotst: int
        Number of plots per pagewidth in A&A template (default = 0).
        If 0, nplots argument to set size of figures.


    - Parameters for controlling the data that is plotted

    istransit: bool
        Indicating whether the phase should be considered from transit time.
        If false, use periastron passage. Default: true

    usehours: bool
        Express time in hours from central time or periastron passage in
        phase-folded plot.

    timeoffset: float
        This number is substracted from all times in plots.

    lcbin: bool
        Control if the phase-folded binned LC data are plotted.

    lcbinnewfig: bool
        Control if the phase-folded binned LC data is plotted in the same
        axes as the unbinned data.

    lcbinsize: float
        Bin size in orbital phase.

    npoints: int
        Number of points of model curves.

    padfrac: float
        Fraction of the entire time span of data to add at each side of the
        time plots.

    mergeLC: bool
        Control if all light curves from same filter are plotted in the same
        axes.

    plotsecondary: bool
        If false, primary transit is at centre of phase array. Otherwise,
        array goes from 0 to 1, leaving phase 0.5 at the centre.

    correct_contam: bool
        Control if contamination is corrected in LC data and models.

    unbinned: bool
        Control if unbinned data is plotted

    rangephLC: list
        2-element list with the phase range of the light curve phase plot
 
    mergeRV: bool
        Control if all radial velocity curves are plotted in the same axes.

    mRVrandomcolors: bool
        random colors in the merged RV plot
    
    mRVrandommarkers: bool
        random markers in the merged RV plot

    histresRV: int
        number of bins of the histogram of the residuals of the radial velocity data. 
        Default is None and nothing is ploted
        
    rangetRV: list
        2-element list with the time range of the radial velocity plot in time 

    considerTriple: bool
        plot Triple orbit

    - Parameters controlling the format of plots

    rvfmt, lcfmt, sedfmt: dict
        Dictionaries containing format codes for RV, LC, and SED plots.
    Possible keys are:

    fmt: str
        Format data points, including bisector, fwhm, etc.
    lfmt: str
            Format for RV model curves
    mec: str
            Color of marker edges in RV plots.
    capsize: int
        Size of errorcaps in RV plots.
    lw: int
            Line weight of model curves.
    binfmt: str
            Format for binned points.
    binmarkersize: str
            Format for binned points.
        binmfc: str
            Marker face color for binned points.
    resmajor: float
        The major tick step size in residual plots
    resminor: float
        The minor tick step size in residual plots
    modelontop: bool
            Control if the model is above or below the data points.
        labelsize: int
            Size of x and y axis labels in points
        reshlcol: str
            color of the horizontal line in the residuals plot
        
    """
    rvdefault = {'fmt': 'or',
         'lfmt': 'k',
         'mec': 'r',
         'capsize': 0,
         'lw' : 2,
         'modelontop' : False,
                 'labelsize' : 16
         }

    lcdefault = {'fmt': ',k',
         'lfmt': 'r',
         'mec': 'r',
         'capsize': 0,
         'lw' : 2,
         'binfmt' : 'go',
         'binmarkersize' : 8,
                 'binmfc' : 'g',
         'resmajor' : 2e3,
         'resminor' : 1e3,
         'modelontop' : False,
                 'labelsize' : 16,
                 'reshlcol' : 'r'
         }

    seddefault = {'fmt': 'or',
          'lfmt': 'k',
          'mec': 'r',
          'capsize': 0,
          'lw' : 2,
          'binfmt' : 'g.',
          'binmarkersize' : 8,
                  'binmfc' : 'g',
          'resmajor' : 5e-2,
          'resminor' : 2.5e-2,
                  'modelontop' : False,
                  'labelsize' : 16
          }

    rvfmt = kwargs.pop('rvfmt', None)
    lcfmt = kwargs.pop('lcfmt', None)
    sedfmt = kwargs.pop('sedfmt', None)

    if rvfmt is not None:
        rvdefault.update(rvfmt)
    if lcfmt is not None:
        lcdefault.update(lcfmt)
    if sedfmt is not None:
        seddefault.update(sedfmt)

    rvfmt = rvdefault.copy()
    lcfmt = lcdefault.copy()
    sedfmt = seddefault.copy()
    

    npoints = kwargs.pop('npoints', 400)
    timeoffset = kwargs.pop('timeoffset', 0.0)
    padfrac = kwargs.pop('padfrac', 0.1)
    mergeLC = kwargs.pop('mergeLC', False)
    plotsecondary = kwargs.pop('plotsecondary', False)
    correct_contam = kwargs.pop('correct_contam', True)
    unbinned = kwargs.pop('unbinned', True)
    rangephLC = kwargs.pop('rangephLC', None)
    mergeRV = kwargs.pop('mergeRV', False)
    mRVrandomcolors = kwargs.pop('mRVrandomcolors', False)
    mRVrandommarkers = kwargs.pop('mRVrandommarkers', False)
    histresRV = kwargs.pop('histresRV', None)
    rangetRV = kwargs.pop('rangetRV', None)
    lcbin = kwargs.pop('lcbin', True)
    lcbinnewfig = kwargs.pop('lcbinnewfig', False)
    lcbinsize = kwargs.pop('lcbinsize', 0.01)
    considerTriple = kwargs.pop('considerTriple', False)
    
    ## SET SIZE OF FIGURES AND AXES
    axsize = kwargs.pop('axsize', [0.125, 0.18+0.125, 0.95-0.125, 0.95-0.305])
    #axsize = kwargs.pop('axsize', [0.0, 0.18, 1.0, 0.82])
    
    ressize = kwargs.pop('ressize', [0.125, 0.125, 0.95-0.125, 0.18])
    #ressize = kwargs.pop('ressize', [0.0, 0.0, 1.0, 0.18])

    fsize = kwargs.pop('figsize', (7.5, 7.5*golden_mean))
    ticksize = kwargs.pop('ticksize', 12)
    nplots = kwargs.pop('nplots', 0)
    nplotst = kwargs.pop('nplotst', 0)

    ## Prepare size of figures
    if nplotst != 0:
        figwidth = aatextwidth/nplotst*inches_per_point # In inches
        figheight = figwidth*golden_mean # In inches
        fsize = ([figwidth, figheight]) #inch
        print(fsize)
    elif nplots != 0:
        figwidth = aacolumnwidth/nplots*inches_per_point # In inches
        figheight = figwidth*golden_mean # In inches
        fsize = ([figwidth, figheight]) #inch
    else:
        pass

    ## ???
    istransit = kwargs.pop('istransit', True)
    usehours = kwargs.pop('usehours', False)

    # Define dictionary of axes for PHOT plots
    axdict = {}

    # Figure list
    figlist = []
    
    f = open(dictionary, 'r')
    dd = pickle.load(f)
    f.close()

    datadict, lc = readdata(dd[2])
    input_dict = dd[1].copy()

    # Construct initial state from info in input_dict
    X, labeldict = state_constructor(input_dict)

    # Build objects
    objects = ObjectBuilder(input_dict)

    # Get ephemeris of all objects
    periods = []
    tps = []

    for obj in objects:
        if isinstance(obj, PlanSys):
            for pl in obj.planets:
                periods.append(pl.orbital_parameters.P)
                if istransit:
                    tps.append(pl.orbital_parameters.T0)
                else:
                    tps.append(pl.orbital_parameters.Tp)

        if isinstance(obj, FitBinary):
            print(obj, istransit, obj.orbital_parameters.Tp, obj.orbital_parameters.T0, obj.orbital_parameters.P, obj.orbital_parameters.ecc, obj.orbital_parameters.omega*180/pi)
            
            periods.append(obj.orbital_parameters.P)
            if istransit:
                tps.append(obj.orbital_parameters.T0)
            else:
                tps.append(obj.orbital_parameters.Tp)

        if isinstance(obj, Triple):
            if considerTriple:
                periods.append(obj.orbital_parameters.P)
                if istransit:
                    tps.append(obj.orbital_parameters.T0)
                else:
                    tps.append(obj.orbital_parameters.Tp)

            for component in (obj.object1, obj.object2):
                if isinstance(component, FitBinary):
                    periods.append(component.orbital_parameters.P)
                    if istransit:
                        tps.append(component.orbital_parameters.T0)
                    else:
                        tps.append(component.orbital_parameters.Tp)
                        

                elif isinstance(component, PlanSys):
                    for pl in component.planets:
                        periods.append(pl.orbital_parameters.P)
                        if istransit:
                            tps.append(pl.orbital_parameters.T0)
                        else:
                            tps.append(pl.orbital_parameters.Tp)



    ########
    ## PLOTS
    ########

    #plot phase merged RV 
    if mergeRV:
        
        keysRV = []
        for key in datadict.keys():
            if datadict[key]['type'] == 'RV':
                keysRV.append(key)
                            
        keysRV.sort()
        for i in range(len(periods)):
 
            if rvfmt['modelontop']:
                zmod = 2
                zdat = 1
            else:
                zmod = 1
                zdat = 2

            # Plot phase
            f2 = p.figure(figsize = fsize)
            figlist.append(f2)
            ax3 = f2.add_axes(axsize)
            ax3.xaxis.set_major_formatter(p.NullFormatter())
                        
            ax3.set_xlim(0.0,1.0)
            
            ax3.set_ylabel('Radial Velocity [km s$^{-1}$]',
                           size = rvfmt['labelsize']
                           )
            ax3.set_title('Merged')
            
            #Residuals
            ax4 = f2.add_axes(ressize)
            
            ax4.axhline(0.0, ls=':', color='0.55', zorder =0)
            ax4.set_xlim(0.0,1.0)
            ax4.set_xlabel('Orbital Phase', size = rvfmt['labelsize'])
            ax4.set_ylabel('O-C [m s$^{-1}$]',
                           size = rvfmt['labelsize']
                           )

            # Plot tiempo
            f1 = p.figure(figsize = fsize)
            figlist.append(f1)
            ax = f1.add_axes(axsize)
            ax.xaxis.set_major_formatter(p.NullFormatter())
            ax2 = f1.add_axes(ressize)

            # Boucle for datasets to get all dates
            rvtime = n.array([])
            for key in keysRV:
                print(key)
            
                # Get information of instrument
                spectro = datadict[key]['instrument']
                mask = datadict[key]['mask']
        
                # Get value of instrument offset
                offset = dd[1][key]['offset'][0]

                # Get time array
                rvtime = n.concatenate((rvtime,datadict[key]['data']['time']))

            
            print('Total number of RV measurements: %d'%len(rvtime))
            # Construct model time array
            dt = rvtime.max() - rvtime.min()
            xx_time_rv = n.arange(rvtime.min() - padfrac*dt,
                                      rvtime.max() + padfrac*dt,
                                      dt/npoints)

            vrad_all, bis_all, fwhm_all, cont_all = PASTIS_RV(xx_time_rv,
                                                              [spectro, mask,
                                                               0],
                                                              *objects)
            xx_time_rv_all = xx_time_rv

            # plot model phase
            phase_all = ((xx_time_rv - tps[i])/periods[i])%1.0
            ax3.plot(n.sort(phase_all), vrad_all[n.argsort(phase_all)],
                     rvfmt['lfmt'], lw =rvfmt['lw'], zorder = zmod)

            # plot model time
            ax.plot(xx_time_rv - timeoffset, vrad_all, rvfmt['lfmt'],
                    lw = rvfmt['lw'], zorder = zmod
                    )
                #ax.set_xlabel('BJD - 2,450,000')
            ax.set_ylabel('Radial Velocity [km s$^{-1}$]',
                          size = rvfmt['labelsize'])
            ax.set_title('Merged')
            ax2.axhline(0.0, ls=':', color='0.55', zorder = 0)

            #Labels
            if timeoffset == 0:
                ax2.set_xlabel('BJD', size = rvfmt['labelsize'])

            else:
                ax2.set_xlabel('BJD - %d'%timeoffset,
                               size = rvfmt['labelsize']
                               )

            ax2.set_ylabel('O-C [m s$^{-1}$]', size = rvfmt['labelsize'])
            
            # Boucle for datasets
            residuals = n.array([])
            if mRVrandommarkers: 
                markers = itertools.cycle(('s', 'o','^','>','v','<','D')) 
            else: markers = itertools.cycle(('o'))
            
            chi2=0.0
            #markers = itertools.cycle(('s','o','D'))
            for key in keysRV:
                #print key
            
                # Get information of instrument
                spectro = datadict[key]['instrument']
                mask = datadict[key]['mask']
        
                # Get value of instrument offset
                offset = dd[1][key]['offset'][0]

                # Get time array
                rvtime = datadict[key]['data']['time']
        
                # Construct model time array
                dt = rvtime.max() - rvtime.min()
                xx_time_rv = n.arange(rvtime[0] - padfrac*dt,
                                      rvtime[-1] + padfrac*dt,
                                      dt/npoints)

                # Compute theoretical RV curves
                vrad, bis, fwhm, cont = PASTIS_RV(rvtime, [spectro, mask, offset],
                                                  *objects)
                                
                phase = ((rvtime - tps[i])/periods[i])%1.0
                    
                
                ### Plot phase
                """
                r=n.random.rand()
                g=n.random.rand()
                b=n.random.rand()
                color=(r,g,b)
                """
                if mRVrandomcolors: 
                    color=colorsys.hsv_to_rgb(n.random.rand(),1.0,0.75+0.25*n.random.rand())
                else: color='r'

                marker=markers.next()

                ax3.errorbar(phase, datadict[key]['data']['vrad']+offset,
                             datadict[key]['data']['svrad'],
                             fmt= rvfmt['fmt'],
                             capsize = rvfmt['capsize'], zorder = zdat,
                             color=color,marker=marker)
                                
                # Residuals
                ax4.errorbar(phase,
                             (datadict[key]['data']['vrad'] - vrad)*1e3,
                             datadict[key]['data']['svrad']*1e3,
                             fmt = rvfmt['fmt'], zorder = zdat,
                             capsize = 0,color=color,marker=marker
                             )
                residuals = n.concatenate((residuals,datadict[key]['data']['vrad'] - vrad))

                chi2=chi2+n.sum(((datadict[key]['data']['vrad'] - vrad)/datadict[key]['data']['svrad'])**2)

                ### Plot tiempo

        # Data
                ax.errorbar(rvtime - timeoffset, datadict[key]['data']['vrad']+offset,
                datadict[key]['data']['svrad'],  fmt = rvfmt['fmt'],
                capsize = rvfmt['capsize'], zorder = zdat,
                            color=color, marker=marker
                )

                

                #Residuals
                ax2.errorbar(rvtime - timeoffset,
                             (datadict[key]['data']['vrad']-vrad)*1e3,
                             datadict[key]['data']['svrad']*1e3,
                 fmt = rvfmt['fmt'], zorder = zdat,
                             capsize = rvfmt['capsize'], color=color,
                             marker=marker)

            # plot limits
            if rangetRV == None:
                ax.set_xlim(xx_time_rv_all.min() - timeoffset,
                            xx_time_rv_all.max() - timeoffset
                            )
                ax2.set_xlim(xx_time_rv_all.min() - timeoffset,
                             xx_time_rv_all.max() - timeoffset
                             )
            else:
                ax.set_xlim(rangetRV[0],rangetRV[1])
                ax2.set_xlim(rangetRV[0],rangetRV[1])
                
            
            # histogram of the residuals
            if histresRV !=None:
                f33 = p.figure()
                figlist.append(f33)
                ax33 = f33.add_axes(axsize)
                ax33.hist(residuals,normed='True',bins=histresRV)
    
                
            print('Chi2: %f'%chi2)
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # If merge data create new dict to hold merged data
    datadictmergeLC = {}

    ## Iterate for all datasets
    for key in datadict.keys():
             

        ######
        # RV #
        ######
        if datadict[key]['type'] == 'RV' and not(mergeRV):
            print(key)
            """try:
                spectro = datadict[key]['instrument']
            except KeyError:
                spectro = key
                
            """
            # Get information of instrument
            spectro = datadict[key]['instrument']
            mask = datadict[key]['mask']
        
            # Get value of instrument offset
            offset = dd[1][key]['RVoffset'][0]
            
            if 'CCF_width' in datadict[key]: 
                CCF_width = datadict[key]['CCF_width']
            else : CCF_width = 2.

            # Get time array
            rvtime = datadict[key]['data']['time']
        
            # Construct model time array
            dt = rvtime.max() - rvtime.min()
            
            if len(periods) == 1:
                ti = rvtime.min()- padfrac*periods[0]
                xx_time0 = n.linspace(ti, ti + periods[0] - periods[0]/npoints, npoints) 
                xx_time_rv = xx_time0

            else:
                xx_time_rv = n.arange(rvtime.min() - padfrac*dt,
                                      rvtime.max() + padfrac*dt,
                                      dt/npoints)
            RVdatadict = {'spectro': spectro, 'mask':mask, 'RV':{'offset':offset}}
            # Compute theoretical RV curves
            RVoutdict = PASTIS_RV(rvtime, RVdatadict, *objects)
            
            RVoutdict0 = PASTIS_RV(xx_time_rv,RVdatadict, *objects)
            
            vrad_all = RVoutdict0['RV']
            
            if len(periods) == 1:
                ## Duplicate array for plot
                jj = 1
                while n.max(xx_time_rv) < (rvtime.max() + padfrac*periods[0]):
                    xx_time_rv = n.concatenate((xx_time_rv, xx_time0 + jj*periods[0]))
                    vrad_all = n.concatenate((vrad_all, RVoutdict0['RV']))
                    #bis_all = n.concatenate((bis_all, bis_all0))
                    #fwhm_all = n.concatenate((fwhm_all, fwhm_all0))
                    #cont_all = n.concatenate((cont_all, cont_all0))
                    jj = jj + 1

            ################
            # Plots in time
            ################
        if 'vrad' in datadict[key]['data']:
            ####### RV ########
            f1 = p.figure(figsize = fsize)
            figlist.append(f1)
            ax = f1.add_axes(axsize)
            ax.xaxis.set_major_formatter(p.NullFormatter())

            if rvfmt['modelontop']:
                zmod = 2
                zdat = 1
            else:
                zmod = 1
                zdat = 2
            
            # Data
            ax.errorbar(rvtime - timeoffset, datadict[key]['data']['vrad'],
            datadict[key]['data']['svrad'],  fmt = rvfmt['fmt'],
            mec = rvfmt['mec'], capsize = rvfmt['capsize'],
            zorder = zdat
            )

            # Model
            ax.plot(xx_time_rv - timeoffset, vrad_all, rvfmt['lfmt'], 
                    lw = rvfmt['lw'], zorder = zmod)
                
            #ax.set_xlabel('BJD - 2,450,000')
            ax.set_ylabel('Radial Velocity [km s$^{-1}$]',
                          size = rvfmt['labelsize'])
            ax.set_title(key)

            #Residuals
            ax2 = f1.add_axes(ressize)
            ax2.errorbar(datadict[key]['data']['time'] - timeoffset,
                         (datadict[key]['data']['vrad']-RVoutdict['RV'])*1e3,
                         datadict[key]['data']['svrad']*1e3,
             fmt = rvfmt['fmt'], mec = rvfmt['mec'],
             capsize = rvfmt['capsize'])

            ax2.axhline(0.0, ls=':', color='0.55', zorder = 0)

            #Labels
            if timeoffset == 0:
                ax2.set_xlabel('BJD', size = rvfmt['labelsize'])
            else:
                ax2.set_xlabel('BJD - %d'%timeoffset,
                               size = rvfmt['labelsize']
                               )
                ax2.set_ylabel('O-C [m s$^{-1}$]', size = rvfmt['labelsize'])


            # plot limits
            if rangetRV == None:
                ax.set_xlim(xx_time_rv.min() - timeoffset,
                            xx_time_rv.max() - timeoffset
                            )
                ax2.set_xlim(xx_time_rv.min() - timeoffset,
                             xx_time_rv.max() - timeoffset
                             )
            else:
                ax.set_xlim(rangetRV[0],rangetRV[1])
                ax2.set_xlim(rangetRV[0],rangetRV[1])
            
            """
            if 'bis' in datadict[key]['data']:

                ####### BIS ########
                f1bis = p.figure(figsize = fsize)
                figlist.append(f1bis)
        ax = f1bis.add_axes(axsize)
                ax.xaxis.set_major_formatter(p.NullFormatter())

        if rvfmt['modelontop']:
            zmod = 2
            zdat = 1
        else:
            zmod = 1
            zdat = 2
            
        # Data
                ax.errorbar(rvtime - timeoffset, datadict[key]['data']['bis'],
                datadict[key]['data']['sbis'],  fmt = rvfmt['fmt'],
                mec = rvfmt['mec'], capsize = rvfmt['capsize'],
                zorder = zdat
                )
                
        # Model
                ax.plot(xx_time_rv - timeoffset, bis_all, rvfmt['lfmt'],
            lw = rvfmt['lw'], zorder = zmod
            )
                
                #ax.set_xlabel('BJD - 2,450,000')
                ax.set_ylabel('BIS [km s$^{-1}$]',
                              size = rvfmt['labelsize'])
                ax.set_title(key)

                #Residuals
                ax2 = f1bis.add_axes(ressize)
                ax2.errorbar(datadict[key]['data']['time'] - timeoffset,
                             (datadict[key]['data']['bis'] - bis)*1e3,
                             datadict[key]['data']['svrad']*1e3,
                 fmt = rvfmt['fmt'], mec = rvfmt['mec'],
                 capsize = rvfmt['capsize'])

        ax2.axhline(0.0, ls=':', color='0.55', zorder = 0)

                #Labels
        if timeoffset == 0:
            ax2.set_xlabel('BJD', size = rvfmt['labelsize'])
        else:
            ax2.set_xlabel('BJD - %d'%timeoffset,
                                   size = rvfmt['labelsize']
                                   )
                ax2.set_ylabel('O-C [m s$^{-1}$]', size = rvfmt['labelsize'])


                # plot limits
                if rangetRV == None:
                    ax.set_xlim(xx_time_rv.min() - timeoffset,
                                xx_time_rv.max() - timeoffset
                                )
                    ax2.set_xlim(xx_time_rv.min() - timeoffset,
                                 xx_time_rv.max() - timeoffset
                                 )
                else:
                    ax.set_xlim(rangetRV[0],rangetRV[1])
                    ax2.set_xlim(rangetRV[0],rangetRV[1])
            """
            """
            if 'fwhm' in datadict[key]['data']:

                ####### FWHM ########
                figure(np)
                clf()
                subplot(211)
                errorbar(datadict[key]['data']['time'], datadict[key]['data']['fwhm'], datadict[key]['data']['sfwhm'], fmt='o')
                plot(datadict[key]['data']['time'], fwhm, 'o')
                xlabel('time')
                ylabel('FWHM [km/s]')
                title(key)
                #Residuals
                subplot(212)
                errorbar(datadict[key]['data']['time'], datadict[key]['data']['fwhm']-fwhm, datadict[key]['data']['sfwhm'], fmt='o')
                axhline(linewidth=1, color='r')
                xlabel('time')
                ylabel('Residuals (data-model)')            
                draw()
                np=np+1
            """

            ################
            # Plots in phase
            ################
            for i in range(len(periods)):
                phase = ((rvtime - tps[i])/periods[i])%1.0
                phase_all = ((xx_time_rv - tps[i])/periods[i])%1.0
                if 'vrad' in datadict[key]['data']:

                    ####### RV ########
                    f2 = p.figure(figsize = fsize)
                    figlist.append(f2)
                    ax3 = f2.add_axes(axsize)
                    ax3.xaxis.set_major_formatter(p.NullFormatter())
                    ax3.errorbar(phase, datadict[key]['data']['vrad'],
                 datadict[key]['data']['svrad'],
                                 fmt= rvfmt['fmt'], mec = rvfmt['mec'],
                 capsize = rvfmt['capsize'], zorder = zdat
                 )

                    phase_all = phase_all[:npoints]
                    vrad_all = vrad_all[:npoints]

                    ax3.plot(n.sort(phase_all), vrad_all[n.argsort(phase_all)],
                 rvfmt['lfmt'], lw =rvfmt['lw'], zorder = zmod)
                    
                    ax3.set_xlim(0.0,1.0)
                    
                    ax3.set_ylabel('Radial Velocity [km s$^{-1}$]',
                                   size = rvfmt['labelsize']
                                   )
                    ax3.set_title(key)

                    #Residuals
                    ax4 = f2.add_axes(ressize)
                    ax4.errorbar(phase,
                 (datadict[key]['data']['vrad'] - RVoutdict['RV'])*1e3,
                 datadict[key]['data']['svrad']*1e3,
                 fmt = rvfmt['fmt'], mec = rvfmt['mec'],
                 capsize = 0
                 )
            
                    ax4.axhline(0.0, ls=':', color='0.55', zorder =0)
                    ax4.set_xlim(0.0,1.0)
                    ax4.set_xlabel('Orbital Phase', size = rvfmt['labelsize'])
                    ax4.set_ylabel('O-C [m s$^{-1}$]',
                                   size = rvfmt['labelsize']
                                   )

                """
                if 'bis' in datadict[key]['data']:

                    ####### BIS ########
                    figure(np)
                    clf()
                    subplot(211)
                    errorbar(phase, datadict[key]['data']['bis'], datadict[key]['data']['sbis'], fmt='o')
                    plot(phase, bis, 'o')
                    xlim(0.0,1.0)
                    xlabel('phase (Tp)')
                    ylabel('BIS [km/s]')
                    title(key+', P = '+str(periods[phases])+' days')
                    #Residuals
                    subplot(212)
                    errorbar(phase, datadict[key]['data']['bis']-bis, datadict[key]['data']['sbis'], fmt='o')
                    axhline(linewidth=1, color='r')
                    xlim(0.0,1.0)
                    xlabel('phase (Tp)')
                    ylabel('Residuals (data-model)')            
                    draw()
                    np=np+1
                    
                if 'fwhm' in datadict[key]['data']:
                    
                    ####### FWHM ########
                    figure(np)
                    clf()
                    subplot(211)
                    errorbar(phase, datadict[key]['data']['fwhm'], datadict[key]['data']['sfwhm'], fmt='o')
                    plot(phase, fwhm, 'o')
                    xlim(0.0,1.0)
                    xlabel('phase (Tp)')
                    ylabel('FWHM [km/s]')
                    title(key+', P = '+str(periods[phases])+' days')
                    #Residuals
                    subplot(212)
                    errorbar(phase, datadict[key]['data']['fwhm']-fwhm, datadict[key]['data']['sfwhm'], fmt='o')
                    axhline(linewidth=1, color='r')
                    xlim(0.0,1.0)
                    xlabel('phase (Tp)')
                    ylabel('Residuals (data-model)')            
                    draw()
                    np=np+1
                """

        ########
        # PHOT #
        ########
        if datadict[key]['type'] == 'PHOT':
            print(key)

        if mergeLC:
            ## Change key to name of filter
            dkey = key[:-1]
               
            if dkey in datadictmergeLC:
                continue
            else:
                # Merge data!
                datadictmergeLC[dkey] = datadict[key].copy()
                ddm = datadictmergeLC[dkey]
                for photband in datadict.keys():
                    if photband[:-1] != dkey: 
                        continue
                        
                    # Get value of contamination and Foot
                    cont = labeldict[photband+'_contamination'].get_value()
                    foot = labeldict[photband+'_foot'].get_value()
                    print(photband,cont,foot)

                    for jj in datadict[photband]['data'].keys():
                        if 'flux' in jj:
                            yi = datadict[photband]['data'][jj]/foot

                            if jj == 'sflux':
                                ## Error decontamination
                                ## Neglects error in cont and foot
                                yic = yi/(1.0 - cont)
                            else:
                                yic = yi/(1.0 - cont) - cont/(1.0 - cont)
                        else:
                            yic = datadict[photband]['data'][jj]
                                
                            ddm['data'][jj] = n.concatenate((ddm['data'][jj],
                                                             yic))
                        if 'overtime' in ddm:
                            ddm['overtime'] = n.concatenate((ddm['overtime'], datadict[photband]['overtime']))

                datadictLC = datadictmergeLC

        else:
        dkey = key
                datadictLC = datadict
         
            try:
                filtre = datadict[dkey]['filter']
            except KeyError:
                filtre = dkey
 
            if not mergeLC:
                # Get value of contamination and Foot
                cont = labeldict[key+'_contamination'].get_value()
                foot = labeldict[key+'_foot'].get_value()

            # If filter is one of CoRoT colors
            if filtre == 'CoRoT-W':  
                lcfmt['binfmt']='ko'
                lcfmt['binmfc']='k'
                lcfmt['lfmt']='r'
                lcfmt['reshlcol']='r'
            if filtre == 'CoRoT-R':  
                lcfmt['binfmt']='ro'
                lcfmt['binmfc']='r'
                lcfmt['lfmt']='k'
                lcfmt['reshlcol']='k'
            if filtre == 'CoRoT-G':  
                lcfmt['binfmt']='go'
                lcfmt['binmfc']='g'
                lcfmt['lfmt']='k'
                lcfmt['reshlcol']='k'
            if filtre == 'CoRoT-B':  
                lcfmt['binfmt']='bo'
                lcfmt['binmfc']='b'
                lcfmt['lfmt']='k'
                lcfmt['reshlcol']='k'

            if (filtre == 'CoRoT-R' or filtre == 'CoRoT-G' or filtre == 'CoRoT-B'):
                compute_global_spectrum(*objects)
                from ..models.SED import global_spectrum

                for key2 in datadictLC.keys():
                    if 'filter' in datadict[key2]:
                        if datadictLC[key2]['filter'] == 'CoRoT-R':
                            meanR = datadict[key2]['MeanFlux']
                            contR = labeldict[key2+'_contamination'].get_value()

                        if datadictLC[key2]['filter'] == 'CoRoT-G':
                            meanG = datadict[key2]['MeanFlux']
                            contG = labeldict[key2+'_contamination'].get_value()

                        if datadictLC[key2]['filter'] == 'CoRoT-B':
                            meanB = datadictLC[key2]['MeanFlux']
                            contB = labeldict[key2+'_contamination'].get_value()
                    
                # Add limbdarkening weights and filters for CoRoT colors
                if 'global_spectrum' not in locals():
                    compute_global_spectrum(*objects)
                phot.corot_colors(meanR, meanG, meanB, contR, contG, contB)
                #has_computed_colors = True

            # Check if need to compute oversampled lightcurve
        lctime = datadictLC[dkey]['data']['time']
            if datadictLC[dkey]['sampling'] > 1.0:

                olctime = datadictLC[dkey]['overtime']

                ksampl = datadictLC[dkey]['sampling']
                
                if False: #len(olctime) < npoints:
                    dt = olctime.max() - olctime.min()
                    olctimet = n.arange(n.min(olctime), n.max(olctime),
                                        dt/(npoints*ksampl)
                                        )
                    kfactor = len(olctimet)/len(olctime)
                else:
                    olctimet = olctime
                    
        # Compute theoretical OVERSAMPLED lightcurve
                osft = PASTIS_PHOT(olctimet, filtre,
                                   datadictLC[dkey]['is_phase'],
                   0.0, 1.0, *objects
                   )
                
                # Rebin model flux to original sampling rate
                ft = rebin(osft, ksampl)
                lctimeb =  rebin(olctimet, ksampl)

                # complete model lightcurve
                dt = n.max(lctime)-n.min(lctime)
                oslctimetAll = n.arange(n.min(olctime), n.max(olctime),olctime[1]-olctime[0])
                osftAll = PASTIS_PHOT(oslctimetAll, filtre, datadictLC[dkey]['is_phase'],
                                    0.0, 1.0, *objects)
                # Rebin model flux to original sampling rate
                ftAll = rebin(osftAll, ksampl)
                lctimetAll =  rebin(oslctimetAll, ksampl)

            else:
                if False: #len(lctime) < npoints:
                    dt = lctime.max() - lctime.min()
                    lctimet = n.arange(n.min(lctime), n.max(lctime),
                                       dt/npoints)

                else:
                    lctimet = lctime

                # Compute theoretical light curves and likelihood
                ft = PASTIS_PHOT(lctimet, filtre, datadictLC[dkey]['is_phase'],
                 0.0, 1.0, *objects
                 )
        lctimeb = lctimet

                # complete model lightcurve
                dt = n.max(lctime)-n.min(lctime)
                lctimetAll = n.arange(n.min(lctime)- padfrac*dt, 
                                        n.max(lctime)+ padfrac*dt, lctime[1]-lctime[0])
                ftAll = PASTIS_PHOT(lctimetAll, filtre, datadictLC[dkey]['is_phase'],
                                    0.0, 1.0, *objects)

            '''
        # Normalize theoretical light curve
        ft = ft/foot
            ftAll = ftAll/foot
        
        if correct_contam:
        ft = ft/(1.0 - cont) - cont/(1.0 - cont)
        ftAll = ftAll/(1.0 - cont) - cont/(1.0 - cont)
            '''

        #res = datadict[dkey]['data']['flux'] - ft

            # Compute light theoretical lightcurve with a
            """
            phtt = n.arange(-0.05, 0.05, 0.001)
            ftt = PASTIS_PHOT(phtt, filtre, True, cont, foot, *objects)
            """

            if not datadictLC[dkey]['is_phase']:
                """
                #LC in time
                f3 = p.figure()
                ax5 = f3.add_axes([0.1, 0.3, 0.85, 0.65])
                ax5.errorbar(datadict[dkey]['data']['time'], datadict[dkey]['data']['flux'], datadict[dkey]['data']['sflux'], fmt=',',color='k')
                ax5.plot(datadict[dkey]['data']['time'],ft,linewidth=2,color='r')

                ax5.xaxis.set_major_formatter(NullFormatter())
                ax5.set_ylabel('Normalized flux')
                ax5.set_title(key)
                
                ax6 = f3.add_axes([0.1, 0.12, 0.85, 0.18])
                ax6.errorbar(datadict[dkey]['data']['time'], datadict[dkey]['data']['flux']-ft, datadict[dkey]['data']['sflux'], fmt=',',color='k')
                ax6.axhline(linewidth=1, color='r')
                ax6.set_xlabel('Time')
                ax6.set_ylabel('Residuals (data-model)')            
                """
        pass


            # Plot PHOT in time, for spots, planetary systems
        
            # Get flux array and normalize it
            if not mergeLC:
                flux = datadictLC[dkey]['data']['flux']/foot
                sflux = datadictLC[dkey]['data']['sflux']/foot
            
                #if correct_contam:
                flux = flux/(1.0 - cont) - cont/(1.0 - cont)
                sflux = sflux/(1.0 - cont)
            else:
                flux = datadictLC[dkey]['data']['flux']
                sflux = datadictLC[dkey]['data']['sflux']

            if lcfmt['modelontop']:
                zmod = 2
                zdat = 1
            else:
                zmod = 1
                zdat = 2

            f6 = p.figure(figsize = fsize)
            figlist.append(f6)
            ax11 = f6.add_axes(axsize)
            ax12 = f6.add_axes(ressize)
            
            #ax11.errorbar(lctime, flux, sflux, fmt = lcfmt['fmt'], mec = lcfmt['mec'],
        #         capsize = lcfmt['capsize'], zorder = zdat)
            ax11.plot(lctime, flux, lcfmt['fmt'],zorder = zmod)
            ax11.plot(lctimetAll, ftAll, lcfmt['lfmt'], lw = lcfmt['lw'], zorder = zmod)
           
            ax12.errorbar(lctime, (flux - ft)*1e6, 1e6*sflux, fmt = lcfmt['fmt'], 
                          zorder = zdat, capsize = lcfmt['capsize'])
            
            ax11.xaxis.set_major_formatter(p.NullFormatter())
            ax11.yaxis.set_major_formatter(p.ScalarFormatter(useOffset = False))
            ax12.yaxis.set_major_formatter(p.ScalarFormatter(useOffset = False))
                        
            ax12.axhline(linewidth=1, color=lcfmt['reshlcol'],zorder = zmod)
            ax12.set_ylabel('O-C [ppm]')
            ax12.yaxis.set_major_locator(p.MultipleLocator(lcfmt['resmajor']))
            ax12.yaxis.set_minor_locator(p.MultipleLocator(lcfmt['resminor']))

            
            ax12.set_xlim(lctimetAll.min(), lctimetAll.max())
            ax12.set_xlabel('Time', size = lcfmt['labelsize'])
        
            ax11.xaxis.set_major_formatter(p.NullFormatter())
            ax11.set_ylabel('Relative flux', size = lcfmt['labelsize'])
            ax11.set_title(dkey)
            ax11.set_xlim(lctimetAll.min(), lctimetAll.max())

            #Phase plots
        for i in range(len(periods)):
        if not datadictLC[dkey]['is_phase']:
                    print('ploting phase for P = %f, T0 = %f'%(periods[i],tps[i]))
            phase = (lctime - tps[i])/periods[i]%1.0
            phaset = (lctimeb - tps[i])/periods[i]%1.0
        else:
            phase = lctime
            phaset = lctimeb
                
        # Create figure and axes or used exisiting ones if merge is
        # True
        if mergeLC and dkey in axdict.keys() and False:
            ax7 = axdict[dkey][0]
            ax8 = axdict[dkey][1]

                    if lcbin and lcbinnewfig:
                        axb = axdict[dkey][2]
                        axbres = axdict[dkey][3]
        else:
            f4 = p.figure(figsize = fsize)
                    figlist.append(f4)

                    ax7 = f4.add_axes(axsize)
            #Residuals
            ax8 = f4.add_axes(ressize)
                    
            # Add axes to dictionary
            axdict[dkey] = []
            axdict[dkey].append(ax7)
            axdict[dkey].append(ax8)

                    if lcbin and lcbinnewfig:
                        f5 = p.figure(figsize = fsize)
                        figlist.append(f5)
                
                        axb = f5.add_axes(axsize)
                #Residuals
                        axbres = f5.add_axes(ressize)
                        
                        axdict[dkey].append(axb)
                        axdict[dkey].append(axbres)

        ax7.xaxis.set_major_formatter(p.NullFormatter())
                ax7.yaxis.set_major_formatter(p.ScalarFormatter(useOffset = False))
                ax8.yaxis.set_major_formatter(p.ScalarFormatter(useOffset = False))

                if plotsecondary:
                    ph = phase
                    pht = phaset
                else:
                    ph = n.where(phase > 0.5, phase - 1., phase)
                    pht = n.where(phaset > 0.5, phaset - 1., phaset)

                '''
        # Get flux array and normalize it
        flux = datadictLC[dkey]['data']['flux']/foot
        sflux = datadictLC[dkey]['data']['sflux']/foot

        if correct_contam:
            flux = flux/(1.0 - cont) - cont/(1.0 - cont)
            sflux = sflux/(1.0 - cont)
                '''

        fluxt = ft

        ind = n.argsort(ph)
        flux2 = flux[ind]; sflux2 = sflux[ind]; ph2 = ph[ind]
        indt = n.argsort(pht)
        fluxt2 = fluxt[indt]; pht2 = pht[indt]


        if lcfmt['modelontop']:
            zmod = 2
            zdat = 1
        else:
            zmod = 1
            zdat = 2

                ## Define array "hours from transit"
                h = ph*periods[i]*24.0

        if usehours:
            #if unbinned: ax7.errorbar(h, flux, sflux, fmt = lcfmt['fmt'], zorder = zdat)
            if unbinned: ax7.plot(h, flux2, lcfmt['fmt'], zorder = zdat)

            # Plot theoretical curve
                    ax7.plot(h, fluxt2, lcfmt['lfmt'], lw = lcfmt['lw'],
                             zorder = zmod)
            
            ax7.set_xlim(-4.0, 4.0,)

            if unbinned: ax8.errorbar(h, (flux2 - fluxt2)*1e6, 1e6*sflux2,
                 fmt = lcfmt['fmt'], zorder = zdat)
            ax8.set_xlim(-4.0, 4.0)
                    
            ax8.set_xlabel('Time from midtransit [hours]',
                                   size = lcfmt['labelsize'])


        else:
            #if unbinned: ax7.errorbar(ph2, flux2, sflux2, fmt=',',color='k')
            if unbinned: ax7.plot(ph2, flux2, lcfmt['fmt'], zorder = zdat)

                    # Plot theoretical curve
                    ax7.plot(pht2, fluxt2, lcfmt['lfmt'], lw = lcfmt['lw'],
                             zorder = zmod)

                    '''
                    f=open('res_norm.pickle','w')
                    res_norm=(flux - fluxt)/sflux
                    pickle.dump(res_norm,f)
                    f.close()
                    '''

            if unbinned: ax8.errorbar(ph2, (flux2 - fluxt2)*1e6, 1e6*sflux2,
                 fmt = lcfmt['fmt'], zorder = zdat)
                        
            #ax8.plot(ph, (flux - fluxt)*1e6, '.k')
            ax8.set_xlim(ph2.min(), ph2.max())
            ax8.set_xlabel('Orbital phase', size = lcfmt['labelsize'])


        if lcbin:
                    # Plot binned points
            bins = n.arange(n.min(ph2), n.max(ph2) + lcbinsize, lcbinsize)

            phb = [] #n.zeros(bins.size)
            fluxb = []#n.zeros(bins.size)
            sfluxb = []#n.zeros(bins.size)
                    fluxtb = []
            ind = n.digitize(ph2, bins)
            for ii in range(1, len(bins) + 1):
            condi = (ind == ii)
            if not n.any(condi): continue
            wi = 1.0/(n.compress(condi, sflux2))**2
            phb.append(n.average(n.compress(condi, ph2),
                         weights = wi)
                   )
            fluxb.append(n.average(n.compress(condi, flux2),
                           weights = wi)
                     )
            sfluxb.append(n.sqrt(1.0/n.sum(wi))
                                      )
                        fluxtb.append(n.average(n.compress(condi, fluxt2))
                                      )
            """
            phb[i - 1] = n.average(n.compress(condi, ph),
                           weights = wi)
            fluxb[i - 1] = n.average(n.compress(condi, flux),
                         weights = wi)
            sfluxb[i - 1] = sqrt(1.0/n.sum(wi))
            """
            phb, fluxb, sfluxb, fluxtb = map(n.array, (phb, fluxb, sfluxb,
                                                               fluxtb))
                    phb = phb[:-1]
                    fluxb =fluxb[:-1]
                    sfluxb = sfluxb[:-1]
                    fluxtb = fluxtb[:-1]

                    if len(phb) == 0:
                        raise Exception('Binning size too big!, no binned points available.')

            if lcbinnewfig:

                        axb = axdict[dkey][2]
                        axbres = axdict[dkey][3]

                        axb.set_title(dkey)
                        axb.set_ylabel('Relative flux')
                        axbres.set_xlabel('Orbital phase')
                        axb.xaxis.set_major_formatter(p.NullFormatter())

                        yformatter = p.ScalarFormatter(useOffset = False)
                        axb.yaxis.set_major_formatter(yformatter)
                        axbres.yaxis.set_major_formatter(yformatter)
                        
                        axbres.axhline(linewidth=1, color=lcfmt['reshlcol'],
                                       zorder = zmod)
                        axbres.set_ylabel('O-C [ppm]')
                        #axbres.yaxis.set_major_locator(p.MultipleLocator(lcfmt['resmajor']))
                        axbres.yaxis.set_minor_locator(p.MultipleLocator(lcfmt['resminor']))

                        # Plot theoretical curve
                        axb.plot(pht2, fluxt2, lcfmt['lfmt'], lw = lcfmt['lw'],
                                 zorder = zmod)
                    
                    else:
                        axb = axdict[dkey][0]
                        axbres = axdict[dkey][1]
                       
                    axb.errorbar(phb, fluxb, sfluxb, fmt = lcfmt['binfmt'],
                 ms = lcfmt['binmarkersize'],
                                 mfc = lcfmt['binmfc'], zorder = zdat
                                 )
                    axbres.errorbar(phb, 1e6*(fluxb - fluxtb), 1e6*sfluxb, 
                                    fmt = lcfmt['binfmt'],
                                    ms = lcfmt['binmarkersize'],
                                    mfc = lcfmt['binmfc'], zorder = zdat
                                    )

                    if rangephLC == None:
                        axb.set_xlim(ph.min(), ph.max())
                        axbres.set_xlim(ph.min(), ph.max())
                    else:
                        axb.set_xlim(rangephLC[0],rangephLC[1])
                        axbres.set_xlim(rangephLC[0],rangephLC[1])

        ax7.xaxis.set_major_formatter(p.NullFormatter())
        ax7.set_ylabel('Relative flux', size = lcfmt['labelsize'])
                ax7.set_title(dkey)

#ax7.set_title('P ~ %.1f days'%periods[phases])
        ax7.set_xlim(ph.min(), ph.max())

        ax8.axhline(linewidth=1, color=lcfmt['reshlcol'], zorder = zmod)
        ax8.set_ylabel('O-C [ppm]', size = lcfmt['labelsize'])
        ax8.yaxis.set_major_locator(p.MultipleLocator(lcfmt['resmajor'])
                                            )
        ax8.yaxis.set_minor_locator(p.MultipleLocator(lcfmt['resminor'])
                                            )

                if rangephLC != None:
                    ax7.set_xlim(rangephLC[0],rangephLC[1])
                    ax8.set_xlim(rangephLC[0],rangephLC[1])

        #######
        # SED #
        #######
        if datadict[key]['type'] == 'SED':
            compute_global_spectrum(*objects)
            from ..models.SED import global_spectrum


            if sedfmt['modelontop']:
                zmod = 2
                zdat = 1
            else:
                zmod = 1
                zdat = 2
                    
            sedmags = datadict[key]['data']['mags']
            # Compute theoretical magnitudes and likelihood
            magst = PASTIS_SED(datadict[key]['data']['band'], *objects)
            photbands=datadict[key]['data']['band']
            magst = []
            for pb in photbands:
                magst.append(phot.flux2mag(phot.get_flux(global_spectrum,
                                                         pb), pb))

            nphotbands=len(photbands)
            xx = [] 
            yy = [] 
            eyy = []
            my = []
            for pb in range(nphotbands):
                xx.append(phot.ZeroMag[photbands[pb]]['lambda_eff'])
                yy.append(phot.mag2flux(sedmags[pb], photbands[pb]))
                argg = []
                for i in range(1000):
                    argg.append(phot.mag2flux((sedmags[pb]) + (datadict[key]['data']['smags'][pb])*n.random.randn(), photbands[pb]))
                eyy.append(n.std(argg))
                my.append(phot.get_flux(global_spectrum, photbands[pb]))


            f5 = p.figure(figsize = fsize)
            figlist.append(f5)

            ax9 = f5.add_axes(axsize)
            ax9.xaxis.set_major_formatter(p.NullFormatter())

            ax9.errorbar(xx, yy, eyy, fmt = sedfmt['fmt'], mec = sedfmt['mec'],
             capsize = sedfmt['capsize'], zorder = zdat)
            ax9.plot(phot.ww, global_spectrum, sedfmt['lfmt'],
             lw = sedfmt['lw'], zorder = zmod)
            ax9.plot(xx, my, sedfmt['binfmt'], ms = sedfmt['binmarkersize'],
                     mfc = sedfmt['binmfc'], zorder = zdat)

            #semilogx()
            ax9.loglog()
            ax9.set_xlim(n.min(xx) - 1e3, n.max(xx) + 1e4)
            ax9.set_ylim(0.8*n.min([yy,my]),1.2*n.max([yy,my]))
            ax9.set_ylabel('Flux [ergs cm$^{-2}$ s$^{-1}$ A$^{-1}$]',
                           size = sedfmt['labelsize'])
            ax9.xaxis.set_major_formatter(p.NullFormatter())
            #ax9.yaxis.set_major_formatter(p.NullFormatter())

            #Residuals
            ax10 = f5.add_axes(ressize)
            ax10.errorbar(xx, sedmags - magst, datadict[key]['data']['smags'],
              fmt = sedfmt['fmt'], mec = sedfmt['mec'],
              capsize = sedfmt['capsize'])
            ax10.axhline(0.0, ls = ':', color = '0.55', zorder = 0)
            
        ax10.semilogx()
        ax10.set_xlim(n.min(xx)-1e3,n.max(xx)+1e4)
            ax10.set_ylim(1.2*n.min(sedmags - magst - datadict[key]['data']['smags']),
              1.2*n.max(sedmags - magst + datadict[key]['data']['smags'])
              )
            ax10.set_xlabel('Wavelength [\AA]', size = sedfmt['labelsize'])
            ax10.set_ylabel('O-C [mag]', size = sedfmt['labelsize'])
        ax10.yaxis.set_major_locator(p.MultipleLocator(sedfmt['resmajor']))
        ax10.yaxis.set_minor_locator(p.MultipleLocator(sedfmt['resminor']))

    if 'global_spectrum' in globals():
        del  globals()['global_spectrum']

    ## Change size of all ticklabels
    # Ticklabels
    for ff in figlist:
        for aa in ff.axes:
            for ticklabel in aa.xaxis.get_majorticklabels():
                ticklabel.set_fontsize(ticksize)
            for ticklabel in aa.yaxis.get_majorticklabels():
                ticklabel.set_fontsize(ticksize)
    
    print('DONE')
    p.show()
    return



