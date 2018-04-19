import os
import numpy as n
from scipy import interpolate

from .paths import idlexec, libpath
from .tools import loadtxt_iter

#Define path
extpath = os.path.join(libpath, 'GIE')
runfilename = os.path.join(extpath, 'run.pro')
runfilename2 = os.path.join(extpath, 'run.pro')


def initialize_extinction(ra, dec, maxdist, extinction_step,
                       EXTpath = extpath,
                       Rv = 3.1):

    # Create file name based on input parameters
    EXTfile = os.path.join(EXTpath, 'ext_ra%.3f_dec%.3f_max%.d_step%.1f.txt'%\
                           (ra, dec, maxdist, extinction_step))

    # If file exists, do not compute the extinction again
    if not os.path.exists(EXTfile):
        ## Compute extinction using IDL
        print('Computing extintion coefficients...')

        f = open(runfilename, 'w')
        f.write('.r '+os.path.join(extpath, 'extinction.pro')+'\n')
        f.write('extinction, ra=%f, dec=%f, equinox = 2000.0, min = 0.0, max = %f, step = %f, outfile = \'%s\''%(ra, dec, maxdist, extinction_step, EXTfile))
        f.close()

        os.system(idlexec+' -quiet < '+runfilename+' > /dev/null')

    else:
        print('Extinction file %s exists. It will not be computed again.'%EXTfile)

    ### Read extinction model for a given coordinate
    dist, EBmV, Av = loadtxt_iter(EXTfile, unpack=True, skiprows=2)  # distance in pc

    # Create interpolated function
    global Avfunc, EBmVfunc
    Avfunc = interpolate.interp1d(dist, Av)
    EBmVfunc = interpolate.interp1d(dist, EBmV)

    print('... DONE! \n')
    return

def get_extinction(ra, dec, dist, EXTpath = extpath, Rv = 3.1):
    """
    Return E(B-V) from the Amores & Lepine - 2005, AJ, 130, 679 model 
    for a given cordinate, distance and R_V (relation between the 
    extinction in visual and the color excess E(B-V))

    Parameters
    ----------
    ra: float
        right ascension in degrees
    
    dec: float
        declination in degrees 

    dist: float
        distance in parsecs

    Other parameters
    ----------------
    EXTpath: str
        Path of the IDL programs

    Rv: float
        relation between the extinction in visual and the color excess E(B-V)

    """

    # Create file name based on input parameters
    EXTfile = os.path.join(EXTpath, 'ext_ra%.3f_dec%.3f_dist%.d.txt'%\
                           (ra, dec, dist))

    # If file exists, do not compute the extinction again
    if not os.path.exists(EXTfile):
        ## Compute extinction using IDL
        print('Computing extintion coefficients...')

        f = open(runfilename2, 'w')
        f.write('.r '+os.path.join(extpath, 'extinction.pro')+'\n')
        f.write('extinction2, ra=%f, dec=%f, equinox = 2000.0, dist=%f, outfile = \'%s\''%(ra, dec, dist, EXTfile))
        f.close()

        os.system(idlexec+' -quiet < '+runfilename+' > /dev/null')

    else:
        print('Extinction file %s exists. It will not be computed again.'%EXTfile)

    ### Read extinction model for a given coordinate
    dist, EBmV, Av = loadtxt_iter(EXTfile, unpack=True, skiprows=2)  # distance in pc

    print('... DONE! \n')
    return EBmV

def ext_law(llambda, Rv = 3.1):
    x = 1.0e4 /llambda
    global extlaw
    extlaw = x*0.0
    extlaw = extlaw.astype('float64')
    
    x0 = 4.596  
    gamma = 0.99	
    c3 = 3.23	
    c4 = 0.41    
    c2 = -0.824 + 4.717/Rv
    c1 = 2.030 - 3.007*c2
    
    
    # Compute UV portion of A(lambda)/E(B-V) curve using FM fitting function and 
    # R-dependent coefficients
    xcutuv = 10000.0/2700.0
    xspluv = 10000.0/n.array([2700.0,2600.0])
    
    iuv = n.where(x >= xcutuv)
    iopir = n.where(x < xcutuv)
    
    #iuv = where(x ge xcutuv, N_UV, complement = iopir, Ncomp = Nopir)

    if len(iuv) > 0:
        xuv = n.concatenate((xspluv, x[iuv]))
    else:
        xuv = xspluv

    yuv = c1  + c2*xuv
    yuv = yuv + c3*xuv**2/((xuv**2-x0**2)**2 +(xuv*gamma)**2)
    xuvv = n.where(xuv > 5.9, xuv, 5.9 + xuv*0.0)
    yuv = yuv + c4*(0.5392*(xuvv-5.9)**2+0.05644*(xuvv-5.9)**3)
    yuv = yuv + Rv
    yspluv  = yuv[0:2]

    if len(iuv) > 0:
        extlaw[iuv] = yuv[2:]

    # Compute optical portion of A(lambda)/E(B-V) curve
    # using cubic spline anchored in UV, optical, and IR

    xsplopir = n.concatenate(([0,],10000.0/n.array([26500.0,12200.0,6000.0,5470.0,4670.0,4110.0])))
    ysplir   = n.array([0.0,0.26469,0.82925])*Rv/3.1 
    ysplop   = [n.poly1d([2.13572e-04, 1.00270, -4.22809e-01 ])(Rv), n.poly1d([-7.35778e-05, 1.00216, -5.13540e-02])(Rv),
                n.poly1d([-3.32598e-05, 1.00184, 7.00127e-01])(Rv),
                n.poly1d([-4.45636e-05, 7.97809e-04, -5.46959e-03, 1.01707, 1.19456])(Rv)
                ]

    ysplopir = n.concatenate([ysplir, ysplop])
    
    if len(iopir) > 0:
        xt = n.concatenate([xsplopir,xspluv])
        yt = n.concatenate([ysplopir,yspluv])
        tck = interpolate.splrep(xt, yt, k = 3)
        extlaw[iopir] = interpolate.splev(x[iopir], tck)
    
    return extlaw





