import sys, os
import numpy as n
from scipy import interpolate
# import pyfits
if sys.version_info < (3, 0):
    import cPickle as pickle
else:
    import pickle


# Intra-package import
from .paths import libpath, filterpath, zeromagfile
from .constants import c
from .exceptions import SpectrumInterpolError
from .tools import loadtxt_iter
from . import tools

Allpbands = ['CoRoT-W', 'Kepler', 'IRAC-I1', 'IRAC-I2', 'IRAC-I3', 'IRAC-I4',
             'SDSS-U', 'SDSS-G', 'SDSS-R', 'SDSS-I', 'SDSS-Z', 'Johnson-U',
             'Johnson-B', 'Johnson-V', 'Johnson-R' ,'Johnson-I', '2MASS-J',
             '2MASS-H', '2MASS-Ks', 'STROMGREN-u', 'STROMGREN-v', 'STROMGREN-b',
             'STROMGREN-y', 'WISE-W1', 'WISE-W2', 'WISE-W3', 'WISE-W4',
             'Bessell-U', 'Bessell-B', 'Bessell-V', 'Bessell-R', 'Bessell-I',
             'MIPS-M1', 'Cousins-R', 'Cousins-I', 'TESS']

#, 'MIPS-M2', 'MIPS-M3'


def initialize_phot(pbands = Allpbands, ZEROMAGfile = zeromagfile,
		    FILTERdir = filterpath, AMmodel = 'BT'):

    ### Read AM spectra and put them in global variable
    global AMspectra, AMspectra01
    if AMmodel == 'BT':
        print('Reading BT spectra...')
        AMspectra, AMspectra01 = read_BTsettl_spectra()
    elif AMmodel == 'CK':
        print('Reading Castelli & Kuruckz spectra...')
        AMspectra, AMspectra01 = read_CK_spectra()
    print('... DONE! \n')
    ###

    print('Loading zero magnitude flux information...')
    global ZeroMag
    ZeroMag = read_zero_magnitude_flux(ZEROMAGfile)
    print('... DONE! \n')


    print('Loading filter transmission curves...')
    global Filters
    global NormalisedFilters 
    Filters = {}
    NormalisedFilters = {}

    for i, pband in enumerate(pbands):

        if i == 0:
            sys.stdout.write('... band: %-12s'%pband)
        else:
            sys.stdout.write('\b'*12+'%-12s'%pband)

        # sys.stdout.flush()
        Filters[pband] = interpol_filter(ww, FILTERdir, pband)
        NormalisedFilters[pband] = Filters[pband]/tools.area(ww, Filters[pband])

    print('... DONE! \n')
    return

def initialize_phot_WD():
    ### Read White Dwarf atmosphere models
    global WDspectra, WDspectra01
    print('Reading WD spectra...')
    WDspectra, WDspectra01 = read_WD_spectra()
    

def read_CK_spectra():
    """
    Read Castelli & Kurucz (2004) stellar atmospheric models.

    Flux units are in erg/s/cm^2/A.
    """
    global AMz
    AMz = n.array([-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.2, 0.5])
    
    global AMteff
    AMteff0 = n.arange(3500.0, 13000., 250.0)
    AMteff1 = n.arange(13000., 50000.1, 1000.0)
    AMteff = n.concatenate([AMteff0, AMteff1])

    global AMlogg
    AMlogg = n.arange(0.0, 5.5, 0.5)
    
    ## Create dictionary with information about the model
    global infoAM
    infoAM = {'lenAMz' : len(AMz),
              'lenAMteff' : len(AMteff),
              'lenAMlogg' : len(AMlogg),
              'minAMz' : n.min(AMz),
              'maxAMz' : n.max(AMz),
              'minAMteff' : n.min(AMteff),
              'maxAMteff' : n.max(AMteff),
              'minAMlogg' : n.min(AMlogg),
              'maxAMlogg' : n.max(AMlogg)
              }


    AMspectra = {}#n.zeros([len(AMz), len(AMteff), len(AMlogg), 1221])
    AMspectra01 = n.zeros([len(AMz), len(AMteff), len(AMlogg)])
    
    for i, zz in enumerate(AMz):

        if zz >= 0:
            metsign = 'p'
        else:
            metsign = 'm'

        parentdirectory = 'ck%s%02d'%(metsign, abs(zz)*10)

        for j, teff in enumerate(AMteff):

            filename = '%s_%d.fits'%(parentdirectory, teff)
            fitsfilename = os.path.join(libpath, 'AM', 'ck04models',
					parentdirectory,filename)

            fitsfile = pyfits.open(fitsfilename)
            table = fitsfile[1]

            # READ WAVELENGTH OF SPECTRA, make it globally available
            if i == 0 and j == 0:
                global ww
                ww = table.data.field(table.columns[0].name)

            # READ INDIVIDUAL SPECTRA; PUT THE IN CKspectra
            for k, logg in enumerate(AMlogg):
                spec = table.data.field('g%02d'%(logg*10))
                if n.max(spec) != 0:
                    AMspectra[(zz, teff, logg)] = spec/teff**4.0
                    AMspectra01[i, j, k] = 1.

            fitsfile.close()
                    
    return AMspectra, AMspectra01


def read_BTsettl_spectra():
    """
    Read BTsettl stellar atmospheric models.

    Flux units are in erg/s/cm^2/A.
    """
    # Read spectra from file
    f = open(os.path.join(libpath, 'AM', 'BT-Settl', 'BTspec_20.pickle'), 'rb')

    global ww, AMz, AMteff, AMlogg
    if sys.version_info > (3, 0):
        ww, AMz, AMteff, AMlogg, AMspectra = pickle.load(f, encoding='bytes')
    else:
        ww, AMz, AMteff, AMlogg, AMspectra = pickle.load(f)
    f.close()

    ## Create dictionary with information about the model
    global infoAM
    infoAM = {'lenAMz' : len(AMz),
              'lenAMteff' : len(AMteff),
              'lenAMlogg' : len(AMlogg),
              'minAMz' : n.min(AMz),
              'maxAMz' : n.max(AMz),
              'minAMteff' : n.min(AMteff),
              'maxAMteff' : n.max(AMteff),
              'minAMlogg' : n.min(AMlogg),
              'maxAMlogg' : n.max(AMlogg)
              }
    
    ## Normalize spectra by effective temperature
    for kk in AMspectra.keys():
        AMspectra[kk] = AMspectra[kk]/kk[1]**4.0
    
    AMspectra01 = n.zeros((len(AMz), len(AMteff), len(AMlogg)), 'int')

    # Fill in the gaps in BTsettl grid. See notes!
    AMspectra[(0.0, 700, 5.5)] = 0.5*(AMspectra[(0.0, 600, 5.5)] + AMspectra[(0.0, 800, 5.5)])
    AMspectra[(0.0, 1000, 3.5)] = 0.5*(AMspectra[(0.0, 900, 3.5)] + AMspectra[(0.0, 1100, 3.5)])
    
    for i in range(len(AMz)):
        for j in range(len(AMteff)):
            for k in range(len(AMlogg)):
                try:
                    if n.max(AMspectra[AMz[i], AMteff[j], AMlogg[k]]) > 0:
                        AMspectra01[i, j, k] = 1
                except:
                    pass
        
    return AMspectra, AMspectra01


def read_WD_spectra():
    """
    Read white dwarf atmospheric models.

    Flux units are in erg/s/cm^2/A.
    """
    # Read spectra from file
    f = open(os.path.join(libpath, 'AM', 'WD', 'WDspec.pickle'), 'rb')

    global WDteff, WDlogg
    WDww, WDteff, WDlogg, WDspectra = pickle.load(f)
    f.close()

    # interpolate spectra to ww and normalize it by effective temperature
    for kk in WDspectra.keys():
        # Prepare flux array
        ff = interpolate.interp1d(WDww, WDspectra[kk], bounds_error = False,
                                            fill_value = 0.0)(ww)

        WDspectra[kk] = ff/kk[0]**4.0
    
    WDspectra01 = n.zeros((len(AMteff), len(AMlogg)), 'int')
    
    for j in range(len(WDteff)):
        for k in range(len(WDlogg)):
            try:
                if n.max(WDspectra[WDteff[j], WDlogg[k]]) > 0:
                    WDspectra01[j, k] = 1
            except:
                pass
        
    return WDspectra, WDspectra01


def read_wwHR():
    """
    Read wavelenght of the high resolution BT-Settl atmosphere model 
    if wwHR is already defined, skip
    """
    if not 'wwHR' in globals():
        f = open(os.path.join(libpath, 'AM', 'BT-Settl', 'Flux_envolvente',
			      'lambda.pickle')
		 )
	
        global wwHR
        wwHR = pickle.load(f)
        f.close()

def read_HR_AM(z, teff, logg):
    """
    Read high resolution BT-Settl atmosphere model 
    only if is available for a given z, teff, logg
    """
    pathHRAM= os.path.join(libpath, 'AM', 'BT-Settl', 'Flux_envolvente')

    try:
        # Try open the file
        f = open(pathHRAM+'/BT'+str(round(z, 1))+'_'+str(int(teff))+'_'+str(round(logg, 1))+'.pickle')
        ffHR = pickle.load(f)
        f.close()
        read_wwHR()
        return ffHR
    except IOError:
        print('get_HR_AM: required spectra not available')

def get_AM(z, teff, logg, HR=False):
    try:
        return get_interpolated_AM(z, teff, logg, HR= HR)
    except SpectrumInterpolError as sie:
        return get_nearest_AM(sie.indz, sie.indteff, sie.indlogg, HR= HR,
                              outsidegrid = sie.outside,
                              isindex = True
                              )
                              
    """
    except SpectrumInterpolError as sie:
        return get_nearest_AM(sie.indz, sie.indteff, sie.logg, HR= HR,
                              outsidegrid = False, isindex = True)
    """
    

    
def get_interpolated_AM(z, teff, logg, HR = False):
    """
    Get interpolated spectrum from AMspectra
    """

    indz = n.searchsorted(AMz, z)
    indteff = n.searchsorted(AMteff, teff)
    indlogg = n.searchsorted(AMlogg, logg)
    
    if indz == 0 or indz == infoAM['lenAMz']:
        raise SpectrumInterpolError('Metallicity (z = %.2f) outside grid.'%z,
                                     indz, indteff, indlogg, outside = True
                                     )
    if indteff == 0 or indteff == infoAM['lenAMteff']:
        raise SpectrumInterpolError('Effective temperature (teff = %d) outside grid.'%teff,
                                     indz, indteff, indlogg, outside = True
                                     )
    if indlogg == 0 or indlogg == infoAM['lenAMlogg']:
        raise SpectrumInterpolError('Log(g) (%.2f) outside grid.'%logg,
                                    indz, indteff, indlogg, outside = True
                                    )

    # Old definitions
    # indz0 = indz - 1; indz1 = indz
    # indteff0 = indteff - 1; indteff1 = indteff
    # indlogg0 = indlogg - 1; indlogg1 = indlogg

    # List of spectra in nodes
    if HR:
        read_wwHR()
        spectra = n.zeros((2, 2, 2, len(wwHR)), 'double')
    else: 
        spectra = n.zeros((2, 2, 2, len(ww)), 'double')

    # Get the spectra from the nodes
    for i, iz in zip(range(2), (indz - 1, indz)):
        for j, iteff in zip(range(2), (indteff - 1, indteff)):
            for k, ilogg in zip(range(2), (indlogg - 1, indlogg)):
                if HR:
                    try:
                        spectra[i, j, k] = read_HR_AM(AMz[iz], AMteff[iteff], AMlogg[ilogg])/AMteff[iteff]**4
                    except TypeError:
                        raise SpectrumInterpolError('Missing node in grid.',
                                                    indz, indteff, indlogg
                                                    )
                    
                else:
                    try:
                        spectra[i, j, k] = AMspectra[(AMz[iz],
                                                  AMteff[iteff],
                                                  AMlogg[ilogg])]
                    except KeyError:
                        raise SpectrumInterpolError('Missing node in grid.',
                                                    indz, indteff, indlogg
                                                    )

    # If all 8 nodes exist, perform a trilinear interpolation
    indzz = (z - AMz[indz-1])/(AMz[indz] - AMz[indz-1])
    indtt = (teff - AMteff[indteff-1])/(AMteff[indteff] - AMteff[indteff-1])
    indgg = (logg - AMlogg[indlogg-1])/(AMlogg[indlogg] - AMlogg[indlogg-1])
    indices = n.array((indzz,indtt,indgg))

    # Get spectrum by trilinear interpolation
    spectrum = tools.trilinear_interpolation(spectra, indices) 
    return spectrum*teff**4.0

def get_nearest_AM(z, teff, logg, HR = False, outsidegrid = None,
                   isindex = False):
    """
    Get nearest spectrum from AMspectra

    Parameters
    ----------
    z, teff, logg: float or int
        Values used as input. They can be either the parameter values or the
        indices obtained using searchsorted in get_interpolated_AM. In this
        case, isindex must be True.
    """
    if isindex:
        indz = z
        indteff = teff
        indlogg = logg

    else:
        indz = n.searchsorted(AMz, z)
        indteff = n.searchsorted(AMteff, teff)
        indlogg = n.searchsorted(AMlogg, logg)
        

    if outsidegrid == None:
        # Check if outside grid
        if indz == 0 or indz == infoAM['lenAMz']:
            outsidegrid = True
        if indteff == 0 or indteff == infoAM['lenAMteff']:
            outsidegrid = True
        if indlogg == 0 or indlogg == infoAM['lenAMlogg']:
            outsidegrid = True

    ## Get index as float
    indz = (z - AMz[indz - 1])/(AMz[indz] - AMz[indz - 1])
    indteff = (teff - AMteff[indteff - 1])/(AMteff[indteff] - AMteff[indteff - 1])
    indlogg = (logg - AMlogg[indlogg - 1])/(AMlogg[indlogg] - AMlogg[indlogg - 1])

    ## Correct if needed
    if outsidegrid:
        ## Error in interpolation. Point is outside of grid
        if z < min(AMz):
            print('get_nearest_AM: Metallicity (z = %.2f) outside grid, using minimum value (%.2f).'%(z, min(AMz)))
            indz = 0

        elif z > max(AMz):
            print('get_nearest_AM: Metallicity (z = %.2f) outside grid, using maximum value (%.2f).'%(z, max(AMz)))
            indz = max(iz)

        if teff < min(AMteff):
            print('get_nearest_AM: Effective temperature (teff = %d) outside grid, using minimum value (%d).'%(teff, min(AMteff)))
            indteff = 0

        elif teff > max(AMteff):
            print('get_nearest_AM: Effective temperature (teff = %d) outside grid, using maximum value (%d).'%(teff, max(AMteff)))
            indteff = max(iteff)

        if logg < min(AMlogg):
             print('get_nearest_AM: Surface gravity (logg = %.2f) outside grid, using minimum value (%.2f).'%(logg, min(AMlogg)))
             indlogg = 0

        elif logg > max(AMlogg):
             print('get_nearest_AM: Surface gravity (logg = %.2f) outside grid, using maximum value (%.2f).'%(logg, max(AMlogg)))
             indlogg = max(ilogg)

    dist = AMspectra01*0+1e200
    for i in range(0, infoAM['lenAMz']):
        for j in range(0, infoAM['lenAMteff']):
            for k in range(0, infoAM['lenAMlogg']):
                if AMspectra01[i, j, k] == 1:
                    dist[i, j, k]=n.sqrt((i-indz)**2+(j-indteff)**2+(k-indlogg)**2)
                        
    indimin = n.argwhere(dist == n.min(dist))[0]
    #indimin = n.argmin(dist)

    zmin = AMz[indimin[0]]
    teffmin = AMteff[indimin[1]]
    loggmin = AMlogg[indimin[2]]

    if HR:
        spectrum = read_HR_AM(zmin, teffmin, loggmin)/teffmin**4
    else:
        spectrum = AMspectra[(zmin, teffmin, loggmin)]

    return spectrum*teff**4.0


def read_zero_magnitude_flux(ZEROMAGfile):
    """
    Read zero magnitude flux for all photometric bands.
    """
    f = open(ZEROMAGfile, 'r')
    lines = f.readlines()
    f.close()

    keys = lines[0].split()[1:-2]
    
    ZeroMag = {}
    for line in lines[2:]:
        band_dict = {}
        band_key = line.split()[0].rstrip()
        #
        ll = line.split()[1:-2]
        for i, key in enumerate(keys):
            band_dict[key] = float(ll[i])

        ZeroMag[line.split()[0].rstrip()] = band_dict

    return ZeroMag


def interpol_filter(ww, FILTERdir, photband = 'CoRoT'):
    """
    Interpolates the filter photband at wavelenghts contained in ww.
    """
    if photband == 'Kepler':
        x, y = loadtxt_iter(os.path.join(FILTERdir, photband+'.dat'), unpack = True, skiprows = 1)
        x = x*10.0

    else:
        x, y = loadtxt_iter(os.path.join(FILTERdir, photband+'.dat'), unpack = True)
  
    # Interpolate filter
    trans = interpolate.interp1d(x, y, bounds_error = False,
                                 fill_value = 0.0)(ww)

    return trans


## Other functions using ZeroMag and Filters
def mag2flux(mag, photband):
    """
    Flux in ergs/cm^2/s/A
    """
    f0 = ZeroMag[photband]['F0'] # Jy
    lambda_eff = ZeroMag[photband]['lambda_eff'] # in Armgstrong
    return (f0*10**(-0.4*mag)*(c*1e10)/lambda_eff**2)/1e23


def get_flux(spectrum, photband):
    """
    Integrate spectrum over normalised photband to get flux in units ergs/cm^2/s/A
    """
    # Recover filter from dictionary and integrate
    return tools.area(ww, NormalisedFilters[photband]*spectrum)


def flux2mag(flux, photband):
    """
    Convert flux density in ergs/cm^2/s/A to magnitude in band photband.
    """
    lambda_eff = ZeroMag[photband]['lambda_eff'] # in Armgstrong
    nu_eff = c*1e10/lambda_eff # in Hz

    flux_nu = flux*c*1e10/nu_eff**2 # in ergs/cm^2/s/Hz
    mag = -2.5*n.log10(1e23*flux_nu/ZeroMag[photband]['F0'])
    return mag

def corot_colors(Rflux, Gflux, Bflux, Rcont, Gcont, Bcont):
    """
    Compute corot RGB filters and limb-darkening coefficients 
    for a given total spectrum and relative flux in each color.
    It's assume that all the flux is inside the mask. 
    Total spectrum: ww,global_spectrum
    Rflux, Gflux, Bflux
    Rcont, Gcont, Bcont
    """
    from models import SED

    # Compute filters

    #Compute total spectrum observed by CoRoT
    CoRoT_W=Filters['CoRoT-W']
    global_spectrum_corot = SED.global_spectrum*CoRoT_W 

    # Decontaminate relative color flux
    Rflux_decont = Rflux*(1.0-Rcont)
    Gflux_decont = Gflux*(1.0-Gcont)
    Bflux_decont = Bflux*(1.0-Bcont)
    total_flux = Rflux_decont+Gflux_decont+Bflux_decont

    R_rf=Rflux_decont/total_flux
    G_rf=Gflux_decont/total_flux
    B_rf=Bflux_decont/total_flux
     
    # Cumultaive total spectrum
    dw = ww[1:] - ww[:-1]
    cumulative = n.concatenate((n.array([0]), n.cumsum((global_spectrum_corot[:-1] + global_spectrum_corot[1:])*0.5*dw)))
    cumulative_norm=cumulative/max(cumulative)
    
    # Find separation wavelenght between filters 
    count=0
    while cumulative_norm[count] == 0: count=count+1.
    ww_1 = ww[count]
    ww_2 = interpolate.interp1d(cumulative_norm, ww)(B_rf)
    ww_3 = interpolate.interp1d(cumulative_norm, ww)(B_rf+G_rf)
    ww_4 = n.array(ww[n.where(cumulative_norm > 0.9999)[0][0]])
    
    # test
    #print ' RED = ',ww_3,ww_4
    #print ' GREEN = ',ww_2,ww_3
    #print ' BLUE = ',ww_1,ww_2

    # Divide CoRoT transmision in 3 filters

    CoRoT_R = n.where(n.logical_and(n.greater(ww, ww_3), n.less_equal(ww, ww_4)),
                      CoRoT_W, 0.0)
    
    CoRoT_G = n.where(n.logical_and(n.greater(ww, ww_2), n.less_equal(ww, ww_3)),
                      CoRoT_W, 0.0)

    CoRoT_B = n.where(n.logical_and(n.greater(ww, ww_1), n.less_equal(ww, ww_2)),
                      CoRoT_W, 0.0)

    Filters['CoRoT-R'] = CoRoT_R
    Filters['CoRoT-G'] = CoRoT_G
    Filters['CoRoT-B'] = CoRoT_B
    NormalisedFilters['CoRoT-R'] = CoRoT_R/tools.area(ww, CoRoT_R)
    NormalisedFilters['CoRoT-G'] = CoRoT_G/tools.area(ww, CoRoT_G)
    NormalisedFilters['CoRoT-B'] = CoRoT_B/tools.area(ww, CoRoT_B)

    # Calculate transformation coefficients between sloan limbdarkening 
    # coefficients and CoRoT colors

    Sloan_u=Filters['SDSS-U']
    Sloan_g=Filters['SDSS-G']
    Sloan_r=Filters['SDSS-R']
    Sloan_i=Filters['SDSS-I']
    Sloan_z=Filters['SDSS-Z']
    
    Sloan_u2=Sloan_u
    Sloan_u2[n.where(Sloan_u2 > 0)]=1.0
    Sloan_g2=Sloan_g
    Sloan_g2[n.where(Sloan_g2 > 0)]=1.0
    Sloan_r2=Sloan_r
    Sloan_r2[n.where(Sloan_r2 > 0)]=1.0
    Sloan_i2=Sloan_i
    Sloan_i2[n.where(Sloan_i2 > 0)]=1.0
    Sloan_z2=Sloan_z
    Sloan_z2[n.where(Sloan_z2 > 0)]=1.0

    global CoRoT_LDC_weights
    CoRoT_LDC_weights=n.zeros([3,5], float)
    for color in ['R','G','B']:
        if color == 'R': 
            color_value=0
            ww_i=ww_3
            ww_e=ww_4
        if color == 'G': 
            color_value=1
            ww_i=ww_2
            ww_e=ww_3
        if color == 'B': 
            color_value=2
            ww_i=ww_1
            ww_e=ww_2
        indl=n.where(n.logical_and(ww >= ww_i,ww <= ww_e))
                
        # first weight have into acount the filter coverage
        wu1=tools.area(ww[indl],Sloan_u[indl])/tools.area(ww,Sloan_u)
        wg1=tools.area(ww[indl],Sloan_g[indl])/tools.area(ww,Sloan_g)
        wr1=tools.area(ww[indl],Sloan_r[indl])/tools.area(ww,Sloan_r)
        wi1=tools.area(ww[indl],Sloan_i[indl])/tools.area(ww,Sloan_i)
        wz1=tools.area(ww[indl],Sloan_z[indl])/tools.area(ww,Sloan_z)

        # second weight have into acount the shape and the emision espectrum 
        wu2=tools.area(ww[indl],global_spectrum_corot[indl]*Sloan_u[indl])/tools.area(ww[indl],global_spectrum_corot[indl]*Sloan_u2[indl])
            
        if n.isnan(wu2) : wu2 = 0.0

        wg2=tools.area(ww[indl],global_spectrum_corot[indl]*Sloan_g[indl])/tools.area(ww[indl],global_spectrum_corot[indl]*Sloan_g2[indl])
        if n.isnan(wg2) : wg2 = 0.0

        wr2=tools.area(ww[indl],global_spectrum_corot[indl]*Sloan_r[indl])/tools.area(ww[indl],global_spectrum_corot[indl]*Sloan_r2[indl])
        if n.isnan(wr2) : wr2 = 0.0

        wi2=tools.area(ww[indl],global_spectrum_corot[indl]*Sloan_i[indl])/tools.area(ww[indl],global_spectrum_corot[indl]*Sloan_i2[indl])
        if n.isnan(wi2) : wi2 = 0.0

        wz2=tools.area(ww[indl],global_spectrum_corot[indl]*Sloan_z[indl])/tools.area(ww[indl],global_spectrum_corot[indl]*Sloan_z2[indl])
        if n.isnan(wz2) : wz2 = 0.0

        # final weight 
        wu=wu1*wu2                 
        wg=wg1*wg2
        wr=wr1*wr2
        wi=wi1*wi2
        wz=wz1*wz2
        
        rwu=wu/(wu+wg+wr+wi+wz)
        rwg=wg/(wu+wg+wr+wi+wz)
        rwr=wr/(wu+wg+wr+wi+wz)
        rwi=wi/(wu+wg+wr+wi+wz)
        rwz=wz/(wu+wg+wr+wi+wz)

        # test
        #print 'weight1=',wu1,wg1,wr1,wi1,wz1
        #print 'weight2=',wu2,wg2,wr2,wi2,wz2
        #print 'weight=',wu,wg,wr,wi,wz
        #print 'relative weight=',rwu,rwg,rwr,rwi,rwz
        #pdb.set_trace()       

        CoRoT_LDC_weights[color_value,0]=rwu
        CoRoT_LDC_weights[color_value,1]=rwg
        CoRoT_LDC_weights[color_value,2]=rwr
        CoRoT_LDC_weights[color_value,3]=rwi
        CoRoT_LDC_weights[color_value,4]=rwz

    return


        

def get_interpolated_WD(teff, logg):
    """
    Get interpolated spectrum from WDspectra
    """
    iteff = range(len(WDteff))
    ilogg = range(len(WDlogg))

    try:
        indteff = interpolate.interp1d(WDteff, iteff)(teff)
    except ValueError:
        raise SpectrumInterpolError('Effective temperature (teff = %d) outside grid.'%teff)
    try:
        indlogg = interpolate.interp1d(WDlogg, ilogg)(logg)
    except ValueError:
        raise SpectrumInterpolError('Log(g) (%.2f) outside grid.'%logg)

    indteff0 = int(indteff); indteff1 = int(indteff + 1.0)
    indlogg0 = int(indlogg); indlogg1 = int(indlogg + 1.0)

    # List of spectra in nodes
    spectra = n.zeros((2, 2, len(ww)), 'double')

    # Get the spectra from the nodes
    for j, iteff in zip(range(2), (indteff0, indteff1)):
        for k, ilogg in zip(range(2), (indlogg0, indlogg1)):
            try:
                spectra[j, k] = WDspectra[(WDteff[iteff],WDlogg[ilogg])]
            except KeyError:
                import pdb
                pdb.set_trace()
                raise SpectrumInterpolError('Missing node in grid.')

    # If all 4 nodes exist, perform a trilinear interpolation
    indices = n.array((indteff - int(indteff),indlogg -int(indlogg)))

    # Get spectrum by trilinear interpolation
    spectrum = tools.bilinear_interpolation(spectra, indices)
    return spectrum*teff**4.0

def get_nearest_WD(teff, logg):
    """
    Get nearest spectrum from WDspectra
    """
    iteff = range(len(WDteff))
    ilogg = range(len(WDlogg))
    
    indcheck=False
    try:
        indteff = interpolate.interp1d(WDteff, iteff)(teff)
    except ValueError:
        indcheck=True
    try:
        indlogg = interpolate.interp1d(WDlogg, ilogg)(logg)
    except ValueError:   
        indcheck=True
        
    if indcheck:
        ## Error in interpolation. Point is outside of grid
        if teff < min(WDteff):
            print('get_nearest_WD: Effective temperature (teff = %d) outside grid, using minimum value (%d).'%(teff, min(WDteff)))
            indteff = 0

        elif teff > max(WDteff):
            print('get_nearest_WD: Effective temperature (teff = %d) outside grid, using maximum value (%d).'%(teff, max(WDteff)))
            indteff = max(iteff)

        if logg < min(WDlogg):
             print('get_nearest_WD: Surface gravity (logg = %.2f) outside grid, using minimum value (%.2f).'%(logg, min(WDlogg)))
             indlogg = 0

        elif logg > max(WDlogg):
             print('get_nearest_WD: Surface gravity (logg = %.2f) outside grid, using maximum value (%.2f).'%(logg, max(WDlogg)))
             indlogg = max(ilogg)

    dist = WDspectra01*0+1e200
    for j in range(0, len(WDteff)):
        for k in range(0, len(WDlogg)):
            if WDspectra01[j, k] == 1:
                dist[j, k]=n.sqrt((j-indteff)**2+(k-indlogg)**2)
                        
    indimin = n.argwhere(dist == n.min(dist))[0]

    teffmin = WDteff[indimin[0]]
    loggmin = WDlogg[indimin[1]]

    spectrum = WDspectra[(teffmin, loggmin)]

    return spectrum*teff**4.0

def get_beaming(sp, photband):
    """
    Compute B factor used in relativistic beaming effect.
    See Bloemen (MNRAS, 2011)
    """
    b = 5.+ tools.derivative(n.log(ww),n.log(sp)) #d ln F / d ln lambda
    return tools.area(ww, NormalisedFilters[photband]*ww*sp*b)/tools.area(ww,NormalisedFilters[photband]*sp*ww)
                                                                           
