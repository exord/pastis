# 2012-03-20: implemented CoRoT colors.
import numpy as n
#import pymacula

from .. import photometry as phot
#import task2_components as task2
from ..exceptions import EBOPparamError

def PASTIS_PHOT(t, photband, isphase, cont, foot, *args):
    """
    Compute SED for all objects in args, for all photbands.
    """
    from .. import AstroClasses as ac
    
    fluxes = []
    lightcurves = []
    for obj in args:
        if isinstance(obj, ac.Star):
            # Get flux from star
            fluxes.append(phot.get_flux(obj.get_spectrum(), photband))
            lightcurves.append(obj.get_spots(t,photband))
            
        elif isinstance(obj, ac.IsoBinary):
            # Add spectrum of each component, and compute flux
            spectrum = obj.star1.get_spectrum() + obj.star2.get_spectrum()
            fluxes.append(phot.get_flux(spectrum, photband))

            # Get lightcurve from binary
            lightcurves.append(obj.get_LC(t, photband, isphase))

        elif isinstance(obj, ac.PlanSys):
            # Get flux from star
            fluxes.append(phot.get_flux(obj.star.get_spectrum(), photband))

            # Get lightcurve from planetary system
            lightcurves.append(obj.get_LC(t, photband, isphase))
            
        elif isinstance(obj, ac.FitBinary):
            # Get flux from star
            fluxes.append(1.)

            # Get lightcurve from planetary system
            lightcurves.append(obj.get_LC(t, photband, isphase))

        elif isinstance(obj, ac.Triple):
            # Get LC for each component of triple system
            for  component in (obj.object1, obj.object2):
                if isinstance(component, ac.Star):
                    fluxes.append(phot.get_flux(component.get_spectrum(), photband))
                    lightcurves.append(n.ones(len(t))) # TEMP, implement spots

                elif isinstance(component, ac.IsoBinary):
                    # Add spectrum of each component, and compute flux
                    spectrum = component.star1.get_spectrum() + \
                               component.star2.get_spectrum()
                    
                    fluxes.append(phot.get_flux(spectrum, photband))

                    # Get lightcurve from binary
                    lightcurves.append(component.get_LC(t, photband, isphase))


                elif isinstance(component, ac.PlanSys):
                    # Get flux from star
                    fluxes.append(phot.get_flux(component.star.get_spectrum(), photband))

                    # Get lightcurve from planetary system
                    lightcurves.append(component.get_LC(t, photband, isphase))


    # Produce contaminating lightcurve and flux
    F = n.sum(fluxes)
    fluxes.append(F*cont/(1 - cont))
    lightcurves.append(n.ones(len(t)))
    
    # Generate global lightcurve: 
    fluxes = n.array(fluxes).reshape((len(fluxes),1))
    lightcurves = n.array(lightcurves)
    #import pdb
    #pdb.set_trace()
    #return fluxes, lightcurves
    #print fluxes
    return foot*n.sum(fluxes * lightcurves, axis = 0)/n.sum(fluxes)
    

def run_EBOP(v, ldtype, x, Nx, Nmax = None, components = False):
    """
    Runs JKTEBOP with variables v
    x in phase
    """
    print_warning = False
    ## Check limits of input parameters
    if v[2-1] > 0.8:
        print_warning = True; wstring = ['Sum of fractional radii', 'maximum', 0.8]
    if v[2-1] < 0.0:
        print_warning = True; wstring = ['Sum of fractional radii', 'minimum', 0.0]
    if v[3-1] < 0.0:
        print_warning = True; wstring = ['Radii ratio', 'minimum', 0.0]

    if v[3-1] > 100:
        print_warning = True; wstring = ['Radii ratio', 'maximum', 100.0]

    if v[6-1] < 50.0:
        print_warning = True; wstring = ['Inclination', 'minimum', 50.0]
                                
    if v[6-1] > 140.0:
        print_warning = True; wstring = ['Inclination', 'maximum', 140.0]

    if v[13-1] < 0.0:
        print_warning = True; wstring = ['Mass ratio', 'minimum', 0.0]

    if v[13-1] > 1e3:
        print_warning = True; wstring = ['Mass ratio', 'maximum', 1000.0]

    if v[7-1] < -1.0 or v[7-1] > 1.0:
        print_warning = True; wstring = ['Eccentricity', 'maximum', 1.0]

    if v[8-1] < -1.0 or v[8-1] > 1.0:
        print_warning = True; wstring = ['Eccentricity', 'maximum', 1.0]

    if v[9-1] < -10.0 or v[10-1] < -10:
        print_warning = True; wstring = ['Gravity darkening', 'minimum', -10.0]
    if v[9-1] > 10.0 or v[10-1] > 10:
        print_warning = True; wstring = ['Gravity darkening', 'maximum', 10.0]
    if v[1-1] > 100.0:
        print_warning = True; wstring = ['Surface brightness ratio', 'maximum', 100.0]
    if v[1-1] < 0.0:
        print_warning = True; wstring = ['Surface brightness ratio', 'minimum', 0.0]
    if v[15-1] < 0.0:
        print_warning = True; wstring = ['Third light', 'minimum', 0.0]
    if v[15-1] > 1.0:
        print_warning = True; wstring = ['Third light', 'maximum', 1.0]

    ### Avoid periapsis distance shorter than stellar radius!
    if v[2-1] > (1 - n.sqrt(v[7-1]**2 + v[8-1]**2)):
        raise EBOPparamError('Periapsis distance shorter than sum of radii.')
        
    if print_warning:
        raise EBOPparamError('%s out of bounds; %s allowed is %f'%(wstring[0], wstring[1], wstring[2]))

    if Nmax == None or Nx > Nmax:
        condg05 = x > 0.5
        condl05 = x < 0.5
        if ~condg05.all() or ~condl05.all():
            # There are not phases > 0.5 or phases < 0.5
            phmin = n.min(x)
            phmax = n.max(x)
        else:
            phmin = n.min(x[condg05] - 1.0)
            phmax = n.max(x[condl05])

        phcov = phmax - phmin
        
    if Nmax == None:
        # Compute apropiate Nmax
        if v[3-1] >= 0.0:
            R1 = v[2-1] / (1.0 + v[3-1])
            R2 = v[2-1] / (1.0 + (1.0/v[3-1]))
        else:
            R1 = v[2-1]
            R2 = n.abs(v[3-1])
            
        # from jktebop
        if R1 < 0.01 or R2 < 0.01:
            Nmax_allph = 1e5 
        else:
            Nmax_allph = 1e4

        # correct for phase coverage
        Nmax = n.int(Nmax_allph*phcov)

    if Nx > Nmax:
        N1 = round(Nmax*abs(phmax)/phcov)
        xx1 = n.linspace(0.0, phmax, N1)
        xx2 = n.linspace(1 + phmin, 1.0, Nmax - N1)
        xx = n.concatenate((xx1, xx2))
        
        if components:
            LPi,LSi,ECLi,REFLi = call_task2(v, ldtype, Nmax, xx, components = components)
            # Interpolate output to original x
            LP = n.interp(x, xx, LPi)
            LS = n.interp(x, xx, LSi)
            ECL = n.interp(x, xx, ECLi)
            REFL = n.interp(x, xx, REFLi)

            return LP,LS,ECL,REFL
        else:
            yy = call_task2(v, ldtype, Nmax, xx)
            # Interpolate output to original x
            return n.interp(x, xx, yy)

    else:
        return call_task2(v, ldtype, Nx, x, components = components)

    
def call_task2(v, ldtype, Npoints, time, components = False):
    import task2_components as task2
    
    # Define output arrays
    outflux = n.zeros(Npoints)
    LP = n.zeros(Npoints)
    LS = n.zeros(Npoints)
    REFL = n.zeros(Npoints)
    ECL = n.zeros(Npoints)
    
    ## Run task 2
    task2.task2(v, ldtype, Npoints, time, outflux, LP, LS, REFL, ECL)

    if components:
        return LP,LS,ECL,REFL
    else:
        return outflux        

def macula(t, star, spots, inst):
    """
    Compute spots model using macula
    https://www.cfa.harvard.edu/~dkipping/macula.html

    ! Istar 	= Theta_star(1)		! Inclination of the star [rads]
    ! Peq 		= Theta_star(2)		! Rot'n period of the star's equator [d]
    ! kappa2 	= Theta_star(3)		! Quadratic differential rotation coeff
    ! kappa4 	= Theta_star(4)		! Quartic differential rotation coeff
    ! c1 		= Theta_star(5)		! 1st of four-coeff stellar LD terms
    ! c2 		= Theta_star(6)		! 2nd of four-coeff stellar LD terms
    ! c3 		= Theta_star(7)		! 3rd of four-coeff stellar LD terms
    ! c4 		= Theta_star(8)		! 4th of four-coeff stellar LD terms
    ! d1 		= Theta_star(9)		! 1st of four-coeff spot LD terms
    ! d2	 	= Theta_star(10)	! 2nd of four-coeff spot LD terms
    ! d3	 	= Theta_star(11)	! 3rd of four-coeff spot LD terms
    ! d4	 	= Theta_star(12)	! 4th of four-coeff spot LD terms
    ! ------------------------------------------------------------------------------
    ! Theta_spot(j,k) = Parameters of the k^th spot
    ! ------------------------------------------------------------------------------
    ! Lambda0(k) 	= Theta_spot(1,k)	! Longitude of spot at time tref(k)
    ! Phi0(k) 	= Theta_spot(2,k)	! Latitude of spot at time tref(k)
    ! alphamax(k)	= Theta_spot(3,k)	! Angular spot size at time tmax(k)
    ! fspot(k)	= Theta_spot(4,k)	! Spot-to-star flux contrast of spot k
    ! tmax(k)	= Theta_spot(5,k)	! Time at which spot k is largest
    ! life(k)	= Theta_spot(6,k)	! Lifetime of spot k (FWFM) [days]
    ! ingress(k)	= Theta_spot(7,k)	! Ingress duration of spot k [days]
    ! egress(k)	= Theta_spot(8,k)	! Egress duration of spot k  [days]
    ! ------------------------------------------------------------------------------
    ! Theta_inst(j,m) = Instrumental/nuisance parameters
    ! ------------------------------------------------------------------------------
    ! U(m) 		= Theta_inst(1,m)	! Baseline flux level for m^th data set
    ! B(m) 		= Theta_inst(2,m)	! Blend factor for m^th data set
    """
    import pymacula
    
    # time range
    t_start = n.min(t)-0.05
    t_stop = n.max(t)+0.05

    derivatives = False
    temporal = False
    tdeltav = False
	
    output = [n.empty(6)]
		
    output = pymacula.maculamod.macula(t, derivatives,temporal,tdeltav,star,spots,inst,t_start,t_stop)
	
    return output[0]
