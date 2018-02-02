def arome(mean_anomaly, sma, inc, lda, ldc, beta0=2.0, Vsini=15.0, sigma0=8.0, zeta=0., Rp=0.1, Kmax=4., units='radian'):
	"""
	AROME: high-level API 
	Computes and return analytical Rossiter-McLaughlin signals
	adapted the CCF or the iodine cell technique.
	
	Inputs:
	-------
	mean_anomaly: float
	Orbital mean anomaly in radian or degree
	
	sma: float
	Semi-major axis in stellar radii
	
	inc: float
	Orbital inclination in radian or degree
	
	lda: float
	Spin-orbit angle (lambda) in radian or degree
	
	ldc: array of 2 or 4 elements
	Limb darkening coefficients (2) for a linear law, (4) for a non-linear law
	
	beta0: float
	width of the non-rotating star in km/s
	
	Vsini: float
	Stellar projected rotational velocity in km/s
	
	sigma0: float
	Width of the best Gaussian fit in km/s
	
	zeta: float
	Macro-turbulence parameter in km/s
	
	Rp: float
	Radius of the planet in solar radius
	
	Kmax: float
	Order of expansion for the Iodine cell technique
	
	units: string
	units of all the input angles (mean_anomaly, inc, lda)
	Possible values: 'degree' or 'radian' (default)
	
	
	Outputs:
	--------
	
	vccf: float
	Value of the RM effect measured by the CCF technique in km/s
	
	
	
	Credits:
	--------
	Author  G. Boue  EXOEarths, Centro de Astrofisica, Universidade do Porto.
	Python translation  A. Santerne EXOEarths, Centro de Astrofisica, Universidade do Porto.
	
	Copyright (C) 2012, CAUP
	email of the author : gwenael.boue@astro.up.pt
	email of the Python translator : alexandre.santerne@astro.up.pt
	
	
	This work has been supported by the European Research Council/European
	Community under the FP7 through Starting Grant agreement number 239953, as
	well as from Fundacao para a Ciencia e a Tecnologia (FCT) through program
	Ciencia 2007 funded by FCT/MCTES (Portugal) and POPH/FSE (EC), and in the
	form of grants reference PTDC/CTE-AST/098528/2008, SFRH/BPD/71230/2010, and
	SFRH/BPD/81084/2011.
	
	
	License of the file :
	---------------------
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
	"""
	import numpy as n

	if units == 'degree':
		# inclination
		inc *= n.pi/180.0           # radian */
		# spin-orbit angle
		lda *= n.pi/180.0           # radian */
		#Mean anomaly
		mean_anomaly *= n.pi/180.0
	
	#planet's coordinates */
	x0 = sma*(-n.cos(lda)*n.cos(mean_anomaly)+n.sin(lda)*n.sin(mean_anomaly)*n.cos(inc))
	y0 = sma*(-n.sin(lda)*n.cos(mean_anomaly)-n.cos(lda)*n.sin(mean_anomaly)*n.cos(inc))
	z0 = sma*n.sin(mean_anomaly)*n.sin(inc)
	#print x0, y0, z0

	EPSILON=1e-20
	rho  = n.sqrt(x0**2+y0**2)
	dmin = rho-Rp
	dmax = rho+Rp
	
	if all(z0 < 0): 
	    if n.isscalar(z0) : return 0.
	    else: return n.zeros(len(z0), float)
	elif all(dmin >= 1.0-EPSILON): 
	    if n.isscalar(dmin): return 0.
	    else: return n.zeros(len(dmin), float)
	else:
		
		rv_rm = n.zeros(len(dmin), float)
		
		parome={'beta0':beta0, 'Vsini':Vsini, 'sigma0':sigma0, 'zeta':zeta, 'Kmax':Kmax, 'Rp':Rp}
		
		# limb darkening
		limb_coef, limb_pow, kern_coef = arome_alloc_LDc(ldc)
		parome['limb_coef'] = limb_coef
		parome['limb_pow'] = limb_pow
		parome['kern_coef'] = kern_coef
		parome['nlimb'] = len(limb_coef)
		
		parome['Gaussfit_a0'] = setGaussfit_a0(parome)
		
		for i in n.where(dmin < 1.0-EPSILON)[0] :
		    parome = arome_calc_fvpbetap(x0[i],y0[i],z0[i], parome)
		    rv_rm[i] = arome_get_RM_CCF(parome)
		
		if n.isscalar(dmin): return rv_rm[0]
		else: return rv_rm


def arome_alloc_LDc(ldc):
	"""
	Define limb darkening coefficient for arome
	
	inputs:
	-------
	
	ldc: array of 2 or 4 floats
	Linear (2) or non-linear (4) limb darkening coefficients
	
	output:
	-------
	
	limb_coef: array
	limb darkening coefficients for arome
	
	limb_pow: array
	limb darkening power for arome
	"""
	import numpy as n
	if len(ldc) == 2:
		limb_coef = n.zeros(3, float)
		limb_pow = n.zeros(3, float)
		denom = n.pi*(1.0-ldc[0]/3.0-ldc[1]/6.0)
		limb_coef[0] = (1.0-ldc[0]-ldc[1])/denom
		limb_coef[1] = (ldc[0]+2.0*ldc[1])/denom
		limb_coef[2] = -1.*ldc[1]/denom
		limb_pow[0]  = 0
		limb_pow[1]  = 2
		limb_pow[2]  = 4
	elif len(ldc) == 4:
		limb_coef = n.zeros(5, float)
		limb_pow = n.zeros(5, float)
		denom = n.pi*(1.0-ldc[0]/5.0-ldc[1]/3.0-3.0*ldc[2]/7.0-ldc[3]/2.0)
		limb_coef[0] = (1.0-ldc[0]-ldc[1]-ldc[2]-ldc[3])/denom
		limb_coef[1] = ldc[0]/denom
		limb_coef[2] = ldc[1]/denom
		limb_coef[3] = ldc[2]/denom
		limb_coef[4] = ldc[3]/denom
		limb_pow[0]  = 0
		limb_pow[1]  = 1
		limb_pow[2]  = 2
		limb_pow[3]  = 3
		limb_pow[4]  = 4
	else: 
		raise IOError("arome_alloc_LDc argument must be a two or four parameters in an array")
	
	kern_coef = setrotkernel(limb_coef, limb_pow)
	#print len(kern_coef)
	return limb_coef, limb_pow, kern_coef
	

def setrotkernel(limb_coef, limb_pow):
	import numpy as n
	Im2 = n.pi # int[-1,1] 1/sqrt(1-x**2) dx
	Im1 = 2.3962804694711844 # int[-1,1] 1/(1-x**2)**(1/4) dx
	Im0 = 2.0 # int[-1,1] dx
	Ip1 = 1.7480383695280799 # int[-1,1] (1-x**2)**(1/4) dx
	
	ntabI = 4
	for k in xrange(len(limb_pow)):
		ntabI = max(ntabI, limb_pow[k]+4)
	
	tabI = n.zeros(ntabI, float)
	tabI[-2] = Im2
	tabI[-1] = Im1
	tabI[0] = Im0
	tabI[1] = Ip1
	for k in xrange(2, int(ntabI)-2): 
		tabI[k] = float(k)/float(k+2)*tabI[k-4] # int[-1,1] (1-x**2)**(k/4) dx
	#print tabI

	kern_coef = n.zeros(len(limb_pow))
	for k in xrange(len(limb_pow)):
		alpha = limb_pow[k]
		kern_coef[k] = limb_coef[k]*tabI[alpha]
	#print kern_coef
	return kern_coef

def setGaussfit_a0(parome):
	import numpy as n
	from scipy import integrate
	res = 4.0*parome['sigma0']*n.sqrt(n.pi)*integrate.quad(funcAmp_a0, 0.0, 1.0, parome, 20)[0]
	#res = 4.0*parome['sigma0']*n.sqrt(n.pi)*qtrap(funcAmp_a0, 0.0, 1.0, parome)
	return res
	
def funcAmp_a0(x, parome):
	import numpy as n
	sig2 = parome['sigma0']**2+parome['beta0']**2+parome['zeta']**2/2.0
	mu   = n.sqrt(1.0-x**2)
	smu  = n.sqrt(mu)
	
	#Rotation kernel
	Rx = 0.0
	for k in xrange(parome['nlimb']):
		Rx += parome['kern_coef'][k]*smu**parome['limb_pow'][k]
	Rx *= mu
	# Gaussian
	
	Gx = 1.0/n.sqrt(2.0*n.pi*sig2)*n.exp(-1.*(x*parome['Vsini'])**2/(2.0*sig2))
	return Rx*Gx


def setIodine_den(parome):
	import numpy as n
	from scipy import integrate
	return 2.0*integrate.fixed_quad(funcIodine_den, 0.0, 1.0, parome, 20)[0]


def funcIodine_den(x, parome):
	import numpy as n
	sig2 = parome['beta0']**2+parome['zeta']**2/2.0
	v = 1.0/x-1.0
	res = dGconvR(n.sqrt(sig2), v*parome['Vsini'], parome)
	res /= x
	res = res**2*parome['Vsini']
	return res


def dGconvR(sig, v, parome):
	"""
	Computes [dG_sig/dv * R](v)
	where R is the rotation kernel
	returns 0 if no problem
	
	inputs:
	-------
	sig   : width of the Gaussian
	v     : coordinates at which the cross-correlation is evaluated
	pdata : astronomical constants containing the limb-darkening coefficients
	
	output:
	-------
	res   : result
	"""

	import numpy as n
	s = sig/parome['Vsini']
	x = v/parome['Vsini']
	x1= -1.0*(1.0+x)/s
	x2=  (1.0-x)/s
	mapp = n.ones(parome['nlimb'], float) # (-alpha)_p/p! = (-alpha)(-alpha+1)...(-alpha+p-1)
	res = 0.0
	fac = 1.0
	for m in xrange(2*parome['Kmax']):
		p = m/2.
		pow2 = 2.0**p
		F1 = gamma_inc(m/2.0+1.0,x1*x1/2.0)
		F2 = gamma_inc(m/2.0+1.0,x2*x2/2.0)
		F1 *= pow2
		F2 *= pow2
		hyp=0.0
		if m%2 == 0:
			for k in xrange(parome['nlimb']):
				if mapp[k]:
					alpha = (parome['limb_pow'][k]+1.0)/2.0
					hyp += parome['kern_coef'][k]*mapp[k]*hypertrunc(-alpha+p, p+0.5, 0.5, x**2, parome['Kmax']-p)
			if p : fac/=p
		elif m%2 == 1:
			if x1 > 0.0 : F1 *= 1.0 
			else : F1 *= -1.0
			if x2 > 0.0 : F2 *= 1.0 
			else : F2 *= -1.0
			for k in xrange(parome['nlimb']):
				if mapp[k]:
					alpha = (parome['limb_pow'][k]+1.0)/2.0
					mapp[k] *= -alpha+p
					hyp += parome['kern_coef'][k]*mapp[k]*2.0*x*hypertrunc(-alpha+p+1.0, p+1.5, 1.5, x**2, parome['Kmax']-p-1.)
		if m : fac*=s
		
		res += hyp*fac*(F2-F1)
	return res / (n.sqrt(2.0*n.pi)*s*parome['Vsini']**2)






def arome_calc_fvpbetap(x,y,z,parome):
	"""
	Computes the flux f, the subplanet velocity vp, and dispersions betapR, betapT
	
	x coordinate of the planet in stellar radius
	y coordinate of the planet in stellar radius
	z coordinate of the planet in stellar radius
	
	r        :     largest radius of the planet                                 
	rx       :     x radius of the rotated projected planet                     
	ry       :     y radius of the rotated projected planet                     
	phi0     :     atan2(x0, y0)                                                
	rho      :     sqrt(x0**2+y0**2)                                            
	psi      :     angle OPA (O=star, P=planet, A=intersection star/planet)     
	phi1     :     limit inf of an arc                                          
	phi2     :     limit sup of an arc                                          
	xbar     :     averaged x coordinate of the planet                          
	ybar     :     averaged y coordinate of the planet                          
	II       :     II = I(xbar, ybar) (mean subplanet intensity)                
	Hxx      :     partial_x partial_x I(x,y)                                   
	Hyy      :     partial_y partial_y I(x,y)                                   
	Hxy      :     partial_x partial_y I(x,y)                                   
	Hxx2     :     partial_x partial_x I(x,y)                                   
	Hyy2     :     partial_y partial_y I(x,y)                                   
	Hxy2     :     partial_x partial_y I(x,y)                                   
	a00      :     area covered by the planet                                   
	axx      :     averaged of (x-xbar)**2                                      
	ayy      :     averaged of (y-ybar)**2                                      
	axy      :     averaged of (x-xbar)(y-ybar)                                 
	ff       :     fraction of flux oculted by the planet                       
	vv       :     subplanet velocity                                           
	v2       :     square of the subplanet velocity                             
	dmin     :     minimal distance between the planet ellipse and the origin   
	dmax     :     maximal distance between the planet ellipse and the origin   
	dbetaR   :     dispersion in the radial direction                           
	dbetaT   :     dispersion in the radial direction                           
	zetaR2   :     radial dispersion due to macro-turbulence                    
	zetaT2   :     tangential dispersion due to macro-turbulence                
	"""
	import numpy as n
	#import math
	EPSILON=1e-20
	
	if z <= 0.0: # the planet is behind the star
		parome['flux']= 0.0
		parome['vp']= 0.0
		parome['betapR'] = parome['beta0']
		parome['betapT'] = parome['beta0']
		return parome


	# parameters of the planet
	rx, ry, r = parome['Rp'],parome['Rp'],parome['Rp']
	phi0 = n.arctan2(y,x)
	rho  = n.sqrt(x**2+y**2)
	dmin = rho-r
	dmax = rho+r

	if dmin >= 1.0-EPSILON: # the planet doesn't overlap the stellar disk, that's life !
		parome['flux']= 0.0
		parome['vp']= 0.0
		parome['betapR'] = parome['beta0']
		parome['betapT'] = parome['beta0']
		return parome
	
	elif dmax<=1.0: #the planet is completely inside the stellar disk, that's cool !
		xbar=x # int x dxdy / int dxdy
		ybar=y # int y dxdy / int dxdy
		a00 = n.pi*rx*ry # int dxdy
		axx = rx**2/4.0  # int (x-x0)**2 dxdy / int dxdy
		ayy = ry**2/4.0  # int (y-y0)**2 dxdy / int dxdy
		axy = 0.0        # int (x-x0)*(y-y0) dxdy / int dxdy
	
	else : #during ingress and egress 
		
		#stellar boundary
		psi   = n.arccos((1.0+rho**2-r**2)/(2.0*rho)) # angle BSP (see Fig. 1)
		phi1  = phi0-psi # angle xSB (see Fig. 1)
		phi2  = phi0+psi # angle xSA (see Fig. 1)
		
		#print phi0, phi1, phi2
		
		a00  = funcF00(0.0,0.0,1.0,1.0,phi2)-funcF00(0.0,0.0,1.0,1.0,phi1)
		xbar = funcF10(0.0,0.0,1.0,1.0,phi2)-funcF10(0.0,0.0,1.0,1.0,phi1)
		ybar = funcF01(0.0,0.0,1.0,1.0,phi2)-funcF01(0.0,0.0,1.0,1.0,phi1)
		axx  = funcF20(0.0,0.0,1.0,1.0,phi2)-funcF20(0.0,0.0,1.0,1.0,phi1)
		ayy  = funcF02(0.0,0.0,1.0,1.0,phi2)-funcF02(0.0,0.0,1.0,1.0,phi1)
		axy  = funcF11(0.0,0.0,1.0,1.0,phi2)-funcF11(0.0,0.0,1.0,1.0,phi1)
		
		# planet boundary
		psi   = n.arccos(-1.*(1.0-rho**2-r**2)/(2.0*r*rho)) # angle APS (see Fig. 1) in [0,pi]
		phi1  = phi0+n.pi-psi # angle xPA (see Fig. 1)
		phi2  = phi0+n.pi+psi # angle xPB (see Fig. 1)
		
		a00  += (funcF00(x,y,rx,ry,phi2)-funcF00(x,y,rx,ry,phi1))
		xbar += (funcF10(x,y,rx,ry,phi2)-funcF10(x,y,rx,ry,phi1))
		ybar += (funcF01(x,y,rx,ry,phi2)-funcF01(x,y,rx,ry,phi1))
		axx  += (funcF20(x,y,rx,ry,phi2)-funcF20(x,y,rx,ry,phi1))
		ayy  += (funcF02(x,y,rx,ry,phi2)-funcF02(x,y,rx,ry,phi1))
		axy  += (funcF11(x,y,rx,ry,phi2)-funcF11(x,y,rx,ry,phi1))
		
		#print a00, xbar, ybar, axx, ayy, axy
		
		xbar /= a00
		ybar /= a00
		axx   = axx/a00 - xbar**2
		ayy   = ayy/a00 - ybar**2
		axy   = axy/a00 - xbar*ybar
		
		#print xbar, ybar, axx, ayy, axy
		
	II = funcIxn(xbar,ybar,parome,0)
	Hxx0,Hyy0,Hxy0 = HessIxn(xbar,ybar,parome,0)
	ff = a00*(II+0.5*(Hxx0*axx+Hyy0*ayy+2.0*Hxy0*axy))
	
	
	Hxx1,Hyy1,Hxy1 = HessIxn(xbar,ybar,parome,1)
	Hxx1 -= xbar*Hxx0
	Hyy1 -= xbar*Hyy0
	Hxy1 -= xbar*Hxy0
	vv    = xbar + 0.5/II*(Hxx1*axx+Hyy1*ayy+2.0*Hxy1*axy)
	
	#print II, Hxx1*axx, Hyy1*ayy, Hxy1*axy, vv
	
	Hxx2,Hyy2,Hxy2 = HessIxn(xbar,ybar,parome,2)

	Hxx2 -= xbar**2*Hxx0
	Hyy2 -= xbar**2*Hyy0
	Hxy2 -= xbar**2*Hxy0
	v2    = xbar**2 + 0.5/II*(Hxx2*axx+Hyy2*ayy+2.0*Hxy2*axy)
	
	#print II, axx, ayy, axy, v2
	#print v2 - vv**2
	
	# results
	
	parome['flux'] = ff
	parome['vp']   = vv
	dbetaR         = n.sqrt(v2-vv**2)
	dbetaT         = dbetaR
	
	#print ff, vv, dbetaR, dbetaT
	
	# set the units
	
	parome['vp'] *= parome['Vsini']
	dbetaR      *= parome['Vsini']
	dbetaT      *= parome['Vsini']

	if parome['zeta']>0.0 : #take into account macro turbulence
		limb_pow = n.zeros(parome['nlimb'], float)
		mu2bar = 1.0-xbar**2-ybar**2
		
		#multiply I(x,y) by cos**2(theta)
		for k in xrange(parome['nlimb']):
			limb_pow[k] = parome['limb_pow'][k]
			parome['limb_pow'][k] += 4
		
		Hxx2,Hyy2,Hxy2 = HessIxn(xbar,ybar,parome,0)
		Hxx2 -= mu2bar*Hxx0
		Hyy2 -= mu2bar*Hyy0
		Hxy2 -= mu2bar*Hxy0
		zetaR2 = mu2bar + 0.5/II*(Hxx2*axx+Hyy2*ayy+2.0*Hxy2*axy)
		zetaT2 = 1.0-zetaR2
		
		zetaR2 *= parome['zeta']**2
		zetaT2 *= parome['zeta']**2
		
		dbetaR = n.sqrt(dbetaR**2+zetaR2)
		dbetaT = n.sqrt(dbetaT**2+zetaT2)
		
		# retrieve the initial limb-darkening law
		for k in xrange(parome['nlimb']):
			parome['limb_pow'][k] = limb_pow[k]
		
	
	# add to the width of the non-rotating star
	parome['betapR'] = n.sqrt(dbetaR**2+parome['beta0']**2)
	parome['betapT'] = n.sqrt(dbetaT**2+parome['beta0']**2)
	
	return parome


def funcF00(x0, y0, rx, ry, phi):
	"""
	Computes F_00(phi) (Boue et al., 2012, Tab1) use for the covered surface
	x_0 : x coordinates of the planet (or of the star)
	y_0 : y coordinates of the planet (or of the star)
	rx  : x radius of the planet (or of the star)
	ry  : y radius of the planet (or of the star)
	phi : limit of the arc 
	"""
	import numpy as n
	return 0.5*(rx*ry*phi+x0*ry*n.sin(phi)-y0*rx*n.cos(phi))


def funcF10(x0,  y0, rx, ry, phi):
	"""
	Computes F_10(phi) (Boue et al., 2012, Tab1) use for the averaged velocity
	x_0 : x coordinates of the planet (or of the star)
	y_0 : y coordinates of the planet (or of the star)
	rx  : x radius of the planet (or of the star)
	ry  : y radius of the planet (or of the star)
	phi : limit of the arc 
	"""
	import numpy as n
	return -x0*y0*rx*n.cos(phi)\
		-0.5*y0*rx**2*(n.cos(phi))**2\
		+0.25*x0*rx*ry*(2.0*phi-n.sin(2.0*phi))\
		+1.0/12.0*rx**2*ry*(3.0*n.sin(phi)-n.sin(3.0*phi))

def funcF01(x0, y0, rx, ry, phi):
	"""
	Computes F_01(phi) (Boue et al., 2012, Tab1) use for the averaged velocity
	x_0 : x coordinates of the planet (or of the star)
	y_0 : y coordinates of the planet (or of the star)
	rx  : x radius of the planet (or of the star)
	ry  : y radius of the planet (or of the star)
	phi : limit of the arc 
	"""
	import numpy as n
	return  x0*y0*ry*n.sin(phi)\
		+0.5*x0*ry**2*(n.sin(phi))**2\
		+0.25*y0*rx*ry*(2.0*phi+n.sin(2.0*phi))\
		-1.0/12.0*rx*ry**2*(3.0*n.cos(phi)+n.cos(3.0*phi))


def funcF20(x0, y0, rx, ry, phi):
	"""
	Computes F_20(phi) (Boue et al., 2012, Tab1) use for the velocity dispersion
	x_0 : x coordinates of the planet (or of the star)
	y_0 : y coordinates of the planet (or of the star)
	rx  : x radius of the planet (or of the star)
	ry  : y radius of the planet (or of the star)
	phi : limit of the arc 
	"""
	import numpy as n
	return -1.*x0**2*y0*rx*n.cos(phi)\
		-x0*y0*rx**2*(n.cos(phi))**2\
		+0.25*x0**2*rx*ry*(2.0*phi-n.sin(2.0*phi))\
		-1.0/12.0*y0*rx**3*(3.0*n.cos(phi)+n.cos(3.0*phi))\
		+1.0/6.0*x0*rx**2*ry*(3.0*n.sin(phi)-n.sin(3.0*phi))\
		+1.0/32.0*rx**3*ry*(4.0*phi-n.sin(4.0*phi))


def funcF02(x0, y0, rx, ry, phi):
	"""
	Computes F_02(phi) (Boue et al., 2012, Tab1) use for the velocity dispersion
	x_0 : x coordinates of the planet (or of the star)
	y_0 : y coordinates of the planet (or of the star)
	rx  : x radius of the planet (or of the star)
	ry  : y radius of the planet (or of the star)
	phi : limit of the arc 
	"""
	import numpy as n
	return  x0*y0**2*ry*n.sin(phi)\
		+x0*y0*ry**2*(n.sin(phi))**2\
		+0.25*y0**2*rx*ry*(2.0*phi+n.sin(2.0*phi))\
		+1.0/12.0*x0*ry**3*(3.0*n.sin(phi)-n.sin(3.0*phi))\
		-1.0/6.0*y0*rx*ry**2*(3.0*n.cos(phi)+n.cos(3.0*phi))\
		+1.0/32.0*ry**3*rx*(4.0*phi-n.sin(4.0*phi))




def funcF11(x0, y0, rx, ry, phi):
	"""
	Computes F_11(phi) (Boue et al., 2012, Tab1) use for the velocity dispersion
	x_0 : x coordinates of the planet (or of the star)
	y_0 : y coordinates of the planet (or of the star)
	rx  : x radius of the planet (or of the star)
	ry  : y radius of the planet (or of the star)
	phi : limit of the arc 
	"""
	import numpy as n
	return 0.25*x0*y0*(2.0*rx*ry*phi+x0*ry*n.sin(phi)-y0*rx*n.cos(phi))\
		+0.125*(x0*ry*n.sin(phi))**2\
		-0.125*(y0**2+ry**2)*(rx*n.cos(phi))**2\
		+1.0/48.0*y0*rx**2*ry*(15.0*n.sin(phi)-n.sin(3*phi))\
		-1.0/48.0*x0*rx*ry**2*(15.0*n.cos(phi)+n.cos(3*phi))



def funcIxn(x0, y0, parome, n):
	"""
	Computes x**n*B(x,y)*I(x,y) at the position (x0,y0)
	x0    : x coordinates of the planet
	y0    : y coordinates of the planet
	pdata : limb-darkening coefficient
	n     : power of x
	"""
	return x0**n*funcLimb(x0,y0,parome)



def HessIxn(x0, y0, parome, n):
	"""
	Computes the Hessian of x**n*I(x,y) at the position (x0,y0)
	
	inputs:
	-------
	x0  :   x coordinates of the planet
	y0  :   y coordinates of the planet
	c1  :   limb-darkening coefficient
	c2  :   limb-darkening coefficient
	n   :   power of x
	
	outputs:
	--------
	Hxx :  partial_x partial_x (x^nI(x,y))
	Hyy :  partial_y partial_y (x^nI(x,y))
	Hxy :  partial_x partial_y (x^nI(x,y))
	"""
	xn  = x0**n
	
	L = funcLimb(x0,y0,parome)
	Lx, Ly = dfuncLimb(x0,y0,parome)
	Lxx, Lyy, Lxy= ddfuncLimb(x0,y0,parome)
	
	#if n in [1,2] : print Lxx, Lyy, Lxy
	
	Hxx = xn*Lxx
	if n>0 : Hxx += 2.0*float(n)*x0**(n-1)*Lx
	if n>1 : Hxx += L*float(n)*float(n-1)*x0**(n-2)
	Hyy = xn*Lyy
	Hxy = xn*Lxy
	if n>0 : Hxy += Ly*float(n)*x0**(n-1)
	
	return Hxx, Hyy, Hxy


def hypertrunc(a, b, c, x, Kmax):
	"""
	Return the hypergeometric series F(a,b,c;x) truncated sum_{k=0}^parome['Kmax']
	"""
	import numpy as n
	res, fac=1.0, 1.0
	aa=a
	bb=b
	cc=c
	for k in xrange(1,int(Kmax)):
		fac*=aa*bb/(cc*k)*x
		res+=fac
		aa+=1.0
		bb+=1.0
		cc+=1.0
	
	return res



def funcLimb(x0, y0, parome):
	"""
	Computes I(x,y)=sum cn*mu**(n/2) at the position (x0,y0)
	x0    : x coordinates of the planet
	y0    : y coordinates of the planet
	pdata : limb-darkening coefficients
	"""
	import numpy as n
	R   = x0**2+y0**2
	mus = (1.0-R)**(1./4.)
	res = 0.0
	for k in xrange(parome['nlimb']):
		res += parome['limb_coef'][k]*mus**parome['limb_pow'][k]
	return res




def dfuncLimb(x0, y0, parome):
	"""
	Computes the first derivatives of I(x,y)=sum cn*mu**(n/2) at the position (x0,y0)
	
	inputs:
	-------
	x0    :  x coordinates of the planet
	y0    :  y coordinates of the planet
	pdata :  limb-darkening coefficients
	
	outputs:
	--------
	Jx    : partial_x I(x,y)
	Jy    : partial_y I(x,y)
	"""
	
	R = x0**2+y0**2
	mu2 = 1.0-R
	mus = mu2**(1./4.)
	dIdR = 0.0
	#print R, mu2, mus, dIdR
	for k in xrange(parome['nlimb']):
		dIdR -= 0.25*parome['limb_pow'][k]*parome['limb_coef'][k]*mus**parome['limb_pow'][k]/mu2
	Jx = 2.0*x0*dIdR
	Jy = 2.0*y0*dIdR
	
	return Jx, Jy



def ddfuncLimb(x0, y0, parome):
	"""
	Computes the second derivatives of I(x,y)=sum cn*mu**(n/2) at the position (x0,y0)
	
	inputs:
	-------
	x0    :   x coordinates of the planet
	y0    :   y coordinates of the planet
	pdata :   limb-darkening coefficient
	
	outputs:
	--------
	Hxx   :  partial_x partial_x (I(x,y))
	Hyy   :  partial_y partial_y (I(x,y))
	Hxy   :  partial_x partial_y (I(x,y))
	"""
	
	R   = x0**2+y0**2
	mu2 = 1.0-R
	mus = mu2**(1./4.)
	IR = 0.0
	IRR = 0.0
	
	#print R, mu2, mus, IR, IRR
	
	for k in xrange(parome['nlimb']):
		var  = 0.25*parome['limb_pow'][k]*parome['limb_coef'][k]*mus**parome['limb_pow'][k]/mu2
		IR  -= var
		IRR += var*(0.25*parome['limb_pow'][k]-1.0)/mu2
	
	
	Hxx = 2.0*IR+4.0*x0**2*IRR
	Hyy = 2.0*IR+4.0*y0**2*IRR
	Hxy = 4.0*x0*y0*IRR
	
	return Hxx, Hyy, Hxy
	



def arome_get_RM_CCF(parome):
	"""
	Computes the RM effect measured by the CCF technique.
	v = 1/den * (2*sig0**2/(sig0**2+betap**2))**(3/2)*f*vp*exp(-vp**2/(2*(sig0**2+betap**2)))
	
	input:
	------
	parome :  simulation structure
	
	output:
	-------
	res    : result
	"""
	import numpy as n
	den   = parome['Gaussfit_a0']
	f     = parome['flux']
	vp    = parome['vp']
	bs2   = parome['sigma0']**2
	bpR2  = parome['betapR']**2
	bpT2  = parome['betapT']**2
	if f<0.0 :
		return 0
	if parome['zeta'] > 0.0: # with macro-turbulence
		return -0.5/den*(2.0*bs2/(bs2+bpR2))**(3./2.)*f*vp*n.exp(-vp**2/(2.0*(bs2+bpR2))) -0.5/den*(2.0*bs2/(bs2+bpT2))**(3./2.)*f*vp*n.exp(-vp**2/(2.0*(bs2+bpT2)))
	else : # without macro-turbulence
		return -1.0/den*(2.0*bs2/(bs2+bpR2))**(3./2.)*f*vp*n.exp(-vp**2/(2.0*(bs2+bpR2)))



def qtrap(func, a, b, pdata):
	"""
	Returns the integral of the function func from a to b. The parameters EPS can be
	set to the desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the
	maximum allowed number of steps. Integration is performed by the trapezoidal rule.
	
	(C) Copr. 1986-92 Numerical Recipes Software
	"""
	EPS = 1e-5
	JMAX = 20
	s,olds=0.0, 0.0

	for j in xrange(1,JMAX):
		if j == 1 : s = None
		s=trapzd(func,a,b,pdata,j,s)
		#print s
		if j > 5:
			if abs(s-olds) < EPS*abs(olds) or (s == 0.0 and olds == 0.0) : return s
			else : olds=s
		else : olds=s
	return s


def trapzd(func, a, b, pdata, n, s):
	"""
	This routine computes the nth stage of refinement of an extended trapezoidal rule.
	func is input as a pointer to the function to be integrated between limits a and b, 
	also input. When called with n=1, the routine returns the crudest estimate of int[a,b] 
	f(x)dx. Subsequent calls with n=2,3,... (in that sequential order) will improve the 
	accuracy by adding 2**(n-2) additional interior points.
	
	(C) Copr. 1986-92 Numerical Recipes Software
	"""
	
	if n==1: 
		return 0.5*(b-a)*(func(a,pdata)+func(b,pdata))
	elif n>1:
		t = 1
		for j in xrange(1,n-1): 
			t = t << 1
		tnm=t
		#print tnm
		del_v=(b-a)/tnm
		#print del_v
		x=a+0.5*del_v
		#print x
		sum_v = 0.
		for j in xrange(1, t):
			x+=del_v
			sum_v += func(x,pdata)
			#print sum_v
		return 0.5*(s+(b-a)*sum_v/tnm)



def gammln(xx):
	"""
	Returns the value ln((Gamma(xx))) for xx > 0.
	
	(C) Copr. 1986-92 Numerical Recipes Software
	"""
	import numpy as n
	cof=[76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5]
	y=xx, x=xx
	tmp=x+5.5
	tmp -= (x+0.5)*n.log(tmp)
	ser=1.000000000190015
	for j in xrange(len(cof)) : ser += cof[j]/(y+j+1.)
	return -tmp+log(2.5066282746310005*ser/x)

def gamma_inc(a, x, parome):
	"""
	Computes the incomplete gamma function P(a,x)
	
	(C) Copr. 1986-92 Numerical Recipes Software
	"""
	import numpy as n
	if x < 0.0 or a <= 0.0:
		raise IOError("Invalid arguments (%lg,%lg) in routine gammp." %(a, x))
	if x < (a+1.0):
		gamser, gln = gser(a,x,parome)
		return gamser*n.exp(gln)
	else :
		gammcf , gln = gcf(a,x,parome)
		return (1.0-gammcf)*n.exp(gln)


def gser(a, x, pdata):
	"""
	Returns the incomplete gamma function P(a,x) evaluated by its series 
	representation as gamser. Also returns ln(Gamma(a)) as gln.
	
	this function has been modified to handle the errors.
	
	(C) Copr. 1986-92 Numerical Recipes Software
	"""
	import numpy as n
	ITMAX = 100.
	EPS = 3e-7
	gln=gammln(a)
	if x == 0.0:
		gamser=0.0
	if x < 0.0:
		raise IOError("x (%lg) less than 0 in routine gser"%x)
	else :
		ap=a
		del_v=1.0/a
		sum_v=1.0/a
		for n in xrange(1,ITMAX+1):
			ap+=1
			del_v *= x/ap
			sum_v += del_v
			if abs(del_v) < abs(sum_v)*EPS :
				gamser=sum_v*n.exp(-x+a*n.log(x)-gln)
	return gamser, gln



def gcf(a, x, pdata):
	"""
	Returns the incomplete gamma function Q(a,x) evaluated by its continued
	fraction representation as gammcf. Also returns ln Gamma(a) as gln.

	this function has been modified to handle the errors.
	
	(C) Copr. 1986-92 Numerical Recipes Software
	"""
	ITMAX = 100.
	EPS = 3e-7
	FPMIN = 1e-30
	gln=gammln(a)
	b=x+1.0-a
	c=1.0/FPMIN
	d=1.0/b
	h=d
	for i in xrange(1, ITMAX+1):
		an = -i*(i-a)
		b += 2.0
		d=an*d+b
		if abs(d) < FPMIN : d=FPMIN
		c=b+an/c
		if abs(c) < FPMIN : c=FPMIN
		d=1.0/d
		del_v=d*c
		h *= del_v
		if abs(del_v-1.0) < EPS: break
	if i > ITMAX : raise ValueError("a (%lg) too large, ITMAX (%d) too small in gcf." %(a, ITMAX))
	gammcf=n.exp(-x+a*n.log(x)-(gln))*h
	return gammcf, gln
