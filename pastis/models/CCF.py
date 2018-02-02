"""
2012-04-16 : program birth
"""


def get_FWHM(vsini, BmV, spectrograph):
	"""
	Return the FWHM of a CCF assuming a given vsini.
	
	Calibrated for HARPS, SOPHIE HE and SOPHIE HR
	"""
	import numpy as n
	if spectrograph == 'HARPS' : 
		sigma0_harps=8.625-20.037*BmV+23.388*BmV**2-10.364*BmV**3+1.273*BmV**4
		A=1.95
		FWHM=n.sqrt((vsini**2.+A**2.*sigma0_harps**2.)/A**2.)*(2.*n.sqrt(2.*n.log(2.)))
	if spectrograph == 'SOPHIE HE' : 
		sigma0_HE=10.52-22.56*BmV+22.37*BmV**2.-6.95*BmV**3.
		A=1.64
		FWHM=n.sqrt((vsini**2.+A**2.*sigma0_HE**2.)/A**2.)*(2.*n.sqrt(2.*n.log(2.)))
	if spectrograph == 'SOPHIE HR' : 
		sigma0_HR=9.90 - 22.56*BmV + 22.37*BmV**2. - 6.95*BmV**3.
		A=1.95
		FWHM=n.sqrt((vsini**2.+A**2.*sigma0_HR**2.)/A**2.)*(2.*n.sqrt(2.*n.log(2.)))
	return FWHM

def get_contrast(FWHM, BmV, z, spectrograph, mask):
	"""
	Return the contrast of a CCF assuming a given FWHM, (B-V) and metallicity z.
	
	Calibrated for SOPHIE and HARPS
	"""
	import numpy as n
	if spectrograph == 'HARPS':

		if mask == 'K5' : c=[-0.4232454896,1.688452363,0.2941150665,-0.7595982552,-0.05445498228,-0.1592344046] 
		if mask == 'G2' : c=[-0.1295745075,1.534696698,0.2792906761,-0.7796853781,-0.05642914772,-0.1540470123]
		
		logW=c[0]+c[1]*BmV+c[2]*z+c[3]*BmV**2+c[4]*z**2+c[5]*BmB*z
		contrast = 10**logW / n.sqrt(2.*n.pi) / (FWHM / (2.*n.log(2.))) * 100. # in % 

	if spectrograph.find('SOPHIE')>-1:

		if mask == 'K5':
			logW = (z + 2.2606*BmV - 0.2615) / 3.9553  #Boisse et al. 2010
			contrast = 10**logW / n.sqrt(2.*n.pi) / (FWHM / (2.*n.log(2.))) * 100. # in % 

		if mask == 'G2':
			logW = (z + 1.4992*BmV + 0.9440) / 3.8807 #Boisse et al. 2010
			contrast = 10**logW / n.sqrt(2.*n.pi) / (FWHM / (2.*n.log(2.))) * 100. # in % 

	return contrast
	
def make_CCF(rv, rv0, FWHM, contrast, flux):
	"""
	Make a synthetic gaussian CCF for a given radial velocity, FWHM, contrast and flux.
	"""
	import numpy as n
	sigma = FWHM / 2.3548
	return (1.-contrast/100.*n.exp(-(rv-rv0)**2/(2.*sigma**2)))*flux


def PASTIS_CCF(t_rv, RV_CCF, spectrograph, *args):
	"""
	return models for RV, bisector, FWHM and contrast
	"""
	import numpy as n
	import pdb
	
	mask = spectrograph[0].split('-')[1]
	spectro = spectrograph[0].split('-')[2]
	offset = spectrograph[1]
	normalized_CCF = n.zeros([r_rv, RV_CCF], float)
	if 1 :
		
		# count number of stars
		nb_star = 0
		for obj in args :
			if obj.type == 'star' : nb_star+=1
			if obj.type == 'binary' : nb_star+=2
			if obj.type == 'triple' : 
				if obj.object1.type == 'star' : nb_star+=1
				if obj.object1.type == 'binary' : nb_star+=2
				if obj.object1.type == 'plansys' : nb_star+=1
				if obj.object2.type == 'star' : nb_star+=1
				if obj.object2.type == 'binary' : nb_star+=2
				if obj.object2.type == 'plansys' : nb_star+=1
			if obj.type == 'plansys' : nb_star+=1
		FWHM = n.zeros(nb_star, float)
		contrast = n.zeros(nb_star, float)
		CCF = n.zeros([nb_star, len(RV_CCF)], float)
		rv0 = n.zeros([nb_star, len(t_rv)], float)
		flux = n.zeros(nb_star, float)
		# define CCF proprieties of each stars
		count = 0
		for obj in args :
			if obj.type == 'star' : 
				FWHM[count] = get_FWHM(obj.vsini, obj.BmV, spectro)
				contrast[count] = get_contrast(FWHM[count], obj.BmV, obj.z, spectro, mask)
				flux[count] = obj.get_flux('Johnson-V')
				rv0[count] = obj.v0*ones(len(t_rv), float)
				count+=1
			if obj.type == 'binary' : 
				FWHM[count] = get_FWHM(obj.star1.vsini, obj.star1.BmV, spectro)
				contrast[count] = get_contrast(FWHM[count], obj.star1.BmV, obj.star1.z, spectro, mask)
				flux[count] = obj.star1.get_flux('Johnson-V')
				rv0[count] = obj.get_RV(t_rv, component = 'primary')
				count+=1
				FWHM[count] = get_FWHM(obj.star2.vsini, obj.star2.BmV, spectro)
				contrast[count] = get_contrast(FWHM[count], obj.star2.BmV, obj.star2.z, spectro, mask)
				flux[count] = obj.star2.get_flux('Johnson-V')
				rv0[count] = obj.get_RV(t_rv, component = 'secondary')
				count+=1
			if obj.type == 'triple' : 
				if obj.object1.type == 'star' : 
					FWHM[count] = get_FWHM(obj.object1.vsini, obj.object1.BmV, spectro)
					contrast[count] = get_contrast(FWHM[count], obj.object1.BmV, obj.object1.z, spectro, mask)
					flux[count] = obj.object1.get_flux('Johnson-V')
					rv0[count] = obj.get_RV(t_rv, component = 'primary')
					count+=1
				if obj.object1.type == 'binary' : 
					FWHM[count] = get_FWHM(obj.object1.star1.vsini, obj.object1.star1.BmV, spectro)
					contrast[count] = get_contrast(FWHM[count], obj.object1.star1.BmV, obj.object1.star1.z, spectro, mask)
					flux[count] = obj.object1.star1.get_flux('Johnson-V')
					rv0[count] = obj.get_RV(t_rv, component = 'primary')
					count+=1
					FWHM[count] = get_FWHM(obj.object1.star2.vsini, obj.object1.star2.BmV, spectro)
					contrast[count] = get_contrast(FWHM[count], obj.object1.star2.BmV, obj.object1.star2.z, spectro, mask)
					flux[count] = obj.object1.star2.get_flux('Johnson-V')
					rv0[count] = obj.get_RV(t_rv, component = 'secondary')
					count+=1
				if obj.object1.type == 'plansys' : 
					FWHM[count] = get_FWHM(obj.object1.star.vsini, obj.object1.star.BmV, spectro)
					contrast[count] = get_contrast(FWHM[count], obj.object1.star.BmV, obj.object1.star.z, spectro, mask)
					flux[count] = obj.object1.star.get_flux('Johnson-V')
					rv0[count] = obj.get_RV(t_rv, component = 'primary')
					count+=1
				if obj.object2.type == 'star' : 
					FWHM[count] = get_FWHM(obj.object2.vsini, obj.object2.BmV, spectro)
					contrast[count] = get_contrast(FWHM[count], obj.object2.BmV, obj.object2.z, spectro)
					flux[count] = obj.object2.get_flux('Johnson-V')
					rv0[count] = obj.get_RV(t_rv, component = 'secondary')
					count+=1
				if obj.object2.type == 'binary' : 
					FWHM[count] = get_FWHM(obj.object2.star1.vsini, obj.object2.star1.BmV, spectro)
					contrast[count] = get_contrast(FWHM[count], obj.object2.star1.BmV, obj.object2.star1.z, spectro, mask)
					flux[count] = obj.object2.star1.get_flux('Johnson-V')
					rv0[count] = obj.get_RV(t_rv, component = 'secondary')
					count+=1
					FWHM[count] = get_FWHM(obj.object2.star2.vsini, obj.object2.star2.BmV, spectro)
					contrast[count] = get_contrast(FWHM[count], obj.object2.star2.BmV, obj.object2.star2.z, spectro, mask)
					flux[count] = obj.object2.star2.get_flux('Johnson-V')
					rv0[count] = obj.get_RV(t_rv, component = 'tertiary')
					count+=1
				if obj.object2.type == 'plansys' : 
					FWHM[count] = get_FWHM(obj.object2.star.vsini, obj.object2.star.BmV, spectro)
					contrast[count] = get_contrast(FWHM[count], obj.object2.star.BmV, obj.object2.star.z, spectro, mask)
					flux[count] = obj.object2.star.get_flux('Johnson-V')
					rv0[count] = obj.get_RV(t_rv, component = 'secondary')
					count+=1
			if obj.type == 'plansys' : 
				FWHM[count] = get_FWHM(obj.star.vsini, obj.star.BmV, spectro)
				contrast[count] = get_contrast(FWHM[count], obj.star.BmV, obj.star.z, spectro, mask)
				flux[count] = obj.star.get_flux('Johnson-V')
				rv0[count] = obj.get_RV(t_rv, component = 'primary')
				count+=1
		
		# Make fixed CCFs
		count = 0
		CCF_sum = n.zeros(len(RV_CCF), float)
		while args[count].type == 'star' :
			CCF[count] = make_CCF(RV_CCF, rv0[count,0]-offset, FWHM[count], contrast[count], flux[count])
			CCF_sum+=CCF[count]
			count+=1	
		CCF_star = CCF_sum.copy()
		
		for i in xrange(len(t_rv)):
			CCF_sum = CCF_star.copy()
			# Make non-fixed CCFs
			for o in xrange(count, nb_star) :
				#print rv0[o,i], FWHM[o], contrast[o], flux[o]
				CCF[o] = make_CCF(RV_CCF, rv0[o,i]-offset, FWHM[o], contrast[o], flux[o])
				CCF_sum+=CCF[o]
			
			normalized_CCF[i] = n.array((CCF_sum-min(CCF_sum))/max(CCF_sum-min(CCF_sum)))
			
	return normalized_CCF[i]


