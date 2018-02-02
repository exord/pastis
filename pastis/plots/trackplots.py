import sys
import numpy as n
import pylab as p

def plot_tracks(input_file, zo, **kwargs):
    """
    WRITE DOC
    """

    AgeUniverse = kwargs.pop('AgeUniverse', 10.14)
    
    sys.stdout.write('Reading evolution tracks file...\n')
    sys.stdout.flush()
    
    global z, mini, mact, logage, logT, logg, logL, phase

    z, mini, mact, logage, logT, logg, logL, phase = \
	n.loadtxt(input_file, usecols = (2,3,4,5,6,7,8, 9), unpack =True,
		  skiprows = 1)

    zu = n.unique(z)
    miniu = n.unique(mini)

    if not zo in zu:
	print('Invalid metallicity!')
	print('Valid values are:')
	for i, zz in enumerate(zu):
	    print('[%d]: %.5f'%(i, zz))
	indz = raw_input('Which one you want? [0 - %d]'%(len(zu)-1))
	zo = zu[int(indz)]

    ## Create figure
    f1 = p.figure()
    ax1 = f1.add_subplot(111)
    f2 = p.figure()
    ax2 = f2.add_subplot(111)
    f3 = p.figure()
    ax3 = f3.add_subplot(111)

    global minz, maxz, minminit, maxminit
    minz = min(zu); maxz = max(zu)
    minminit = min(miniu); maxminit = max(miniu)

    for mm in n.unique(mini):
	
	condTrack = n.logical_and(mini == mm, z == zo)
	### Find condition for main sequence. This is model dependent.
	if 'STAREVOL' in input_file or 'Geneva' in input_file:
	    ## Only keep points from main sequence and below age of Universe
	    condMS = n.logical_and(phase > 0.0, logage < AgeUniverse)
	
	if 'Dartmouth' in input_file:
	
	    ## Only keep points from main sequence and below age of Universe
	    ## To do this, find for each track the maximum difference in the
	    ## time array. This is the moment when the MS kicks off.
	
	    tt = n.compress(condTrack, logage)
	
	    if len(tt) == 0:
		continue
	
	    # Get the time of MS beginning.
	    tZAMS = tt[n.argmax(n.diff(tt)) + 1]

	    # Write the condition
	    condMS = n.logical_and(logage >= tZAMS, logage < AgeUniverse)

	condd = n.logical_and(condMS, condTrack)

	tt = n.compress(condd, logage)
	if len(tt) < 2:
	    continue

	logTi = n.compress(condd, logT)
	loggi = n.compress(condd, logg)
	macti = n.compress(condd, mact)
	logLi = n.compress(condd, logL)

	ax1.plot(10**logTi, loggi,  **kwargs)
	ax1.annotate(str(mm), (10**logTi[0], loggi[0]), ha = 'left')

	# Compute radius
	Ri = n.sqrt(10**logLi*(5777.0/10**logTi)**4.0)

	ax2.plot(Ri, logLi,  **kwargs)
	ax2.annotate(str(mm), (Ri[0], logLi[0]), ha = 'left')
	
	# Compute density
	densi = macti/Ri**3
	
	ax3.plot(10**logTi, densi,  **kwargs)
	ax3.annotate(str(mm), (10**logTi[0], densi[0]), ha = 'left')

    for aa in (ax1, ax3):
	aa.set_xlabel('Effective Temperature [K]', fontsize = 16)

    ax2.set_xlabel('Radius [solar radii]', fontsize = 16)

    ax1.set_ylabel('logg [cgs]', fontsize = 16)
    ax2.set_ylabel('logL [solar unit]', fontsize = 16)
    ax3.set_ylabel('Density [solar units]', fontsize = 16)

    ax1.invert_xaxis()
    ax3.invert_xaxis()
    ax1.invert_yaxis()

    return
