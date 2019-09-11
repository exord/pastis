"""
Module to deal with isochrones.
"""
import sys
import numpy as n
from math import log10
from scipy import interpolate

# Intra-package imports
from .exceptions import EvolTrackError, OutofIsochroneError
from .constants import Msun, Rsun, G
from .tools import loadtxt_iter


def interpol_tracks(input_file):
        
    sys.stdout.write('Reading evolution tracks file...\n')
    sys.stdout.flush()
    global z, mini, mact, logage, logT, logg, logL
    z, mini, mact, logage, logT, logg, logL = \
          loadtxt_iter(input_file, usecols = (2,3,4,5,6,7,8), unpack =True,
                  skiprows = 1)

    zu = n.unique(z)

    global miniu
    miniu = n.unique(mini)

    global minz, maxz, minminit, maxminit
    minz = min(zu); maxz = max(zu)
    minminit = min(miniu); maxminit = max(miniu)
    
    global Y
    Y = {}

    Ntotal = len(zu)*len(miniu)
     
    for i, zz in enumerate(zu):

        for j, mm in enumerate(miniu):
             
            ni = i*len(miniu) + j + 1
            nip = 1e2*float(ni)/Ntotal

            if j==0 and i==0:
                sys.stdout.write('Interpolating tracks... %02d %%'%nip)
            elif n.round(nip)%1 == 0:
                sys.stdout.write('\b'*4+'%02d %%'%nip)
            sys.stdout.flush()
             
            cond = n.logical_and(n.equal(z, zz), n.equal(mini, mm))

            # Check if there is only one node in the grid with mm or zz
             
            tt = n.compress(cond, logage)
            if len(tt) == 0: continue
            # Define normalised lifetime of star to used as interpolating
            # variable
            xx = (tt - min(tt))/(max(tt) - min(tt))
            
            Y[(mm, zz)] = []
            Y[(mm, zz)].extend([tt.min(), tt.max()])
            for y in (logT, logg, logL, mact):

                y0 = n.compress(cond, y)
                yi = interpolate.interp1d(xx, y0)
                Y[(mm, zz)].append(yi)

    global verticesY
    verticesY = n.array(list(Y.keys()))
    
    print('... DONE! \n')
    return


def prepare_tracks_target(input_file, AgeUniverse = 10.4):
        
    sys.stdout.write('Reading evolution tracks file...\n')
    sys.stdout.flush()
    
    global z, mini, mact, logage, logT, logg, logL, phase

    z, mini, mact, logage, logT, logg, logL, phase = \
	loadtxt_iter(input_file, usecols = (2,3,4,5,6,7,8, 9), unpack =True,
		  skiprows = 1)

    # Compute radius and density
    global dens
    R = n.sqrt(10**logL*(5777.0/10**logT)**4.0)
    dens = mact/R**3.0

    zu = n.sort(n.unique(z))

    global miniu
    miniu = n.unique(mini)

    global minz, maxz, minminit, maxminit
    minz = min(zu); maxz = max(zu)
    minminit = min(miniu); maxminit = max(miniu)

    global T
    T = {}

    Ntotal = (len(zu) - 1)*len(miniu)

    if 'Starevol' in input_file:
	## Only keep points from main sequence and below age of Universe
        condMS = n.logical_and(phase > 0.0, logage < AgeUniverse)

    elif 'Geneva' in input_file:
        condMS = (logage < AgeUniverse)
        
    elif 'Parsec' in input_file:
        condMS = n.logical_and(phase > 0.0, logage < AgeUniverse)
        
    for i, zz in enumerate(zu[:-1]):

        for j, mm in enumerate(miniu):
             
            ni = i*len(miniu) + j + 1
            nip = 1e2*float(ni)/Ntotal

            if j==0 and i==0:
                sys.stdout.write('Preparing tracks... %02d %%'%nip)
            elif n.round(nip)%1 == 0:
                sys.stdout.write('\b'*4+'%02d %%'%nip)
            sys.stdout.flush()
             
	    ## Keep only one z and Mini
            condTrack = n.logical_and(n.equal(z, zz), n.equal(mini, mm))
            condTrack2 = n.logical_and(n.equal(z, zu[i + 1]), n.equal(mini, mm))
            
	    ### Find condition for main sequence. This is model dependent.
            if 'Dartmouth' in input_file:

		## Only keep points from main sequence and below age of Universe
		## To do this, find for each track the maximum difference in the
		## time array. This is the moment when the MS kicks off.

                tt = n.compress(condTrack, logage)

                if len(tt) == 0:
                    continue

		# Get the time of MS beginning
                tZAMS = tt[n.argmax(n.diff(tt)) + 1]

		# Write the condition
                condMS = n.logical_and(logage >= tZAMS, logage < AgeUniverse)

            condd = n.logical_and(condMS, condTrack)
            condd2 = n.logical_and(condMS, condTrack2)
	    
            tt = n.compress(condd, logage)
            tt2 = n.compress(condd2, logage)

            if len(tt) < 2 or len(tt2) < 2:
                continue

            # Define normalised lifetime of star to used as interpolating
            # variable
            xx = (tt - min(tt))/(max(tt) - min(tt))
            xx2 = (tt2 - min(tt2))/(max(tt2) - min(tt2))

            # Concatenate life arrays
            xxc = n.sort(n.concatenate((xx, xx2)))
            
            # Concatenate age arrays
            ttc = n.sort(n.concatenate((tt, tt2)))

            T[(mm, zz, zu[i+1])] = []
            T[(mm, zz, zu[i+1])].extend([(tt.min(), tt2.min()),
                                (tt.max(), tt2.max())]
                               )

	    # Add concatenated life array to list
            T[(mm, zz, zu[i+1])].append(xxc)

            
            for y in (logT, logg, dens, mact, logL):

                y0 = n.compress(condd, y)
                y02 = n.compress(condd2, y)

                yi = interpolate.interp1d(xx, y0)(xxc)
                yi2 = interpolate.interp1d(xx2, y02)(xxc)

                T[(mm, zz, zu[i+1])].append(n.array((yi, yi2)))


    global verticesT
    verticesT = n.array(list(T.keys()))

    print('... DONE! \n')
    return


def get_stellarparams(z, logage, minit, method = 'linear'):

    ## Check that input z and minit are within absolute limits of grid
    if n.logical_or(n.greater(z, maxz), n.less(z, minz)):
        raise EvolTrackError('Metallicity out of range')
        return

    if n.logical_or(n.greater(minit, maxminit), n.less(minit, minminit)):
        raise EvolTrackError('Initial mass out of range')
        return

    ## Compute distance to _all_ nodes in the grid
    dx = (verticesY[:, 0] - minit)/minit
    dy = (verticesY[:, 1] - z)

    vertdist = dx*dx + dy*dy

    ## Evaluate if closest point in the grid corresponds to a dead star
    tmax = Y[tuple(verticesY[n.argmin(vertdist)])][1]
    tmin = Y[tuple(verticesY[n.argmin(vertdist)])][0]

    xx = (logage - tmin)/(tmax - tmin)
    
    if n.logical_or(n.less_equal(logage, tmin), n.greater_equal(logage, tmax)):
        raise EvolTrackError('Logage out of range for closest track node.')
        return

    ## Check if requested point corresponds to a node; in that case, simply
    ## return values for that node
    if min(vertdist) < 1e-10:
        vals = n.zeros((4,), 'd')
        interpfuncs = Y[tuple(verticesY[n.argmin(vertdist)])][2:]

        # Define normalised stellar lifetime
        xx = (logage - tmin)/(tmax - tmin)
        for i in range(4):
            vals[i] = interpfuncs[i](xx)
        return vals

    if method == 'linear':
        
        ## List of conditions for being in each quadrant
        condsq = []

        # Write conditions for node to be in each quadrant
        xy=dx*dy
        cond1 = n.logical_and(n.greater_equal(xy, 0), n.greater(dx, 0))
        if n.any(cond1):
            condsq.append(cond1)
        cond2 = n.logical_and(n.less(xy, 0), n.less_equal(dx, 0))
        if n.any(cond2):
            condsq.append(cond2)
        cond3 = n.logical_and(n.greater_equal(xy, 0), n.less(dx, 0))
        if n.any(cond3):
            condsq.append(cond3)
        cond4 = n.logical_and(n.less(xy, 0), n.greater_equal(dx, 0))
        if n.any(cond4):
            condsq.append(cond4)

        if len(condsq) < 3:
            raise EvolTrackError('Nodes missing from more than one quadrant.')
            return

        # List with vertices on each quadrant
        verticesq = []
        for cond in condsq:
            verticesi = n.compress(cond, verticesY, axis = 0)
            verticesq.append(verticesi)

        # List of vertices to use for interpolation
        vs = []
        for i in range(len(condsq)):
            ## Get closest grid node per quadrant
            vi = verticesq[i][n.argmin(n.compress(condsq[i], vertdist))]

            ## Check vertice correspond to existing tracks
            vertinfo = Y[tuple(vi)]
            tmin = vertinfo[0]
            tmax = vertinfo[1]
            if n.logical_or(n.less_equal(logage, tmin),
                            n.greater_equal(logage, tmax)):
                continue

            vs.append(vi)

        if len(vs) < 3:
            raise EvolTrackError('Less than three nodes.')
            return

    elif method == 'cubic':
        vs = verticesY.copy()
    
    # Array that will hold extreme track ages for vertices
    tlims0 = n.zeros((2, len(vs)), dtype='d')
    # Array that will hold parameters for vertices
    vals0 = n.zeros((4, len(vs)), dtype='d')
    # Array that will hold interpolated parameters
    vals  = n.zeros((4,), dtype='d')

    for ii, vv in enumerate(vs):
        vertinfo = Y[tuple(vv)]

        # Get tmin and tmax for each node
        tlims0[:, ii] = [vertinfo[0], vertinfo[1]]
        
    # X and Y position of nodes
    invalues = n.array(vs)

    ## Define normalised lifetime for interpolated track
    # Get interpolated tmin and tmax
    aa = n.array([minit, z]).reshape((1,2))
    
    tmini = interpolate.griddata(invalues, tlims0[0, :], aa,
                                 method = method)
    tmaxi = interpolate.griddata(invalues, tlims0[1, :], aa,
                                 method = method)

    ### USE THE SAME xx FOR ALL VERTICES!
    xx = (logage - tmini)/(tmaxi - tmini)

    # Get parameters evaluated at input age (normalised to star lifetime).
    for ii, vv in enumerate(vs):
        vertinfo = Y[tuple(vv)]
        interpfuncs = vertinfo[2:]
        for param in range(4):
            vals0[param, ii] = interpfuncs[param](xx)
        
    # Interpolate using griddata
    for ii in range(4):
        vals[ii] = interpolate.griddata(invalues, vals0[ii],
                                        (minit, z), method = method
                                        )
        if n.isnan(vals[ii]):
            raise EvolTrackError('Impossible to triangulate!')

    return vals


def get_stellarparams_target(z, y, logT, N = 4, Nt = 10, planethost = False):
    """
    Compute the stellar parameters from an input set given by either z, logg,
    and logT, or z, density, logT. The parameters are obtained by 'cleverly'
    interpolating within the evolution tracks.

    It is required that the global variables T, vertices, maxz, and minz be
    defined. This is simply done by running the function prepare_tracks_target.

    The first step is to interpolate the tracks by mass to the desired
    metallicity. This is done by interpolating tracks by pairs having the same
    mass and metallicities at both sides of the requested one. The interpolation
    is done for each star 'life' (i.e. age normalised to maximum and minimum
    age values in the track). The sampling is given by the concatenation of the
    age arrays of both original tracks.

    The second step consists in computing the distance from the desired point
    to all the points in the track, and then interpolate linearly between the
    Nt closest points of the N closest tracks

    Output: Mass, logL, logAge
    """

    # Check that input z and minit are within absolute limits of grid
    if n.logical_or(n.greater(z, maxz), n.less(z, minz)):
        raise EvolTrackError('Metallicity out of range')
        return

    nodes = []
    Tmin = []
    Tmax = []

    # Get number of parameters to interpolate
    nparams = len(T[list(T.keys())[0]][3:])

    # Start interpolation in metallicity; miniu is a global variable.
    for mm in miniu:
        zvi = n.sort(n.compress(verticesT[:, 0] == mm, verticesT[:, 1]))
        zvf = n.sort(n.compress(verticesT[:, 0] == mm, verticesT[:, 2]))

        if len(zvi) == 0 or z < n.min(zvi) or z > n.max(zvf):
            # z outside metallicity range, skip mass
            continue
        
        indz = n.searchsorted(zvi, z)

        z0 = n.array([zvi[indz - 1], zvf[indz - 1]])

        # Get factor for linear interpolation
        fmin = (z - z0[0])/(z0[1] - z0[0])
        fmax = 1 - fmin


        # Get values for the required key of dictionary
        TT = T[(mm, z0[0], z0[1])]

	# Obtain for each parameter of interest, values at the nodes
        Tmin0 = n.zeros(2)
        Tmax0 = n.zeros(2)
        logages0 = []
        lifes0 = []

        # Get concatenated life array for this pair of tracks
        xx = TT[2]
        
        # Define array that will contain intepolated tracks + life array used.
        valsmm = n.zeros((nparams + 1, len(xx)))

        # Assign interpolated values
        for ip in range(nparams):
            valsmm[ip] = TT[3+ip][0]*fmax + TT[3+ip][1]*fmin
        
        # Interpolate Minimum and Maximum age
        Tmin.append(fmax*TT[0][0] + fmin*TT[0][1])
        Tmax.append(fmax*TT[1][0] + fmin*TT[1][1])

        # Include life as last element of list
        valsmm[nparams] = xx

        nodes.append(valsmm)
        
    # If less than two tracks remain, raise Exception
    if len(nodes) < 2:
        raise EvolTrackError('Not enough tracks to do interpolation.')

    Tmin = n.array(Tmin)
    Tmax = n.array(Tmax)
    ### End of interpolation in metallicity

    
    dists = []
    ## Compute distance from point to each node
    ## Iterate for all interpolated tracks.
    for ii, node in enumerate(nodes):

	# Get input parameters for this node
        logTi = node[0, :]

        if planethost:
            rhoi = node[2, :]
            ## Compute a normalized distance
            dist2 = (1 - 10**(logTi - logT))**2 + (1 - rhoi/y)**2
        else:
            loggi = node[1, :]
            ## Compute a normalized distance
            dist2 = (1 - 10**(logTi - logT))**2 + (1 - 10**(loggi - y))**2
	
        # Keep shortest distance
        dists.append(n.min(dist2))
	
        # Get output parameter for this node
        macti = node[-2, :]
        logLi = node[-1, :]

    dists = n.array(dists)
    isort = n.argsort(dists)

    # Get only N closest tracks
    closesttracks = []
    for ii in isort[:N]:
        closesttracks.append(nodes[ii])
    Tminclosest = Tmin[isort[:N]]
    Tmaxclosest = Tmax[isort[:N]]

    ## Define array that will contain, for each node, Mact, logL, life,
    ## Tmin, Tmax of closest point in track
    vals0 = n.zeros((N*Nt, 3))
    x0 = n.zeros((N*Nt, 2))
    
    valsint = n.zeros(5)

    #Iterate over the N closest tracks
    for ii, node in enumerate(closesttracks):

	# Get input parameters for this node
        logTi = node[0, :]

        if planethost:
            yi = node[2, :]
            ## Compute a normalized distance
            dist2 = (1 - 10**(logTi - logT))**2 + (1 - yi/y)**2
        else:
            yi = node[1, :]
            ## Compute a normalized distance
            dist2 = (1 - 10**(logTi - logT))**2 + (1 - 10**(yi - y))**2
	
	# Get output parameters for this node
        macti = node[-3, :]
        logLi = node[-2, :]

	# Recover log(age) of node i
        lifei = node[-1, :]
        logagei = lifei*(Tmaxclosest[ii] - Tminclosest[ii]) + Tminclosest[ii]


	## Keep Nt closest points per track
        indt = n.argsort(dist2)[: Nt]


        ## Get node points for interpolation
        x0[ii*Nt : ii*Nt + len(indt), 0] = logTi[indt]
        x0[ii*Nt : ii*Nt + len(indt), 1] = yi[indt]

        ## Get value points for interpolation
        vals0[ii*Nt : ii*Nt + len(indt), 0] = macti[indt]
        vals0[ii*Nt : ii*Nt + len(indt), 1] = logLi[indt]
        vals0[ii*Nt : ii*Nt + len(indt), 2] = logagei[indt]

    # To deal with the case where Nt is smaller than the number of points
    # in a given track
    indices = n.argwhere(x0[:, 0] != 0)
    indices = indices.reshape(indices.shape[:-1])
    xx0  = n.take(x0, indices, axis = 0)
    vals0 = n.take(vals0, indices, axis = 0)

    ## Interpolate using all tracks
    for ip in range(vals0.shape[1]):
        valsint[ip] = interpolate.griddata(xx0, vals0[:, ip], (logT, y),
					  method = 'linear')

    if valsint[2] > n.max(Tmaxclosest) or valsint[2] < n.min(Tminclosest):
        raise EvolTrackError('Requested point is outside limits of tracks.')
    
    if n.any(n.isnan(valsint)):
        raise EvolTrackError('Impossible to interpolate!')

    return valsint[0], valsint[1], valsint[2]


def get_track(z, minit, t):
    teff = n.zeros(t.shape)
    logg = n.zeros(t.shape)
    logL = n.zeros(t.shape)
    
    for i, tt in enumerate(t):
        try:
            vals = get_stellarparams(z, tt, minit)
        except EvolTrackError:
            continue
        teff[i] = 10**vals[0]
        logg[i] = vals[1]
        logL[i] = vals[2]

    teff = n.compress(n.not_equal(teff, 0.0), teff)
    logg = n.compress(n.not_equal(logg, 0.0), logg)
    logL = n.compress(n.not_equal(logL, 0.0), logL)

    # Radius in Solar Radii
    R = n.sqrt(10**logL*(5777.0/teff)**4.0)
    #Rcgs = n.sqrt(G*minit*Msun/10**logg*1e6)
    #R = Rcgs/(Rsun*1e2)
    #Dens = 10**logg/Rcgs/(G*Msun/Rsun**3.0) # G*M/R^3
    
    return teff, logg, logL, R


def plot_maps(Y, t, istime=True):
    import pylab as p
    mini = []
    z = []
    logT = []
    logL = []
    mact = []
    logg = []
    
    for key in Y.keys():
        mini0 = key[0]
        z0 = key[1]

        if istime:
            if t > Y[key][1]: continue
            xx = (t - Y[key][0])/(Y[key][1] - Y[key][0])
        else:
            if t > 1.0:
                print('Error! Lifetime has to be < 1.')
            xx = t
        logT.append(Y[key][2](xx))
        logg.append(Y[key][3](xx))
        logL.append(Y[key][4](xx))
        mact.append(Y[key][5](xx))
        mini.append(mini0)
        z.append(z0)

    mini, z, logT, logg, logL, mact = map(n.array, (mini, z, logT, logg,
                                                    logL, mact))


    for zz, zlabel in zip((logT, logg, logL, mact),('logT', 'logg', 'logL', 'Mact')):
        f1 = p.figure()
        ax = f1.add_subplot(111)
        ax.scatter(mini, z, c = zz, edgecolor = 'none')
        ax.set_title(zlabel+'; LogAge : %.4f'%t)
        ax.set_xlabel('Minit', fontsize = 16)
        ax.set_ylabel('[Fe/H]', fontsize = 16)
        #p.colorbar()

    f = p.figure()
    ax = f.add_subplot(111)
    ax.semilogy(logg, 10**logT, 'o')

    """
    x, y, nlines = loadtxt_iter('/data/PASTIS/lib/Dartmouth/IDL/nlines.txt',
                             unpack = True)
    for i in range(len(x)):
        ax.annotate(str(int(nlines[i])), (x[i], y[i]))
    """
    p.show()
    return

    
def interpol_WD(input_file):
    """
    Interpol Bergeron's WD tables
    Teff, logg -> mass, logage
    """

    sys.stdout.write('Reading Bergeron WD tables...\n')
    sys.stdout.flush()
    
    teff, logg, mact, age = \
          loadtxt_iter(input_file, usecols = (0,1,2,27), unpack =True,
                  skiprows = 2)
                  
    logage=n.log10(age)

    invalues = n.array((teff, logg)).transpose()

    global WDmact, WDlogage
    #WDmact = interpolate.interp2d(teff,logg,mact,kind='cubic')
    #WDlogage = interpolate.interp2d(teff,logg,logage,kind='cubic')
    WDmact = interpolate.LinearNDInterpolator(invalues,mact)
    WDlogage = interpolate.LinearNDInterpolator(invalues,logage)

    return

def get_WD_mass(teff,logg):
    return WDmact(teff,logg)

def get_WD_logage(teff,logg):
    return WDlogage(teff,logg)
