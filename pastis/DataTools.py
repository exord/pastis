import numpy as n


def oversample_time_notexp(v, oversampling_factor, dt=None):
    """
    Given a sorted array v, produces a oversampled vector with
    oversampling_factor elements for each element in v
    
    Inputs:
    -------
    
    v: list array
    Array to oversample
    
    oversampling_factor: float or list or array
    Oversampling factor
    
    dt: float
    no idea...
    
    Output:
    -------
    
    t: array
    oversampled array of v
    
    """
    v = n.array(v, float)
    if len(v) == 1:
        return v
    if dt is None:
        dt = v[1] - v[0]

    """
    ## Produce dt from mode of v difference
    dv = n.diff(v)
    
    # Compute reasonable size of bin
    mindiff = n.min(dv)
    bbins = n.arange(n.min(dv), n.max(dv) + mindiff/5000.0, mindiff/5000.0)
    hist, bins_edges = n.histogram(dv, bbins)
    
    indmax = n.argmax(hist)
    dt = 0.5*(bins_edges[indmax] + bins_edges[indmax + 1.0])
    """
    if oversampling_factor == 1:
        return v
    
    elif n.isscalar(oversampling_factor):
        oversampling_factor = n.zeros(len(v), float) + oversampling_factor

    t = n.array([])
    for jj in range(len(v)):
        if dt is None:
            dt = v[jj+1] - v[jj]
        ti = n.arange(v[jj] - dt * (oversampling_factor[jj] - 1) /
                      (2 * oversampling_factor[jj]),
                      v[jj] + dt * (oversampling_factor[jj] - 1) /
                      (2 * oversampling_factor[jj]) +
                      dt / (2 * oversampling_factor[jj]),
                      dt / oversampling_factor[jj])
        t = n.concatenate((t, ti))

    return t


def oversample_time(t, texp=None, dt=None, oversampling_factor=None):
    """
    Given a sorted array t with texp time integration, produces an oversampled
    vector with dt time-fixed elements for each element in t
    
    Inputs
    ------
    
    t: array
    Input time vector. The time defined in this vector should be the center of
    the exposure.
    
    texp: float, list or array
    Integration time (or exposure time) for each elements of the vector t
    (in seconds).
    
    dt: float
    Oversampling timescale. dt should be equal or smaller than any elements
    of texp (in seconds).
    
    oversampling_factor: integer
    Oversampling factor
    
    Output
    ------
    
    tt: array
    oversampled time vector with timescale of dt or by a factor of
    oversampling_factor
    """
    t = n.array(t, float)
    if texp is not None and n.isscalar(texp):
        texp = n.zeros(len(t), float) + texp
    elif texp is not None and len(texp) != len(t):
        raise IOError('Inputs arguments t and texp should have the same length')
    elif texp is not None:
        texp = n.array(texp, float)
    elif texp is None and oversampling_factor is not None:
        return oversample_time_notexp(t, oversampling_factor)
    else:
        raise IOError('Input argument texp or oversampling_factor should be at '
                      'least provided')

    if oversampling_factor is not None and n.isscalar(oversampling_factor):
        oversampling_factor = n.zeros(len(texp), int) + oversampling_factor
    elif oversampling_factor is not None and len(oversampling_factor) != len(t):
        raise IOError('Inputs arguments t and oversampling_factor should have '
                      'the same length')
    elif oversampling_factor is not None:
        oversampling_factor = n.array(oversampling_factor, int)

    if oversampling_factor is not None and any(oversampling_factor) < 1:
        raise IOError('Any elements of input argument oversampling_factor '
                      'should be larger than 1')
    elif dt is not None and dt >= any(texp):
        raise IOError('Oversampling timescale dt ({0}) should be equal or '
                      'smaller than any elements of '
                      'texp ({1})'.format(dt, n.min(texp)))

    if oversampling_factor is None and dt is None:
        raise IOError('Either oversampling factor or oversampling timescale dt '
                      'should be provided')
    elif oversampling_factor is not None and dt is not None:
        raise IOError('Oversampling factor and oversampling timescale dt '
                      'should not be provided together. Choose one !')

    tt = []
    if dt is not None:
        for i in range(len(t)):
            n_sample = int(round(texp[i] / dt))
            if n_sample < 1:
                raise ValueError('Oversampling timescale (dt) should be '
                                 'equal or smaller than any elements of texp')
            elif n_sample < 2:
                tt.append(t[i])
            else:
                tt.append(n.linspace(t[i] - texp[i] / 2. / 86400.,
                                     t[i] + texp[i] / 2. / 86400., n_sample))
    elif oversampling_factor is not None:
        for i in range(len(t)):
            if oversampling_factor[i] < 1:
                raise ValueError('Oversampling factor should be at least 1.')
            elif oversampling_factor[i] == 1:
                tt.append(t[i])
            else:
                tt.append(n.linspace(t[i]-texp[i]/2./86400.,
                                     t[i]+texp[i]/2./86400.,
                                     oversampling_factor[i]))
    return n.hstack(tt)


def readdata(input_dict):
    """
    Read data from files as indicated in datadict.

    Return dictionary with same keys, continaing data for each instrument.
    """
    global datadict
    datadict = input_dict.copy()

    global lightcurves
    lightcurves = []

    for key in datadict.keys():
        # By default, don't check oversampling
        checksampling = False

        if datadict[key]['type'] == 'PHOT':
            lightcurves.append(key)
            # If photometric, check if oversampling need to be done (but later)
            checksampling = True

        if datadict[key]['type'] == 'RV':
            checksampling = True

        f = open(datadict[key]['datafile'], 'r')
        lines = f.readlines()
        f.close()

        # Read header from datafile
        header = lines[0].split()
        fmtall = dict((header[i], i) for i in range(len(header)))

        # Prepare arrays for data
        datadict[key]['data'] = {}

        for line in lines[2:]:
            if line.startswith('#') or len(line) < 2:
                continue

            for fmt in fmtall.keys():
                elem = line.split()[fmtall[fmt]]
                try:
                    elem = float(elem)
                except ValueError:
                    pass

                if fmt in datadict[key]['data']:
                    datadict[key]['data'][fmt].append(elem)
                else:
                    datadict[key]['data'][fmt] = [elem, ]

        for dd in datadict[key]['data']:
            datadict[key]['data'][dd] = n.array(datadict[key]['data'][dd])

        if checksampling:
            if 'texp' in datadict[key]:
                texp = datadict[key]['texp']
            elif 'texp' in datadict[key]['data']:
                texp = datadict[key]['data']['texp']
            else:
                texp = None

            if 'sampling' in datadict[key]:
                sampling = datadict[key]['sampling']
            else:
                sampling = None

            if 'dt' in datadict[key]:
                dt = datadict[key]['dt']
            else:
                dt = None

            if sampling is None and dt is None and texp is None:
                sampling = 1

            # Generate oversampled time array and save it in dictionary
            tt = oversample_time(datadict[key]['data']['time'],
                                 oversampling_factor=sampling, texp=texp,
                                 dt=dt)

            datadict[key]['overtime'] = tt
    return datadict, lightcurves