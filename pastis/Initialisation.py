import os, sys
import numpy as n
import cPickle as pickle
from scipy import interpolate

def bin_BTsettl_spectra(redfactor = 4.0, node = ''):
    """
    Read BTsettl stellar atmospheric models.

    Flux units are in erg/s/cm^2/A.
    """

    specdir = '/data/PASTIS/lib/BT-Settl/Flux'+node

    ### Define nodes in grid to use.
    global BTteff
    BTteff0 = n.arange(4, 71)
    BTteff1 = n.arange(72, 121, 2)
    BTteff2 = n.arange(125, 200, 5)
    BTteff3 = n.arange(200, 710, 10)
    BTteff = n.concatenate([BTteff0, BTteff1, BTteff2, BTteff3])

    global BTlogg
    BTlogg = n.arange(-0.5, 6.1, 0.5)

    global BTz
    BTz1 = n.arange(-4.0, 0.1, 0.5)
    BTz = n.concatenate([BTz1, n.array([0.3, 0.5])])


    ### Define resolution array
    ww1 = n.arange(10.0, 3600.0, 0.1)
    ww2 = n.arange(3600.0, 10600.0, 0.05)
    ww3 = n.arange(10600.0, 25000.0, 0.2)
    ww4 = n.arange(2.5e4, 5.2e4, 0.5)
    ww5 = n.arange(5.2e4, 2.7e5, 10.0)
    ww6 = n.arange(2.7e5, 8.0e5, 100.0)
    ww7 = n.arange(8.0e5, 2.0e6, 1000.0)
    ww8 = n.arange(2.0e6, 1e7, 5000.0)   
    wwHR = n.concatenate((ww1, ww2, ww3, ww4, ww5, ww6, ww7, ww8))
    
    ### To define wavelength array, read smallest spectrum
    smallest = 'lte010-4.5-0.0.BT-Settl.dat.bz2'
    smallfile = os.path.join(specdir, smallest)
    
    os.system('bunzip2 %s'%smallfile) 
    f = open(smallfile[:-4])
    lines = f.readlines()
    f.close()
    os.system('bzip2 %s'%smallfile[:-4]) 

    wwLR = n.zeros((len(lines),), 'double')

    for i, ll in enumerate(lines):
        wwLR[i] = float(ll.replace('D','e').split()[0])

    ###
    # Keep only up to 300,000 A and from 2900 A
    wwLR = n.compress(n.logical_and(wwLR <= 3.0e5,wwLR >= 2900),wwLR)

    ### Define global wavelength array
    ww = wwLR[::redfactor]

    ### Define deltas for bins
    deltau = 0.5*(ww[1:] - ww[:-1])
    deltal = 0.5*(ww[:-1] - ww[1:])

    global indiu, indil
    indiu = n.zeros((len(ww),), 'int')
    indil = n.zeros((len(ww),), 'int')

    ### Compute indices of wwHR corresponding to each bin
    print('Indices')
    for jj in range(len(ww)):

        if jj == 0:
            indil[jj] = 0
            cond = n.less_equal(wwHR, ww[jj] + deltau[jj])
            indiu[jj] = n.argwhere(cond)[-1, 0]
            
        elif jj == len(ww) - 1:
            indiu[jj] = len(wwHR) - 1
            cond = n.greater(wwHR, ww[jj] + deltal[jj - 1])
            indil[jj] = n.argwhere(cond)[0, 0]

        else:
            cond = n.logical_and(n.greater(wwHR, ww[jj] + deltal[jj - 1]),
                                 n.less_equal(wwHR, ww[jj] + deltau[jj]))
            ind = n.argwhere(cond)
            indil[jj] = ind[0, 0]
            indiu[jj] = ind[-1, 0]
    print('Indices DONE')
    ###
            
    """
    ### Define array to contain all binned spectra
    BTspectra = n.zeros([len(BTz), len(BTteff), len(BTlogg), len(ww)],
                        dtype ='double')
    """
    BTspectra = {}

    ####
    # START LOOP OVER ALL SPECTA
    ####
    Ni = 0
    Ntotal = 12891.0
    # Iterate over z
    for i, zz in enumerate(BTz):

        if zz >= 0.0:
            alpha = '0.0'
        elif zz == -0.5:
            alpha = '0.2'
        elif zz <= -1.0:
            alpha = '0.4'

        # To get its sign correctly in the filename
        if zz == 0.0: zz = -0.0001

        
        #Iterate over Teff
        for j, teff in enumerate(BTteff):

            #Iterate over logg
            for k, logg in enumerate(BTlogg):

                nip = 100*Ni/Ntotal
                if i==0 and j==0 and k == 0:
                    sys.stdout.write('Progress... %02d %%'%nip)
                elif n.round(nip)%1 == 0:
                    sys.stdout.write('\b'*4+'%02d %%'%nip)
                sys.stdout.flush()

                
                # To get its sign correctly in the filename
                logg = -1*logg
                if logg == 0.0:
                    logg = -0.0001

                #if logg < 0.0: signlogg = '+'
                #else: signlogg = '-'
               
                if teff < 26:
                    filename = 'lte%03d%+.1f%+.1f.BT-Settl.dat.bz2'%(teff,
                                                                     logg,
                                                                     zz)

                else:
                    filename = 'lte%03d%+.1f%+.1fa+%s.BT-Settl.dat.bz2'%(teff,
                                                                         logg,
                                                                         zz,
                                                                         alpha)

                filepath = os.path.join(specdir, filename)

                # If file does not exist, continue
                if not os.path.exists(filepath): continue

                # Uncompress file, read it and recompress it
                Ni = Ni + 1
                os.system('bunzip2 %s'%filepath)
                f = open(filepath[:-4], 'r')
                lines = f.readlines()
                f.close()
                os.system('bzip2 %s'%filepath[:-4])
                
                #
                wave = n.zeros( (len(lines),), 'double')
                flux = n.zeros( (len(lines),), 'double')

                for ii, ll in enumerate(lines):
                    ll = ll.replace('D', 'e')
                    wave[ii] = float(ll.split()[0])
                    flux[ii] = float(ll.split()[1])

                ## Interpolate flux to the wavelengths of the HR upper envelope.
                ffHR = interpolate.interp1d(wave, flux, bounds_error = False,
                                            fill_value = 0.0)(wwHR)

                # Prepare flux array
                ff = n.zeros(len(ww))

                for jj in xrange(len(ww)):
                    ## Bin spectrum
                    ff[jj] = n.mean(ffHR[indil[jj]: indiu[jj] + 1])

                # Put mean spectrum in BTspectra dict
                #BTspectra[i, j, k] = 10**(ff - 8.0)
                if zz == -0.0001:
                    zzk = 0.0
                else:
                    zzk = round(zz, 1)
                teffk = int(100*teff)
                loggk = -round(logg, 1)
                
                BTspectra[(zzk, teffk, loggk)] = 10**(ff - 8.0)

    fout = open('/data/PASTIS/lib/BT-Settl/BTspec_%d.pickle'%redfactor, 'w')
    #fout = open('/data/PASTIS/lib/BT-Settl/BTspec_AllLambda_%d.pickle'%redfactor, 'w')
    pickle.dump([ww, BTz, 100*BTteff, BTlogg, BTspectra], fout)
    fout.close()
    return

def envolvente_BTsettl_spectra():
    """
    Read BTsettl stellar atmospheric models.

    Flux units are in erg/s/cm^2/A.
    """

    specdir = '/data/PASTIS/lib/BT-Settl/Flux'
    outdir = '/data/PASTIS/lib/BT-Settl/Flux_envolvente/'
 
    ### Define nodes in grid to use.
    BTteff0 = n.arange(4, 71)
    BTteff1 = n.arange(72, 121, 2)
    BTteff2 = n.arange(125, 200, 5)
    BTteff3 = n.arange(200, 710, 10)
    BTteff = n.concatenate([BTteff0, BTteff1, BTteff2, BTteff3])

    BTlogg = n.arange(-0.5, 6.1, 0.5)

    BTz1 = n.arange(-4.0, 0.1, 0.5)
    BTz = n.concatenate([BTz1, n.array([0.3, 0.5])])


    ### Define resolution array
    ww1 = n.arange(10.0, 3600.0, 0.1)
    ww2 = n.arange(3600.0, 10600.0, 0.05)
    ww3 = n.arange(10600.0, 25000.0, 0.2)
    ww4 = n.arange(2.5e4, 5.2e4, 0.5)
    ww5 = n.arange(5.2e4, 2.7e5, 10.0)
    ww6 = n.arange(2.7e5, 8.0e5, 100.0)
    ww7 = n.arange(8.0e5, 2.0e6, 1000.0)
    ww8 = n.arange(2.0e6, 1e7, 5000.0)   
    wwHR = n.concatenate((ww1, ww2, ww3, ww4, ww5, ww6, ww7, ww8))

    # Keep only HARPS-SOPHIE-ESPaDOns range
    # HARPS:    3780 -  6910 A
    # SOPHIE:   3872 -  6943 A
    # ESPaDOns: 3700 - 10050 A
    wwHR = n.compress(n.logical_and(wwHR >= 3700,wwHR <= 10050),wwHR)

    fout = open(outdir+'lambda.pickle','w')
    pickle.dump(wwHR,fout)
    fout.close()

    ####
    # START LOOP OVER ALL SPECTA
    ####
    Ni = 0
    Ntotal = 12891.0
    # Iterate over z
    for i, zz in enumerate(BTz):

        if zz >= 0.0:
            alpha = '0.0'
        elif zz == -0.5:
            alpha = '0.2'
        elif zz <= -1.0:
            alpha = '0.4'

        # To get its sign correctly in the filename
        if zz == 0.0: zz = -0.0001

        
        #Iterate over Teff
        for j, teff in enumerate(BTteff):

            #Iterate over logg
            for k, logg in enumerate(BTlogg):

                nip = 100*Ni/Ntotal
                if i==0 and j==0 and k == 0:
                    sys.stdout.write('Progress... %02d %%'%nip)
                elif n.round(nip)%1 == 0:
                    sys.stdout.write('\b'*4+'%02d %%'%nip)
                sys.stdout.flush()

                
                # To get its sign correctly in the filename
                logg = -1*logg
                if logg == 0.0:
                    logg = -0.0001

                #if logg < 0.0: signlogg = '+'
                #else: signlogg = '-'
               
                if teff < 26:
                    filename = 'lte%03d%+.1f%+.1f.BT-Settl.dat.bz2'%(teff,
                                                                     logg,
                                                                     zz)

                else:
                    filename = 'lte%03d%+.1f%+.1fa+%s.BT-Settl.dat.bz2'%(teff,
                                                                         logg,
                                                                         zz,
                                                                         alpha)

                filepath = os.path.join(specdir, filename)

                # If file does not exist, continue
                if not os.path.exists(filepath): continue

                # Uncompress file, read it and recompress it
                Ni = Ni + 1
                os.system('bunzip2 %s'%filepath)
                f = open(filepath[:-4], 'r')
                lines = f.readlines()
                f.close()
                os.system('bzip2 %s'%filepath[:-4])
                
                #
                wave = n.zeros( (len(lines),), 'double')
                flux = n.zeros( (len(lines),), 'double')

                for ii, ll in enumerate(lines):
                    ll = ll.replace('D', 'e')
                    wave[ii] = float(ll.split()[0])
                    flux[ii] = float(ll.split()[1])

                ## Interpolate flux to the wavelengths of the HR upper envelope.
                ffHR = interpolate.interp1d(wave, flux, bounds_error = False,
                                            fill_value = 0.0)(wwHR)

                if zz == -0.0001:
                    zzk = 0.0
                else:
                    zzk = round(zz, 1)
                teffk = int(100*teff)
                loggk = -round(logg, 1)
                
                fileout=os.path.join(outdir, 'BT'+str(zzk)+'_'+str(teffk)+'_'+str(loggk))
                #fileout=os.path.join(outdir, filename[:-4])
                
                fout = open(fileout+'.pickle','w')
                pickle.dump(10**(ffHR - 8.0),fout)
                fout.close()

    return

def WD_spectra():
    """
    Read Koester white dwarf atmospheric models.

    Flux units are in erg/s/cm^2/A.
    """

    specdir = '/data/PASTIS/lib_v1/WD_Koester/'

    ### Define nodes in grid to use.
    global WDteff
    WDteff0 = n.arange(10,12.1)*1e3
    WDteff1 = n.arange(14,19.1)*1e3
    WDteff2 = n.arange(20,40.1,2)*1e3
    WDteff = n.concatenate([WDteff0, WDteff1, WDteff2])

    global WDlogg
    WDlogg = n.arange(7.0, 9.1,0.25)

    WDspectra = {}

    ####
    # START LOOP OVER ALL SPECTRA
    ####
    Ni = 0
    Ntotal = 180.0
    # Iterate over teff
    for j, teff in enumerate(WDteff):
        
        #Iterate over logg
        for k, logg in enumerate(WDlogg):
            
            nip = 100*Ni/Ntotal
            if j==0 and k == 0:
                sys.stdout.write('Progress... %02d %%'%nip)
            elif n.round(nip)%1 == 0:
                sys.stdout.write('\b'*4+'%02d %%'%nip)
                sys.stdout.flush()

            filename = str('da'+str(long(WDteff[j])).zfill(5)+'_'+str(long((WDlogg[k]*100.))).zfill(3)+'.dk')

            filepath = os.path.join(specdir, filename)

            # If file does not exist, continue
            if not os.path.exists(filepath): continue

            # read file
            Ni = Ni + 1
            f = open(filepath, 'r')
            print(filepath)
            wave,flux = n.loadtxt(filepath, usecols = (0,1), unpack =True,skiprows = 33)
            if Ni == 1: WDww = wave
            ff = interpolate.interp1d(wave,flux, bounds_error = False,
                                      fill_value = 0.0)(WDww)

            WDspectra[(WDteff[j], WDlogg[k])] = ff/1e8 #erg/s/cm^2/cm -> erg/s/cm^2/A

    fout = open('/data/PASTIS/lib/AM/WD/WDspec.pickle', 'w')
    
    pickle.dump([WDww, WDteff, WDlogg, WDspectra], fout)
    fout.close()
    return
