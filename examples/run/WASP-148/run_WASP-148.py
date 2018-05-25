objectname = 'WASP-148'
version = 'pastis'
submit = False

PASTIS = __import__(version+'.run')

datapath = __import__(version).paths.datapath
configpath = __import__(version).paths.configpath

import os
import pickle
import numpy as n

def create_pastis_file(commentstr):

    infodict = {'name': objectname,
        'comment': commentstr,
                #'EvolModel' : 'Parsec',
                'EvolModel': 'Dartmouth',
                #'EvolModel': 'Geneva',
                #'EvolModel': 'StarEvol',
        'SAM': 'BT-settl',
        'LDC': 'Claret2011',
        'email': 'rodrigo.diaz@oamp.fr',
                #
        'alpha': 284.87997,
        'delta': 49.26699,
        'EXTstep': 10.0,
        'MaxDist': 10000.0,
        #
        'PCA': 1,
        'Min_PCA': 10000,
        'BI_PCA': 5000,
        'Max_PCA': 100000000,
        'N_update_PCA': 5000,
        #
        'beta': 1.0,
        'Nbeta': 1,
        'Nchain': 1,
        'Nmcmc': int(5e4),
        'Nsave': int(5e5),
        'randomstart': 1,
        'beaming': False,
        'Nmax': int(1e4),
        'walltime' : 100
        }

    ddir = os.path.join(datapath, objectname)
    target = objectname.replace('-','')
    
    datafile_list = {'WASP1': {'MeanFlux': 0.0,
                    'datafile': os.path.join(ddir,
                                            '{}_SuperWASP1.rdb'.format(target)),
                    'filter': 'Johnson-V',
                    'is_phase': 0,
                    'sampling': 1,
                    'type': 'PHOT'
                    },

                    'WASP2': {'MeanFlux': 0.0,
                    'datafile': os.path.join(ddir,
                                             '{}_SuperWASP2.rdb'.format(target)),
                    'filter': 'Johnson-V',
                    'is_phase': 0,
                    'sampling': 1,
                    'type': 'PHOT'
                    },
                
                    'SOPHIEp': {'datafile': os.path.join(ddir,
                                                        '{}_SOPHIE2.rdb'.format(target)),
                                'instrument': 'SOPHIE HE',
                                'mask': 'G2',
                                'is_BIS': 0,
                                'is_CONT': 0,
                                'is_FWHM': 0,
                                'is_RV': 1,
                                'type': 'RV'
                                },

                    'SOPHIE': {'datafile': os.path.join(ddir,
                                                        '{}_SOPHIE1.rdb'.format(target)),
                                'instrument': 'SOPHIE HE',
                                'mask': 'G2',
                                'is_BIS': 0,
                                'is_CONT': 0,
                                'is_FWHM': 0,
                                'is_RV': 1,
                                'type': 'RV'
                                }

                }
        

    #---------------------
    #        OBJECTS
    #--------------------
    hostdict = {
        'dens': [0.0, 1, 'Normal', 1.79 , 0.3, 0.00, 0.0, ''],
        'teff': [4850.0, 1, 'Normal', 4850.0, 100.0, 0.0, 0.0, ''],
        'z': [-0.07, 1, 'Normal', -0.07, 0.18, 0.0, 0.0, ''],
        'dist': [1000.0, 1, 'Uniform', 100.0, 2000.0, 0.0, 0.0, ''],
        'ebmv': [0.1, 1, 'Uniform', 0.0, 0.2, 0.0, 0.0, '', None],
        'v0': [-74.0, 1, 'Uniform', -78.0, -70.0, 0.0, 0.0, ''],
        'vsini': [10.0, 0, 'Normal', 0.0, 0.0, 0.0, 0.0, ''],
        'gd': [1.0, 0, 'Normal',0.0,0.0,0.0,0.0,'']
        }

    fp1dict = {'P': [8.803611, 2, 'Normal', 8.803611, 0.014, 0.0, 0.0, ''],
               'T0': [53133.0505, 1, 'Normal', 53133.0505, 0.015, 0.0, 0.0, ''], 
               'K1': [6.0, 1, 'Uniform', 0.0, 0.5, 0.0, 0.0, ''],
               'albedo2': [0.0, 0, 'Normal', 0.0, 0.0, 0.0, 0.0, ''],
               'ecc': [0.2, 1, 'Uniform', 0.0, 1.0, 0.0, 0.0, ''],
               'b' : [0.5, 1, 'Uniform', 0.0, 1.0, 0.0, 0.0, ''],
               # 'incl' : [90.0, 1, 'Sine', 60.0, 90.0, 0.0, 0.0, ''],
               'kr': [0.1, 1, 'Jeffreys', 0.01, 0.5, 0.0, 0.0, ''],
               'ar': [5, 1, 'Jeffreys', 5.0, 20.0, 0.0, 0.0, ''],
               'omega': [240.7, 1, 'Uniform', 0.0, 360.0, 0.0, 0.0, ''],
               'ua1': [0, 0, 'Uniform',  0.0, 360.0, 0.0, 0.0, ''],
               'ub1': [0, 0, 'Uniform',  0.0, 360.0, 0.0, 0.0, ''],
               'gd1': [0, 0, 'Uniform',  0.0, 360.0, 0.0, 0.0, ''],
               'v0': [0, 0, 'Uniform',  0.0, 360.0, 0.0, 0.0, ''],
               'q': [0.0, 0, 'Uniform', 0.0, 0.2, 0.0, 0.0, '']
              }

    fp2dict = {'P': [34.737, 2, 'Normal', 30.360447, 5.0, 0.0, 0.0,
                    ''],
              'T0': [54981.09095, 1, 'Normal', 54981.09095, 0.0015,
                     0.0, 0.0, ''],
              'K1': [6.0, 1, 'Uniform', 0.0, 0.5, 0.0, 0.0, ''],
              'albedo2': [0.0, 0, 'Normal', 0.0, 0.0, 0.0, 0.0, ''],
              'ecc': [0.2, 1, 'Uniform', 0.0, 1.0, 0.0, 0.0, ''],
              'b' : [0.5, 1, 'Uniform', 0.0, 1.0, 0.0, 0.0, ''],
              # 'incl' : [90.0, 1, 'Sine', 60.0, 90.0, 0.0, 0.0, ''],
              'kr': [0.1, 1, 'Jeffreys', 0.01, 0.5, 0.0, 0.0, ''],
              'ar': [5, 1, 'Jeffreys', 5.0, 20.0, 0.0, 0.0, ''],
              'omega': [240.7, 1, 'Uniform', 0.0, 360.0, 0.0, 0.0, ''],
              'ua1': [0, 0, 'Uniform',  0.0, 360.0, 0.0, 0.0, ''],
              'ub1': [0, 0, 'Uniform',  0.0, 360.0, 0.0, 0.0, ''],
              'gd1': [0, 0, 'Uniform',  0.0, 360.0, 0.0, 0.0, ''],
              'v0': [0, 0, 'Uniform',  0.0, 360.0, 0.0, 0.0, ''],
              'q': [0.0, 0, 'Uniform', 0.0, 0.2, 0.0, 0.0, '']
          }

    waspdict = {'contamination': [0.0, 0, 'Normal', 0.0, 0.0, 0.0, 0.0, ''],
                'foot': [1.0, 1, 'Uniform', 0.98, 1.02, 0.0, 0.0, ''],
                'jitter': [0.0001, 1, 'Uniform', 0.0, 0.08, 0.0, 0.0, '']
               }

    sophiedict = {'RVoffset': [0.0, 1, 'Uniform', 4.0, 6.0, 0.0, 0.0, ''],
                  'RVjitter': [0.01, 1, 'Uniform', 0.0, 0.120, 0.0, 0.0, '']
                 }

    sophiepdict = {'RVoffset': [0.0, 1, 'Uniform', 4.0, 6.0, 0.0, 0.0, ''],
                   'RVjitter': [0.01, 1, 'Uniform', 0.0, 0.120, 0.0, 0.0, '']
                  }


    list_objects = {'WASP1' : waspdict,
                    'WASP2' : waspdict.copy(),
                    #'OVERSKY' : seddict,
                    'SOPHIE': sophiedict,
                    'SOPHIEp': sophiepdict,
                    'FitPlanet1' : fp1dict,
                    'FitPlanet2' : fp2dict}

    list_prior = {}

    pastisfile = os.path.join(configpath, infodict['name'],
                  infodict['name']+'_'+commentstr+'.pastis'
                  )

    f = open(pastisfile, 'wb')
    pickle.dump([infodict,list_objects,datafile_list,list_prior], f)
    f.close()

    return pastisfile

if __name__ == '__main__':
    print('STARTING ITERATION. Fingers crossed!')

    ## Iteration over these files
    pastisfile = create_pastis_file('test')

    PASTIS.run.run_sim(pastisfile, pastisversion = version,
                       submit = submit)

           
