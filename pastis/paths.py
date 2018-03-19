#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 11:44:57 2018

@author: rodrigo

Define useful paths and modify PYTHONPATH
"""
import os
import sys

# Read environment variables
PASTISroot = os.getenv('PASTISPATH')
assert PASTISroot is not None, 'PASTISPATH env variable not defined.'

libpath = os.getenv('PASTISLIB')
assert libpath is not None, 'PASTISLIB env variable not defined.'

idlexec = os.getenv('IDLEXEC')


# Define other useful paths
datapath = os.path.join(PASTISroot, 'datafiles')
configpath = os.path.join(PASTISroot, 'configfiles')
resultpath = os.path.join(PASTISroot, 'resultfiles')
runpath = os.path.join(PASTISroot, 'run')

# Add directory to pythonpath to read modules from fortran
fortranpath = os.path.join(os.getenv('PASTISLIB'), 'fortran')
if not fortranpath in sys.path:
    sys.path.append(fortranpath)

# Define PATHS and FILES
filterpath = os.path.join(libpath, 'Filters/')
zeromagfile = os.path.join(filterpath, 'ZeroFlux.dat')
setpath = os.path.join(libpath, 'SET')
ldpath = os.path.join(libpath, 'LD')
priorspath = os.path.join(libpath, 'Priors')
minsizefile = os.path.join(priorspath, 'MinimumSize_arcsec2.dat')

