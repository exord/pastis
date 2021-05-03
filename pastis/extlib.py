#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 12:01:32 2018

@author: rodrigo

Define dictionaries redirecting to external files in the PASTIS lib.
"""
import os

from .paths import ldpath, setpath, priorspath

# Dictionaries
LDdict = {'Claret2011': os.path.join(ldpath, 'claret2011ab.txt'),
          'Claret2011_wTESS': os.path.join(ldpath, 'claret2011abTESS.txt')}

EMdict = {'Dartmouth': os.path.join(setpath, 'Dartmouth.trk'),
          'Geneva': os.path.join(setpath, 'Geneva.trk'),
          'StarEvol': os.path.join(setpath, 'Starevol.trk'),
          'Parsec': os.path.join(setpath, 'Parsec.trk')}

SAMdict = {'CastelliKurucz': 'CK', 'BT-settl': 'BT'}

PRIORdict = {'Mayor': [os.path.join(priorspath, 'PDF_SmallPlanets_Mayor.dat'),
                       os.path.join(priorspath, 'CDF_GiantPlanets_Mayor.dat')],
             'Fressin': [os.path.join(priorspath,
                                      'Fressin_cumulative_occurrence.txt'),
                         os.path.join(priorspath,
                                      'Fressin_period_bin_edges.txt')]}
# GIEdict = {'AmoresLepine2005': ''}
