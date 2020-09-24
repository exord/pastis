#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 00:19:12 2020

@author: rodrigo

Prepare a Kernel Density Estimate of the orbital period,
transit depth and duration to use when simulating light curves.
"""
import numpy as np
import pandas as pd
import scipy.stats as st

# home = os.getenv('HOME')

def prepare_toi_dist(toicsv):

    # Read TOI list
    dd = pd.read_csv(toicsv, skiprows=4)
    
    # Read values and compute log (much nicer distributions)
    lp = np.log10(dd['Orbital Period Value'])
    ld = np.log10(dd['Transit Depth Value'])
    ldur = np.log10(dd['Transit Duration Value'])

    ii = ~pd.isnull(lp) & ~pd.isnull(ld) & ~pd.isnull(ldur)
    
    return st.kde.gaussian_kde([lp[ii], ld[ii], ldur[ii]])

if __name__ == '__main__':
    kk = prepare_toi_dist('/Users/rodrigo/ExP/pastisML/'
                          'csv-file-toi-catalog.csv')
    