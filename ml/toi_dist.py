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

class TOIdist(object):
    
    def __init__(self, toicsv):
        
        self.csvfile = toicsv
        
        # Read TOI list
        dd = pd.read_csv(toicsv, skiprows=4)
            
        self.toitable = dd
        
        # Read values and compute log (much nicer distributions)
        lp = np.log10(self.toitable['Orbital Period Value'])
        ld = np.log10(self.toitable['Transit Depth Value'])
        ldur = np.log10(self.toitable['Transit Duration Value'])
        
        # Keep all columns without nulls
        ii = ~pd.isnull(lp) & ~pd.isnull(ld) & ~pd.isnull(ldur)
    
        self.kde = st.kde.gaussian_kde([lp[ii], ld[ii], ldur[ii]])
        
        return
    
    def sample(self, size=1):
        """
        Return sample from toi distribution
        """
        
        assert self.kde is not None, "Cannot sample from uninstanciated"
        
        x = self.kde.resample(size=size)
        
        xdict = dict((['logP', x[0]], ['logDepth', x[1]], 
                      ['logDur', x[2]]))
        
        return x, xdict


if __name__ == '__main__':
    toidist = TOIdist('/Users/rodrigo/ExP/pastisML/'
                      'csv-file-toi-catalog.csv')
    