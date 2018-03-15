#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 11:40:35 2018

@author: rodrigo

Define numerical constants.
"""
# Define physical constants
G = 6.67428e-11  # m3 kg-1 s-2
c = 299792458.0  # m/s
SB = 5.67051e-8  # W m-2 K-4

# Astronomical constants
pc = 3.08568025e16  # parsec [m]
# au = 149597870691  # astronomical unit [m]
au = 149597870700  # New definition (IAU 2012) of the astronomical unit [m]
Msun = 1.98842e30  # kg
Mjup = 1.89852e27  # kg
Mearth = 5.9736e24  # kg
Rsun = 6.95508e8  # equatorial radius [m]
Rjup = 71492000.  # equatorial radius [m]
Rearth = 6378137.0  # equatorial radius [m]
R_V = 3.1  # normal extinction of the galaxy

# Bolometric solar magnitude for prior computation
Mbolsun = 4.74

# Define useful unit conversions
# Msun2Mjup = 1.047348644e3 # From IAU
# (http://maia.usno.navy.mil/NSFA/NSFA_cbe.html)
Mjup2Mearth = Mjup/Mearth
Mearth2Mjup = Mearth/Mjup
Msun2Mjup = Msun/Mjup
Mjup2Msun = Mjup/Msun
Rjup2Rearth = Rjup/Rearth
Rearth2Rjup = Rearth/Rjup
Rsun2Rjup = Rsun/Rjup
Rjup2Rsun = Rjup/Rsun
pc2Rsun = pc/Rsun
