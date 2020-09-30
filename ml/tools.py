#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 23:50:39 2020

@author: rodrigo
"""
import numpy as np

def invert_duration(dur_hour, per_day, ecc, omega_deg,
                    kr, b):
    """
    The inverse of planet.transit_duration_ecc for a/R*
    """
    
    # Compute auxiliary variables
    ume2 = 1 - ecc**2
    umco = 1 + ecc * np.cos(np.pi/2 - omega_deg * np.pi / 180.0)
    tau_hour = per_day * 24 / np.pi
    
    ar2 = (tau_hour/dur_hour)**2 * (ume2/umco) * ((1 + kr)**2 - b**2)
        
    return np.sqrt(ar2)
    
def transit_duration_ecc(per_day, ecc, omega_deg, kr, ar, b):
    """
    Returns the duration (in hours) of a transit for a planet in an 
    eccentric orbit, under the assumptions in Tingley & Sackett (2005).
    
    The input parameters can be either float numbers or arrays
    
    :param float or np.array per_day: orbital period in days
    :param float or np.array ecc: orbital eccentricity
    :param float or np.array omega_deg: argument of pericenter in degrees
    :param float or np.array kr: radius ratio, Rp/Rs
    :param float or np.array ar: normalised semi-major axis, a/Rs
    :param float or np.array b: impact parameter, rt/a * cos(incl)
    
    """
    # Compute distance at inferior conjunction
    rt = r_infconj(ecc, omega_deg, ar)
    cosi = b/rt
    # Useful renames
    per_hour = per_day * 24.0
    tau = per_hour / np.pi
    
    ume2 = 1 - ecc**2
    umco = 1 + ecc * np.cos(np.pi/2 - omega_deg * np.pi / 180.0)
    rho = ar/(1 + kr)
    
    d_dcirc = np.sqrt(ume2)/umco * np.sqrt(1 - ( rho * ume2 * cosi / umco)**2 )
    dcirc = tau / rho
    return dcirc * d_dcirc
