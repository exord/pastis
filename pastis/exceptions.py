#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 11:34:35 2018

@author: rodrigo

Define Exception classes
"""
class EvolTrackError(Exception):
    pass


class OutofIsochroneError(Exception):
    pass


class SpectrumInterpolError(Exception):
    def __init__(self, mesg, indz, indteff, indlogg, outside=False):
        self.mesg = mesg
        self.indz = indz
        self.indteff = indteff
        self.indlogg = indlogg
        self.outside = outside

    def __str__(self):
        return self.mesg

"""
class SpectrumGridError(Exception):
    def __init__(self, mesg, indz, indteff, indlogg):
        self.mesg = mesg
        self.indz = indz
        self.indteff = indteff
        self.indlogg = indlogg
"""

class EBOPparamError(Exception):
    pass


class GlobalSpectrumError(Exception):
    pass


__all__ = ['EvolTrackError', 'SpectrumInterpolError']