#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 11:34:35 2018

@author: rodrigo

Define Exception classes
"""


class EvolTrackError(BaseException):
    pass


class OutofIsochroneError(BaseException):
    pass


class SpectrumInterpolError(BaseException):
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


class EBOPparamError(BaseException):
    pass


class GlobalSpectrumError(BaseException):
    pass


__all__ = ['EvolTrackError', 'SpectrumInterpolError']
