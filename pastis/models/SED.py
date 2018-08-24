import numpy as np

from ..AstroClasses import Star, IsoBinary, PlanSys, Triple
from .. import photometry as phot
from ..exceptions import GlobalSpectrumError


def PASTIS_SED(photbands, *args):
    """
    Compute SED for all objects in args, for all photbands.
    """

    ### COMMENTED WHEN TRANSLATED TO PACKAGE ####
    ### global_spectrum WAS NOT DELETED IN PASTIS_MCMC ###
    """
    if 'global_spectrum' in globals():
        pass
    else: compute_global_spectrum(*args)
    """

    compute_global_spectrum(*args)
    mags = []
    mags = [phot.flux2mag(phot.get_flux(global_spectrum, pb), pb) \
            for pb in photbands]
#    for pb in photbands:
        #mags[pb] = flux2mag(get_flux(spectra, pb), pb)
#        mags.append(phot.flux2mag(phot.get_flux(global_spectrum, pb), pb))

    return mags

def compute_global_spectrum(*args):
    spectra = 0.0
    for obj in args:
        if isinstance(obj, Star):
            # Get spectrum from star
            spectra = spectra + obj.get_spectrum()

        if isinstance(obj, IsoBinary):
            # Add spectrum of each component
            spectra = spectra + obj.star1.get_spectrum()
            spectra = spectra + obj.star2.get_spectrum()

        if isinstance(obj, PlanSys):
            # Get spectrum from host star
            spectra = spectra + obj.star.get_spectrum()
            

        if isinstance(obj, Triple):
            # Get LC for each component of triple system
            for component in (obj.object1, obj.object2):
                if isinstance(component, Star):
                    spectra = spectra + component.get_spectrum()
                elif isinstance(component, IsoBinary):
                    spectra = spectra + component.star1.get_spectrum()
                    spectra = spectra + component.star2.get_spectrum()
                elif isinstance(component, PlanSys):
                    spectra = spectra + component.star.get_spectrum()

    global global_spectrum
    global_spectrum = spectra

    if np.max(global_spectrum) == 0:
        raise GlobalSpectrumError('Global Spectrum is zero!')

    return 
