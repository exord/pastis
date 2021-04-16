from math import pi, log10
import numpy as n

# Intra-package imports
# from . import *
from . import TrefRV, Nmax, beaming, spotmodel
from .constants import G, c, Msun, Rsun, pc2Rsun, Mjup2Msun, Rjup2Rsun, Mjup
from .exceptions import SpectrumInterpolError
from .models.PHOT import run_EBOP
from .models.PHOT import macula
from . import tools
from . import photometry as phot

try:
    from .photometry import ww
except ImportError:
    pass
    # Commented (see issue #334)
    print('Could not import stellar atmosphere models. Are you only using '
          'FitObs objects?')

from . import isochrones as iso
from . import limbdarkening as ldark
from . import gravitydarkening as gdark
from . import albedo as albedo

try:
    from .extinction import EBmVfunc
except ImportError:
    # Commented (see issue #334)
    # print('Could not load extinction functions. Will E(B-V) be fitted?')
    # In this case, Absorption might be fitted.
    pass

from .extinction import ext_law

quiet = True


class Star(object):
    """
    Star

    age: log10(Age [yr])
    z: Metallicity ([M/H] = log(Z/Zsun); abundance of metals with respect to H;
    Sun = 0.019)
    Minit: Initial mass [solar masses]

    """

    def __init__(self, **kwargs):

        # Stellar track parameters
        self.minit = kwargs.pop('minit', None)
        self.z = kwargs.pop('z', None)
        self.logage = kwargs.pop('logage', None)
        self.age = kwargs.pop('age', None)
        if self.logage is None and self.age is not None:
            self.logage = n.log10(self.age * 1e9)

        self.dens = kwargs.pop('dens', None)

        # Stellar atmospheric parameters
        self.teff = kwargs.pop('teff', None)
        self.logg = kwargs.pop('logg', None)
        self.zeta = kwargs.pop('zeta', 0.)

        # Direct parameters
        self.mact = kwargs.pop('mact', None)
        self.R = kwargs.pop('R', None)

        # Other parameters
        self.dist = kwargs.pop('dist', 0.001)
        self.vsini = kwargs.pop('vsini', 10.0)
        self.v0 = kwargs.pop('v0', 0.0)
        self.drift1 = kwargs.pop('drift1', 0.0)
        self.drift2 = kwargs.pop('drift2', 0.0)
        self.drift3 = kwargs.pop('drift3', 0.0)
        self.ebmv = kwargs.pop('ebmv', None)

        # Wavelength dependent parameters
        self.albedo = kwargs.pop('albedo', None)
        self.ua = kwargs.pop('ua', None)  # Linear coefficient LD
        self.ub = kwargs.pop('ub', None)  # Quad coefficient LD
        self.gd = kwargs.pop('gd', None)  # GravityDarkening
        self.B = kwargs.pop('B', None)  # Beaming effect
        self.alphaS = kwargs.pop('alphaS',
                                 0.0)  # Effective filling factor due to spots

        if spotmodel == 'Macula':
            # Parameters for spot modeling
            for spotstarparams in ['rotangle', 'period', 'kappa2', 'kappa4',
                                   'c1', 'c2', 'c3', 'c4', 'd1', 'd2', 'd3',
                                   'd4']:
                self.__setattr__(spotstarparams,
                                 kwargs.pop(spotstarparams, None))

            if self.c1 is None or self.c2 is None or self.c3 is None or \
                    self.c4 is None:
                if self.ua is not None and self.ub is not None:
                    # Impose quadratic limb-darkening from star
                    self.c1 = 0
                    self.c2 = self.ua + 2 * self.ub
                    self.c3 = 0
                    self.c4 = -self.ub

            # Number of spots
            self.Nspots = 0
            for key in kwargs.keys():
                if 'life' in key:
                    self.Nspots += 1

            if self.Nspots != 0:
                for i in range(self.Nspots):
                    for spotparam in ['lambda0', 'phi0', 'alphamax', 'fspot',
                                      'tmax', 'life', 'ingress', 'egress']:
                        self.__setattr__('spot' + str(i) + spotparam,
                                         kwargs.pop('spot'+str(i)+spotparam,
                                                    None))

        # To avoid nan in spectrum
        if self.dist == 0.0:
            self.dist = 0.001

        self._parent = None

        # Check that critical parameter z was provided
        self.test_critical(['z'], ['Metallicity'])

    def test_critical(self, attrs, names):

        missing = []
        for critical, name in zip(attrs, names):
            # Get value of critical parameter
            x = getattr(self, critical)

            # If not defined add to missing list
            if x is None:
                missing.append([critical, name])

        # If any critical parameter is missing, print info and exit.
        if len(missing) != 0:
            s = ['{} ({}) not provided'.format(nam, cri) for
                 cri, nam in missing]

            raise ValueError('Undefined critical arguments\n '
                             ''+str.join('\n ', s))

        return

    def get_parent(self):
        return self._parent

    def get_spectrum(self, force=True):
        """
        Get the spectrum of the star

        Needs global variable: CKspectra, CKz, CKteff, CKlogg
        """
        if not force:
            try:
                return self.spectrum
            except AttributeError:
                pass

        teff = self.teff
        logg = self.logg
        z = self.z

        # Extinction
        if self.ebmv is None:
            # Extinction comes from a model
            self.ebmv = self.get_extinction()

        try:
            spectrum = phot.get_interpolated_AM(z, teff, logg)
        except SpectrumInterpolError:
            spectrum = phot.get_nearest_AM(z, teff, logg)

        # Correct from distance, extinction, and radius
        spectrum = spectrum * (self.R / (self.dist * pc2Rsun)) ** 2
        spectrum = spectrum * 10 ** (-0.4 * ext_law(ww) * self.ebmv)
        self.spectrum = spectrum

        return spectrum

    def get_flux(self, photband):
        """
        Get flux of star in a given photometric band
        """
        spectrum = self.get_spectrum()

        flux = phot.get_flux(spectrum, photband)
        return flux

    def get_mag(self, photband):
        """
        Get magnitude of Star in photometric band photband.
        """
        spectrum = self.get_spectrum()

        flux = phot.get_flux(spectrum, photband)
        mag = phot.flux2mag(flux, photband)

        return mag

    def get_LD(self, photband):
        """
        Return quadratic limb-darkening parameters in a given photometric band.
        """
        if (isinstance(self.ua, float) or isinstance(self.ua, int)) and \
                (isinstance(self.ub, float) or isinstance(self.ub, int)):
            # Return same coefficients independetly of photband demanded
            return n.array([self.ua, self.ub])

        elif isinstance(self.ua, dict) and isinstance(self.ub, dict):
            try:
                uai, ubi = self.ua[photband], self.ub[photband]
                if (isinstance(uai, float) or isinstance(uai, int)) and \
                        (isinstance(ubi, float) or isinstance(ubi, int)):
                    return n.array([uai, ubi])
                else:
                    ldark.get_LD(self.teff, self.logg, self.z, photband)

            except KeyError:
                return ldark.get_LD(self.teff, self.logg, self.z, photband)

        else:
            return ldark.get_LD(self.teff, self.logg, self.z, photband)

    def get_GD(self, photband):
        """
        Return gravity darkening in a given photometric band.
        """
        if isinstance(self.gd, float) or isinstance(self.gd, int):
            # Return same coefficient independetly of photband demanded
            return self.gd

        elif isinstance(self.gd, dict):
            try:
                gdi = self.gd[photband]
                if isinstance(gdi, float) or isinstance(gdi, int):
                    return gdi
                else:
                    return gdark.get_gd(self.teff)

            except KeyError:
                return gdark.get_gd(self.teff)

        else:
            return gdark.get_gd(self.teff)

    def get_albedo(self, photband):
        """
        Return albedo in a given photometric band.
        """
        if self.albedo is not None:
            if isinstance(self.albedo, float) or isinstance(self.albedo, int):
                # Return same coefficient independetly of photband demanded
                return self.albedo

            elif isinstance(self.albedo, dict):
                try:
                    albi = self.albedo[photband]
                    if isinstance(albi, float) or isinstance(albi, int):
                        return albi
                    else:
                        return albedo.get_albedo()

                except KeyError:
                    return albedo.get_albedo()

        else:
            return albedo.get_albedo()

    def get_beamingB(self, photband):
        """
        Compute B factor used in relativistic beaming effect.
        See Bloemen (MNRAS, 2011)
        """
        if isinstance(self.B, float) or isinstance(self.B, int):
            # Return same coefficient independetly of photband demanded
            return self.B

        elif isinstance(self.B, dict):
            try:
                Bi = self.B[photband]
                if isinstance(Bi, float) or isinstance(Bi, int):
                    return Bi
                else:
                    return phot.get_beaming(self.get_spectrum(), photband)

            except KeyError:
                return phot.get_beaming(self.get_spectrum(), photband)

        else:
            return phot.get_beaming(self.get_spectrum(), photband)

    def get_extinction(self):
        """
        Return EBmV for a given distance.
        """
        return EBmVfunc(self.dist)

    def get_BmV(self):
        """
        Compute B-V
        """
        B = self.get_mag('Johnson-B')
        V = self.get_mag('Johnson-V')
        return B - V

    def get_spots(self, t, photband='Kepler'):
        # Star
        """
        Compute spots model using macula
        https://www.cfa.harvard.edu/~dkipping/macula.html
        """
        if spotmodel == 'Macula' and self.Nspots > 0:

            star = n.array(
                [self.rotangle * n.pi / 180.0, self.period, self.kappa2,
                 self.kappa4, self.c1, self.c2, self.c3, self.c4, self.d1,
                 self.d2, self.d3, self.d4])
            inst = n.array([[1.0], [1.0]])
            spot = n.array([n.zeros(self.Nspots)] * 8)
            for i in range(self.Nspots):
                for j, spotparam in enumerate(
                        ['lambda0', 'phi0', 'alphamax', 'fspot', 'tmax',
                         'life', 'ingress', 'egress']):
                    spot[j][i] = self.__getattribute__(
                        'spot' + str(i) + spotparam)
                    if spotparam == 'lambda0' or spotparam == 'phi0' or \
                                    spotparam == 'alphamax':
                        spot[j][i] = spot[j][i] * n.pi / 180.0
                        # longitud, latitud and spot size in radians

            lcspots = macula(t, star, spot, inst)

            return lcspots

        else:
            return n.ones(len(t))


class Target(Star):
    def __init__(self, **kwargs):
        Star.__init__(self, **kwargs)

        self.test_critical(['teff', 'logg'], ['Effective temperature',
                                              'Surface gravity'])

        # Get parameters from tracks
        self.get_stellarparameters()

    def get_stellarparameters(self, output=False):
        """Get parameters from evolution models."""
        # Get Target Star parameters from evolution tracks
        mact, logL, logage = iso.get_stellarparams_target(self.z,
                                                          self.logg,
                                                          n.log10(self.teff))
        self.L = 10 ** logL
        self.logage = logage
        self.mact = mact

        # Compute radius
        self.R = n.sqrt(
            self.L * (5777.0 / self.teff) ** 4.0 / (1 - self.alphaS))

        if output:
            return mact, self.R, self.L, logage


class Blend(Star):
    def __init__(self, **kwargs):
        Star.__init__(self, **kwargs)

        self.test_critical(['minit', 'logage'], ['Effective temperature',
                                                 'Age'])

        # Get parameters from tracks
        self.get_stellarparameters()

    def get_stellarparameters(self, output=False):
        """
        Evaluate for input age, z, Minit
        OUTPUT: Mact, R, L, Teff, logg
        """

        # Get Blend Star parameters from evolution tracks
        logT, logg, logL, mact = iso.get_stellarparams(self.z, self.logage,
                                                       self.minit)

        self.teff = 10 ** logT
        self.logg = logg
        self.mact = mact
        self.L = 10 ** logL

        # Compute radius
        self.R = n.sqrt(
            self.L * (5777.0 / self.teff) ** 4.0 / (1 - self.alphaS))

        if output:
            return mact, self.R, self.L, self.teff, self.logg


class PlanetHost(Star):
    def __init__(self, **kwargs):
        Star.__init__(self, **kwargs)

        self.test_critical(['teff', 'dens'], ['Effective temperature',
                                              'Mean stellar density'])

        # Get parameters from tracks
        self.get_stellarparameters()

    def get_stellarparameters(self, output=False):
        """
        OUTPUT: Mact, R, L, logage, logg
        """
        mact, logL, logage = iso.get_stellarparams_target(self.z, self.dens,
                                                          log10(self.teff),
                                                          planethost=True)
        self.mact = mact
        self.logage = logage
        self.L = 10 ** logL

        # Compute radius
        self.R = n.sqrt(
            self.L * (5777.0 / self.teff) ** 4.0 / (1 - self.alphaS))

        # Compute logg from mact, and R
        self.logg = log10(
            G * 1e3 * (self.mact * Msun * 1e3) / (self.R * Rsun * 1e2) ** 2.0)

        if output:
            return mact, self.R, self.L, logage, self.logg


class WhiteDwarf(Star):
    def __init__(self, **kwargs):

        Star.__init__(self, **kwargs)

        self.test_critical(['teff', 'logg'], ['Effective temperature',
                                              'Surface gravity'])

        self.mact = iso.get_WD_mass(self.teff, self.logg)
        self.logage = iso.get_WD_logage(self.teff, self.logg)
        # self.R = kwargs.pop('R', n.sqrt((G*1e3*self.mact*Msun*1e3)/
        # (10**self.logg))/(Rsun*1e2))
        self.R = n.sqrt(
            (G * 1e3 * self.mact * Msun * 1e3) / (10 ** self.logg)) / (
                     Rsun * 1e2)

    def get_spectrum(self):
        """
        Get the spectrum of the white dwarf
        """
        # print('get_spectrum white dwarf')
        # print(self.teff,self.logg,self.ebmv,self.dist)
        try:
            return self.spectrum
        except AttributeError:
            pass

        teff = self.teff
        logg = self.logg

        try:
            spectrum = phot.get_interpolated_WD(teff, logg)
        except SpectrumInterpolError:
            spectrum = phot.get_nearest_WD(teff, logg)

        # Correct from distance, extinction, and radius
        spectrum = spectrum * n.pi * (self.R / (self.dist * pc2Rsun)) ** 2
        spectrum = spectrum * 10 ** (-0.4 * ext_law(ww) * self.ebmv)
        self.spectrum = spectrum

        return spectrum


class Triple(object):
    """
    Class for triple systems.

    This class will contain the following attributes and methods.

    Attributes
    ----------

    1) floats

    P : orbital period of the binary
    incl: orbital inclination (in degrees?)
    ecc: orbital eccentricity
    omega: argument of passage through periastron (in degrees?)
    T0: cental time of transit
    Tp: time of passage through periastron

    d: distance (should be the same as primary.d)
    age: the age of the system (same as primary.age)
    z: metallicity (same as primary.z)
    Av: extinction (same as primary.Av)
    ebmv: colour excess (same as primary.ebmv)
    v0 : systemic radial velocity

    2) objects

    primary: the primary star. A Star object.
    secondary: the secondary star. A Star object.


    Methods
    -------

    K1: the semi-amplitude of the primary star [in km/s??].
    K2: the semi-amplitude of the secondary star [in km/s??].
    RV1(t): the radial velocity of the primary star at time t [in km/s?]
    RV2(t): the radial velocity of the secondary star at time t [in km/s?]
    LC(t, band): the light curve of the binary on filter 'band'.


    """

    def __init__(self, orbital_parameters, Object1, Object2):

        self.dist = Object1.dist
        self.z = Object1.z
        self.logage = Object1.logage
        self.v0 = Object1.v0
        self.drift1 = Object1.drift1
        self.drift2 = Object1.drift2
        self.drift3 = Object1.drift3
        self.ebmv = Object1.ebmv

        self.object1 = Object1
        self.object2 = Object2

        self.orbital_parameters = orbital_parameters

        Object1._parent = self
        Object2._parent = self
        self._parent = None

    def get_parent(self):
        return self._parent

    def get_K(self, component='object1'):
        """
        Return the Radial Velocity semi-amplitude in km/s
        """
        if component == 'object1':
            obj1 = self.object1
            obj2 = self.object2
            ind = 0.
        elif component == 'object2':
            obj1 = self.object2
            obj2 = self.object1
            ind = 1.
        else:
            print('Component must be \"object1\" or \"object2\"')
            return

        orbital_parameters = self.orbital_parameters

        K = obj2.mact * Msun / (
            ((obj2.mact + obj1.mact) * Msun) ** (2. / 3.)) * (
                (2. * n.pi * G) ** (1. / 3.)) / (
                (orbital_parameters.P * 86400.) ** (1. / 3.)) * n.sin(
                    orbital_parameters.incl) / n.sqrt(
                    1. - orbital_parameters.ecc ** 2) / 1000.

        return K * (-1) ** ind

    def get_phase_transit(self, t):
        """
        Return the orbital phase, in the transit definition
        """
        phase_0 = 2. * n.pi * ((t - self.orbital_parameters.T0) %
                               self.orbital_parameters.P /
                               self.orbital_parameters.P)
        return phase_0

    def get_phase_periastron(self, t):
        """
        Return the orbital phase, in the periastron passage definition.
        """
        phase_p = 2. * n.pi * ((t - self.orbital_parameters.Tp) %
                               self.orbital_parameters.P /
                               self.orbital_parameters.P)
        return phase_p

    def get_RV(self, t, isphase=False, component='object1'):

        if self.orbital_parameters is None:
            if component == 'object1':
                return n.ones(len(t)) * self.object1.v0 + \
                       (self.drift1 * (t - TrefRV) +
                        self.drift2 * (t - TrefRV) ** 2 +
                        self.drift3 * (t - TrefRV) ** 3)
            elif component == 'object2':
                return n.ones(len(t)) * self.object2.v0 + \
                       (self.drift1 * (t - TrefRV) +
                        self.drift2 * (t - TrefRV) ** 2 +
                        self.drift3 * (t - TrefRV) ** 3)

        omega = self.orbital_parameters.omega
        ecc = self.orbital_parameters.ecc

        if isphase:
            # Convert transit phase to periastron phase
            Mc = self.orbital_parameters.get_E0() - ecc * n.sin(
                self.orbital_parameters.get_E0())
            M = 2. * pi * t - Mc
        else:
            M = self.get_phase_periastron(t)

        K = self.get_K(component=component)

        nu = tools.trueanomaly(M, ecc)

        rv = self.v0 + K * (n.cos(nu + omega) + ecc * n.cos(omega)) + (
             self.drift1 * (t - TrefRV) + self.drift2 * (t - TrefRV) ** 2 +
             self.drift3 * (t - TrefRV) ** 3)
        return rv


class Planet(object):
    """
    Class for planetary systems.

    This class will contain the following attributes and methods.

    Attributes
    ----------

    1) floats

    P : orbital period of the planet
    incl: orbital inclination (in radians)
    ecc: orbital eccentricity
    omega: argument of passage through periastron (in radians)
    T0: cental time of transit
    Tp: time of passage through periastron

    2) objects

    primary: the primary star. A Star object.

    Methods
    -------

    K1: the semi-amplitude of the primary star [in km/s??].
    RV1(t): the radial velocity of the primary star at time t [in km/s?]
    LC(t, band): the light curve of the binary on filter 'band'.

    """

    def __init__(self, orbital_parameters, Mp, Rp, f=0., **kwargs):

        self.Mp = Mp  # in Jupiter mass
        self.mact = Mp * Mjup2Msun  # in Solar masses

        self.Rp = Rp  # in Jupiter radii
        self.R = Rp * Rjup2Rsun  # in Solar radii

        self.orbital_parameters = orbital_parameters
        self.f = f

        self.albedo = kwargs.pop('albedo', None)
        self.gd = kwargs.pop('gd', 0.0)
        self.ua = kwargs.pop('ua', 0.0)
        self.ub = kwargs.pop('ub', 0.0)

        self._parent = None

    def get_parent(self):
        return self._parent

    def get_phase_transit(self, t):
        """
        Return the orbital phase, in the transit definition
        """
        orbital_parameters = self.orbital_parameters
        phase_0 = 2*n.pi * (t - orbital_parameters.T0)/orbital_parameters.P % 1
        return phase_0

    def get_phase_periastron(self, t):
        """
        Return the orbital phase, in the periastron passage definition.
        """
        orbital_parameters = self.orbital_parameters
        phase_p = 2*n.pi * (t - orbital_parameters.Tp)/orbital_parameters.P % 1
        return phase_p

    def get_LD(self, photband):
        return [self.ua, self.ub]

    def get_GD(self, photband):
        return self.gd

    def get_albedo(self, photband):
        """
        Return albedo in a given photometric band.
        """
        if isinstance(self.albedo, float) or isinstance(self.albedo, int):
            # Return same coefficient independently of photband demanded
            return self.albedo

        elif isinstance(self.albedo, dict):
            try:
                albi = self.albedo[photband]

                if isinstance(albi, float) or isinstance(albi, int):
                    return albi

                else:
                    raise TypeError('Invalid albedo for photband '
                                    '%s' % photband)

            except KeyError:
                raise KeyError('No albedo for photband %s' % photband)

        else:
            raise TypeError('Invalid albedo for photband %s' % photband)

    def get_LC(self, t, photband='Kepler'):
        """
        Return the normalized LC of the planet.
        NOT TO BE DONE ANYTIME SOON
        """
        return t * 0.0 + 1.0

    def get_spectrum(self):
        """
        Desactivated for the moment.
        When a library of planet spectra becomes available....
        """
        return ww * 0.0


class PlanSys(object):
    def __init__(self, star, *planets, **kwargs):

        # Include star
        self.star = star

        self.dist = star.dist
        self.z = star.z
        self.logage = star.logage
        self.v0 = star.v0
        self.drift1 = star.drift1
        self.drift2 = star.drift2
        self.drift3 = star.drift3
        self.ebmv = star.ebmv

        star._parent = self

        mact = self.star.mact

        # Include planets
        self.planets = []
        for planet in planets:
            if not (isinstance(planet, Planet) or
                    isinstance(planet, FitPlanet)):
                print('PlanetHost only takes Planet or FitPlanet instances '
                      'as planets.')
                continue

            orbit = planet.orbital_parameters

            if isinstance(planet, Planet):
                planet.q = planet.mact / star.mact
                mact = mact + planet.mact

            elif isinstance(planet, FitPlanet):
                if planet.q is None:
                    # If no q is given, then compute mass iteratively
                    orbit = planet.orbital_parameters

                    if planet.orbital_parameters.incl is not None:

                        Mp = tools.iterative_mass(planet.K1, orbit.P,
                                                  self.star.mact,
                                                  orbit.ecc,
                                                  orbit.incl * 180.0 / pi)

                    elif planet.orbital_parameters.b is not None:

                        Mp = tools.iterative_mass(planet.K1, orbit.P,
                                                  self.star.mact,
                                                  orbit.ecc, orbit.b,
                                                  omega=orbit.omega,
                                                  Rs=self.star.R,
                                                  use_b=True)

                    else:
                        raise ValueError

                    planet.q = Mp / self.star.mact
                    mact = mact + Mp

                else:
                    # If given, use to compute total mass of the system
                    mact = mact + planet.q * self.star.mact

            if not isinstance(star, PlanetHost):
                star.dens = star.mact / star.R ** 3.0

            # Impose stellar density
            # Compute a/R* for this planet based on stellar density,
            # and period.
            Ps = planet.orbital_parameters.P * 86400  # in seconds
            ar3 = (0.25 * G / pi**2 * (Msun / Rsun ** 3) * (1 + planet.q) *
                   self.star.dens * Ps ** 2)
            planet.ar = ar3 ** (1. / 3.)

            if (planet.orbital_parameters.bprime is not None and
                    planet.ar is not None):
                cosi = planet.orbital_parameters.bprime / planet.ar
                incl = n.arccos(cosi)
                planet.orbital_parameters.incl = incl

            elif (planet.orbital_parameters.b is not None and
                  planet.ar is not None):
                ecc = planet.orbital_parameters.ecc
                nu0 = pi / 2.0 - planet.orbital_parameters.omega
                bprime = planet.orbital_parameters.b / (1 - ecc ** 2) * (
                         1 + ecc * n.cos(nu0))
                planet.orbital_parameters.bprime = bprime
                cosi = bprime / planet.ar
                incl = n.arccos(cosi)
                planet.orbital_parameters.incl = incl

            self.planets.append(planet)
            planet._parent = self

        # Third light (secret) WARNING! May not work!
        self.f3 = kwargs.pop('f3', 0.0)

        self._parent = None

        self.mact = mact

    def get_parent(self):
        return self._parent

    def get_K(self, planet):
        """
        Return the Radial Velocity semi-amplitude in km/s
        """
        op = planet.orbital_parameters

        if isinstance(planet, Planet):
            K = planet.Mp * Mjup / (
                (planet.Mp * Mjup + self.star.mact * Msun)**(2./3.)) * (
                (2. * n.pi * G) ** (1. / 3.)) / (
                (op.P * 86400.) ** (1. / 3.)) * n.sin(op.incl) / n.sqrt(
                1. - op.ecc ** 2) / 1000.
            return K

        else:
            print('Something went wrong. One of the orbiting objects in this '
                  'planet system is not a planet!')

    def get_phase_transit(self, t, orbital_parameters):
        """
        Return the orbital phase, in the transit definition
        """
        phase_0 = 2.*n.pi*(t - orbital_parameters.T0)/orbital_parameters.P % 1
        return phase_0

    def get_phase_periastron(self, t, orbital_parameters):
        """
        Return the orbital phase, in the periastron passage definition.
        """
        phase_p = 2.*n.pi*(t - orbital_parameters.Tp)/orbital_parameters.P % 1
        return phase_p

    def get_true_lat(self, t, orbital_parameters, isphase=False):
        # P = orbital_parameters.P
        # Tp = orbital_parameters.Tp
        ecc = orbital_parameters.ecc
        omega = orbital_parameters.omega

        if isphase:
            # Convert transit phase to periastron phase
            Mc = orbital_parameters.get_E0() - ecc * n.sin(
                orbital_parameters.get_E0())
            M = 2. * pi * t - Mc
        else:
            M = self.get_phase_periastron(t, orbital_parameters)
        nu = tools.trueanomaly(M, ecc)
        return omega + nu

    def get_RV(self, t, isphase=False):
        """
        Return the Radial Velocity curve of Star in km/s
        """
        rv = n.zeros(len(t), float) + self.star.v0 + (
             self.star.drift1 * (t - TrefRV) + self.star.drift2 *
             (t - TrefRV)**2 + self.star.drift3 * (t - TrefRV)**3)

        for pl in self.planets:
            orbit = pl.orbital_parameters

            if not isinstance(pl, FitPlanet):
                pl.K1 = self.get_K(pl)

            if isphase:
                # Convert transit phase to periastron phase
                Mc = orbit.get_E0() - orbit.ecc * n.sin(orbit.get_E0())
                M = 2 * pi * t - Mc
            else:
                M = self.get_phase_periastron(t, orbit)

            nu = tools.trueanomaly(M, orbit.ecc)
            # Add RV of current planet
            rv += pl.K1 * (n.cos(nu + orbit.omega) +
                           orbit.ecc * n.cos(orbit.omega))
        return rv

    def get_sbr(self, planet, photband):

        if isinstance(planet, Planet):
            sp1 = self.star.get_spectrum()
            sp2 = planet.get_spectrum()

            f1 = phot.get_flux(sp1, photband)
            f2 = phot.get_flux(sp2, photband)
            ll1 = f1 / (f1 + f2)
            ll2 = 1. - ll1

            ldc1 = self.star.get_LD(photband)
            ldc2 = planet.get_LD(photband)

            return n.double((ll2 / ll1) * (
                   (self.star.R / (planet.Rp * Rjup2Rsun)) ** 2.0) * (
                            1. - ldc1[0] / 3. - ldc1[1] / 6.) / (
                            1. - ldc2[0] / 3. - ldc2[1] / 6.))

        elif isinstance(planet, FitPlanet):
            return planet.get_sbr(photband)

    def get_LD(self, planet, photband):
        if isinstance(planet, Planet):
            return planet.get_LD(photband)
        elif isinstance(planet, FitPlanet):
            return planet.get_LD(photband, 'secondary')

    def get_albedo(self, planet, photband):
        if isinstance(planet, Planet):
            return planet.get_albedo(photband)
        elif isinstance(planet, FitPlanet):
            return planet.get_albedo(photband, 'secondary')

    def get_GD(self, planet, photband):
        if isinstance(planet, Planet):
            return planet.get_GD(photband)
        elif isinstance(planet, FitPlanet):
            return planet.get_GD(photband, 'secondary')

    def get_LC(self, t, photband='Kepler', isphase=False, dt0=0.0):
        # PlanSys

        if len(self.planets) > 1 and isphase:
            raise Exception(
                'Can\'t work on phase if more than one planet is present.')

        # Get limbdarkening coefficients of central star
        ldc1 = self.star.get_LD(photband)
        gd1 = self.star.get_GD(photband)
        albedo1 = self.star.get_albedo(photband)

        LCS = []

        # Prepare to run EBOP

        # Fix quadratic limbdarkening and get coefficients
        ldtype = n.empty((2,), dtype='i')
        ldtype[0] = int(4)  # quad LD law type for star A
        ldtype[1] = int(4)  # quad LD law type for star B

        # START ITERATION OVER PLANETS IN THE SYSTEM
        for planet in self.planets:

            # Prepare input parameters array
            v = n.zeros(67, 'd')

            # Prepare phase or time array
            if not isphase:
                t0 = planet.orbital_parameters.T0 + dt0
                # Compute phase
                ph = (t - t0) / planet.orbital_parameters.P % 1.0
            else:
                ph = t
                ph = ph.astype('d')

            # Prepare output array
            nph = int(len(ph))
            y = n.zeros(ph.shape, 'd')

            # Useful orbital parameters
            incl = planet.orbital_parameters.incl
            ecc = planet.orbital_parameters.ecc
            omega = planet.orbital_parameters.omega  # in radians

            if isinstance(planet, Planet):  # IS CLASSICAL PLANET

                # Ratio of radii
                planet.kr = planet.R / self.star.R

            # Get wavelength dependent parameters
            ldc2 = self.get_LD(planet, photband)
            albedo2 = self.get_albedo(planet, photband)
            gd2 = self.get_GD(planet, photband)
            sbr = self.get_sbr(planet, photband)

            # Get surface brightness ratio
            # Get fractional fluxes
            kll = sbr * planet.kr**2 * (1 - ldc2[0] / 3.0 - ldc2[1] / 6.0) / (
                  1 - ldc1[0] / 3.0 - ldc1[1] / 6.0)
            ll1 = 1.0 / (1.0 + kll)
            ll2 = 1.0 - ll1

            r1 = 1 / planet.ar
            r2 = r1 * planet.kr

            v[18 - 1] = n.double(1.0)  # Integ. ring size (deg)
            v[2 - 1] = n.double(r1 + r2)  # Sum of radii (normalised to sma)
            v[3 - 1] = n.double(planet.kr)  # Ratio of the radii
            v[6 - 1] = n.double(incl * 180.0 / pi)  # Orbital inclination (deg)
            v[13 - 1] = n.double(planet.q)  # Mass ratio of system
            v[7 - 1] = n.double(
                ecc * n.cos(omega))  # e # e cos(omega) OR ecentricity
            v[8 - 1] = n.double(
                ecc * n.sin(omega))  # omega  #e sin(omega) OR omega
            v[9 - 1] = n.double(gd1)  # Gravity darkening (star A)
            v[10 - 1] = n.double(gd2)  # Grav darkening (star B)
            v[1 - 1] = n.double(sbr)  # Surface brightness ratio
            v[15 - 1] = n.double(0.0)  # Amount of third light
            v[4 - 1] = n.double(ldc1[0])  # LD star A (linear coeff)
            v[5 - 1] = n.double(ldc2[0])  # LD star B (linear coeff)
            v[21 - 1] = n.double(ldc1[1])  # LD star A (nonlin coeff)
            v[22 - 1] = n.double(ldc2[1])  # LD star B (nonlin coeff)
            # Reflection effect star A
            v[11 - 1] = n.double(0.5 * albedo1 * ll2 * r1 ** 2.0)
            # Reflection effect star B
            v[12 - 1] = n.double(0.5 * albedo2 * ll1 * r2 ** 2.0)
            v[16 - 1] = n.double(0.0)  # Phase shift of primary min
            v[17 - 1] = n.double(0.0)  # Light scale factor (mag)

            y = run_EBOP(v, ldtype, ph, nph, Nmax=Nmax, components=beaming)

            if beaming:
                # Include Beaming effect
                # y = [[LP],[LS],[ECL],[REFL]]

                # get componet radial velocity
                VR1 = (self.get_RV(t) - self.star.v0) * 1e3  # m/s

                B1 = self.star.get_beamingB(photband)

                # neglecting beaming effect in reflection
                y = y[0] * (1. - B1 * VR1 / c) + y[1] - y[2] + y[3]

            # Dilute EBOP output using third light
            y = y * (1 - self.f3) + self.f3

            LCS.append(y)

        return (1 - n.sum((1 - n.array(LCS)), axis=0)) * (
                self.star.get_spots(t, photband) * (1 - self.f3) + self.f3)
        # the flux of the planets are not taken into account


class FitBinary(object):
    def __init__(self, **kwargs):

        self._parent = None

        self.kr = kwargs.pop('kr', None)
        self.sumr = kwargs.pop('sumr', None)
        self.K1 = kwargs.pop('K1', None)
        self.q = kwargs.pop('q', None)  # Mass ratio M2/M1
        self.v0 = kwargs.pop('v0', None)
        self.drift1 = kwargs.pop('drift1', 0.0)
        self.drift2 = kwargs.pop('drift2', 0.0)
        self.drift3 = kwargs.pop('drift3', 0.0)
        self.vsini1 = kwargs.pop('vsini1', None)
        self.vsini2 = kwargs.pop('vsini2', None)
        self.zeta1 = kwargs.pop('zeta1', 0.0)
        self.zeta2 = kwargs.pop('zeta2', 0.0)

        # Set orbital parameters in two steps to avoid spurious
        # warning message
        self.orbital_parameters = kwargs.pop('orbital_parameters', None)
        if self.orbital_parameters is None:
            self.orbital_parameters = orbital_parameters(**kwargs)

        self.ar = kwargs.pop('ar', None)

        if self.ar is not None and self.kr is not None:
            self.sumr = (1 + self.kr) / self.ar

        # Modification to use impact parameter
        # self.b = kwargs.pop('b', None)
        # self.bprime = kwargs.pop('bprime', None)

        # Compute inclination in case a/Rs is an additional free parameter.
        if self.orbital_parameters.bprime is not None and self.ar is not None:
            cosi = self.orbital_parameters.bprime / self.ar
            incl = n.arccos(cosi)
            self.orbital_parameters.incl = incl

        elif self.orbital_parameters.b is not None and self.ar is not None:
            ecc = self.orbital_parameters.ecc
            nu0 = pi / 2.0 - self.orbital_parameters.omega
            bprime = self.orbital_parameters.b / (1 - ecc ** 2) * (
                1 + ecc * n.cos(nu0))
            self.orbital_parameters.bprime = bprime
            cosi = bprime / self.ar
            incl = n.arccos(cosi)
            self.orbital_parameters.incl = incl

        # Compute additional parameters for IsoBinary Class
        if isinstance(self, IsoBinary):
            # Compute semi-major axis of binary.
            self.sma = (((G * (self.mact * Msun) * (
                        self.orbital_parameters.P * 24.0 * 3600.0) ** 2.0) / (
                        4.0 * (n.pi ** 2))) ** (1.0 / 3.0)) / (Rsun)

            # Compute sum of radii
            self.sumr = (self.star1.R + self.star2.R) / self.sma

            # Compute semi-amplitude of RV for primary
            orbit = self.orbital_parameters

            self.K1 = ((2.0*n.pi*G / (orbit.P*86400))**(1. / 3.) *
                       (self.star2.mact * Msun * n.sin(orbit.incl) *
                       (self.mact * Msun)**(-2./3.)) *
                       (1 / n.sqrt(1. - orbit.ecc ** 2)) * 1e-3)

        # Third light (secret)
        self.f3 = kwargs.pop('f3', 0.0)

        # Wavelength dependent parameters
        self.sbr = kwargs.pop('sbr', None)  # Surface brightness ratio
        self.ua1 = kwargs.pop('ua1', None)
        self.ub1 = kwargs.pop('ub1', None)
        self.ua2 = kwargs.pop('ua2', None)
        self.ub2 = kwargs.pop('ub2', None)
        self.gd1 = kwargs.pop('gd1', None)
        self.gd2 = kwargs.pop('gd2', None)
        self.albedo1 = kwargs.pop('albedo1', None)
        self.albedo2 = kwargs.pop('albedo2', None)
        self.B1 = kwargs.pop('B1', None)  # For beaming effect
        self.B2 = kwargs.pop('B2', None)

        if spotmodel == 'Macula':
            # Parameters for spot modeling
            if not isinstance(self, IsoBinary):
                # For IsoBinary is handle inside stars objects
                for comp in ['comp1', 'comp2']:
                    for spotstarparams in ['rotangle', 'period', 'kappa2',
                                           'kappa4', 'c1', 'c2', 'c3', 'c4',
                                           'd1', 'd2', 'd3', 'd4']:
                        self.__setattr__(comp + spotstarparams,
                                         kwargs.pop(comp + spotstarparams,
                                                    None))

                    if self.__getattribute__(comp+'c1') is None or \
                            self.__getattribute__(comp+'c2') is None or \
                            self.__getattribute__(comp+'c3') is None or \
                            self.__getattribute__(comp+'c4') is None:
                        if self.__getattribute__('ua'+comp[-1]) is not None \
                         and self.__getattribute__('ub'+comp[-1]) \
                         is not None:
                            # Impose quadratic limb-darkening from star
                            self.__setattr__(comp+'c1', 0)
                            self.__setattr__(comp+'c2', self.__getattribute__(
                                'ua'+comp[-1]) + 2 * self.__getattribute__(
                                'ub'+comp[-1]))
                            self.__setattr__(comp+'c3', 0)
                            self.__setattr__(comp+'c4',
                                             -self.__getattribute__(
                                                 'ub' + comp[-1]))

                    Nspots = 0
                    for key in kwargs.keys():
                        if comp in key and 'life' in key:
                            Nspots += 1

                    self.__setattr__(comp + 'Nspots', Nspots)

                    if self.__getattribute__(comp + 'Nspots') != 0:
                        for i in range(self.__getattribute__(comp + 'Nspots')):
                            for spotparam in ['lambda0', 'phi0', 'alphamax',
                                              'fspot', 'tmax', 'life',
                                              'ingress', 'egress']:
                                self.__setattr__(
                                    comp + 'spot' + str(i) + spotparam,
                                    kwargs.pop(
                                        comp + 'spot' + str(i) + spotparam,
                                        None))

    def get_parent(self):
        return self._parent

    def get_phase_transit(self, t):
        """
        Return the orbital phase, in the transit definition
        """
        P = self.orbital_parameters.P
        T0 = self.orbital_parameters.T0
        phase_0 = ((t - T0) / P) % 1
        return 2. * n.pi * phase_0

    def get_phase_periastron(self, t):
        """
        Return the orbital phase, in the periastron passage definition.
        """

        t = n.atleast_1d(t)
        # FUTURE: compute for series of parameters t = n.atleast_2d(t)

        P = self.orbital_parameters.P
        Tp = self.orbital_parameters.Tp
        phase_p = (t - Tp)/P % 1
        return 2. * n.pi * phase_p

    def get_true_lat(self, t, isphase=False):
        # P = self.orbital_parameters.P
        # Tp = orbital_parameters.Tp
        ecc = self.orbital_parameters.ecc
        omega = self.orbital_parameters.omega

        if isphase:
            # Convert transit phase to periastron phase
            Mc = self.orbital_parameters.get_E0() - ecc * n.sin(
                self.orbital_parameters.get_E0())
            M = 2. * pi * t - Mc
        else:
            M = self.get_phase_periastron(t)
        nu = tools.trueanomaly(M, ecc)
        return omega + nu

    def get_RV(self, t, isphase=False, component='primary', istransit=True):

        t = n.atleast_1d(t)
        # FUTURE: perform computation for list of params
        # t = n.atleast_1d(t)[:, n.newaxis]

        # All others, converted to 2-D
        omega = self.orbital_parameters.omega
        ecc = self.orbital_parameters.ecc

        if isphase and istransit:
            # Convert transit phase to periastron phase
            Mc = self.orbital_parameters.get_E0() - ecc * n.sin(
                self.orbital_parameters.get_E0())
            M = 2. * pi * t - Mc
        elif isphase:
            M = 2. * pi * t
        else:
            M = self.get_phase_periastron(t)

        if component == 'primary':
            K = self.K1
        elif component == 'secondary':
            K = -self.K1 / self.q

        # Compute true anomaly
        nu = tools.trueanomaly(M, ecc)

        if not isphase:
            rv = self.v0 + K * (n.cos(nu + omega) + ecc * n.cos(omega)) + (
                self.drift1 * (t - TrefRV) + self.drift2 * (
                    t - TrefRV) ** 2 + self.drift3 * (t - TrefRV) ** 3)

        else:
            rv = self.v0 + K * (n.cos(nu) + ecc * n.cos(omega)) + (
                self.drift1 * (t - TrefRV) + self.drift2 * (
                    t - TrefRV) ** 2 + self.drift3 * (t - TrefRV) ** 3)

        return rv

    def get_LD(self, photband, component='primary'):
        """
        Return quadratic limb-darkening parameters in a given photometric band.
        """
        if component == 'primary':
            ua = self.ua1
            ub = self.ub1
        elif component == 'secondary':
            ua = self.ua2
            ub = self.ub2

        if (isinstance(ua, float) or isinstance(ua, int)) and (
             isinstance(ub, float) or isinstance(ub, int)):
            # Return same coefficients independetly of photband demanded
            return n.array([ua, ub])

        elif isinstance(ua, dict) and isinstance(ub, dict):
            try:
                uai, ubi = ua[photband], ub[photband]

                if (isinstance(uai, float) or isinstance(uai, int)) and (
                     isinstance(ubi, float) or isinstance(ubi, int)):
                    return n.array([uai, ubi])
                else:
                    raise TypeError(
                        'Invalid limbdarkening for photband %s' % photband)

            except KeyError:
                raise KeyError('No limbdarkening provided for photband '
                               '%s' % photband)

        else:
            raise TypeError('Invalid limbdarkening for photband %s' % photband)

    def get_GD(self, photband, component='primary'):
        """
        Return gravity darkening in a given photometric band.
        """
        if component == 'primary':
            gd = self.gd1
        elif component == 'secondary':
            gd = self.gd2

        if isinstance(gd, float) or isinstance(gd, int):
            # Return same coefficient independently of photband demanded
            return gd

        elif isinstance(gd, dict):
            try:
                gdi = gd[photband]

                if isinstance(gdi, float) or isinstance(gdi, int):
                    return gdi

                else:
                    raise TypeError(
                        'Invalid gravity darkening for photband %s' % photband)

            except KeyError:
                raise KeyError(
                    'No gravity darkening for photband %s' % photband)

        else:
            raise TypeError(
                'Invalid gravity darkening for photband %s' % photband)

    def get_albedo(self, photband, component='primary'):
        """
        Return albedo in a given photometric band.
        """
        if component == 'primary':
            albedo = self.albedo1
        elif component == 'secondary':
            albedo = self.albedo2

        if isinstance(albedo, float) or isinstance(albedo, int):
            # Return same coefficient independently of photband demanded
            return albedo

        elif isinstance(albedo, dict):
            try:
                albi = albedo[photband]

                if isinstance(albi, float) or isinstance(albi, int):
                    return albi

                else:
                    raise TypeError('Invalid albedo for photband '
                                    '%s' % photband)

            except KeyError:
                raise KeyError('No albedo for photband %s' % photband)

        else:
            raise TypeError('Invalid albedo for photband %s' % photband)

    def get_sbr(self, photband):
        """
        Return surface brightness ratio in a given photometric band.
        """
        if isinstance(self.sbr, float) or isinstance(self.sbr, int):
            # Return same coefficient independetly of photband demanded
            return self.sbr

        elif isinstance(self.sbr, dict):
            try:
                sbr = self.sbr[photband]

                if isinstance(sbr, float) or isinstance(sbr, int):
                    return sbr
                else:
                    raise TypeError('Invalid sbr for photband %s' % photband)

            except KeyError:
                raise KeyError('No sbr for photband %s' % photband)

        else:
            raise TypeError('Invalid sbr for photband %s' % photband)

    def get_beamingB(self, photband, component='primary'):
        """
        Return beaming coefficient in a given photometric band.
        """
        if component == 'primary':
            B = self.B1
        elif component == 'secondary':
            B = self.B2

        if isinstance(B, float) or isinstance(B, int):
            # Return same coefficient independently of photband demanded
            return B

        elif isinstance(B, dict):
            try:
                Bi = B[photband]

                if isinstance(Bi, float) or isinstance(Bi, int):
                    return Bi

                else:
                    raise TypeError('Invalid B for photband %s' % photband)

            except KeyError:
                raise KeyError('No B for photband %s' % photband)

        else:
            raise TypeError('Invalid B for photband %s' % photband)

    def get_spots(self, t, photband, component='primary'):
        # FitBinary
        """
        get spots of one of the binary components
        """
        if component == 'primary':
            comp = 'comp1'
        elif component == 'secondary':
            comp = 'comp2'
        else:
            raise NameError('Component not correctly specified.')

        if spotmodel == 'Macula' and \
           self.__getattribute__(comp + 'Nspots') != 0:

            star = n.array(
                [(self.__getattribute__(comp + 'rotangle')) * n.pi / 180.0,
                 self.__getattribute__(comp + 'period'),
                 self.__getattribute__(comp + 'kappa2'),
                 self.__getattribute__(comp + 'kappa4'),
                 self.__getattribute__(comp + 'c1'),
                 self.__getattribute__(comp + 'c2'),
                 self.__getattribute__(comp + 'c3'),
                 self.__getattribute__(comp + 'c4'),
                 self.__getattribute__(comp + 'd1'),
                 self.__getattribute__(comp + 'd2'),
                 self.__getattribute__(comp + 'd3'),
                 self.__getattribute__(comp + 'd4')])
            inst = n.array([[1.0], [1.0]])
            spot = n.array(
                [n.zeros(self.__getattribute__(comp + 'Nspots'))] * 8)
            for i in range(self.__getattribute__(comp + 'Nspots')):
                for j, spotparam in enumerate(
                        ['lambda0', 'phi0', 'alphamax', 'fspot', 'tmax',
                         'life', 'ingress', 'egress']):
                    spot[j][i] = self.__getattribute__(
                        comp + 'spot' + str(i) + spotparam)
                    if spotparam == 'lambda0' or spotparam == 'phi0' or \
                            spotparam == 'alphamax':
                        # longitude, latitude and spot size in radians
                        spot[j][i] = spot[j][i] * n.pi / 180.0

            lcspots = macula(t, star, spot, inst)
            return lcspots

        else:
            return n.ones(len(t))

    def get_LC(self, t, photband='Kepler', isphase=False, dt0=0.0):
        # FitBinary
        """
        Return the normalized LC of the planetary system.

        if dt0 is given, it is used to correct the transit ephemeris. This is
        used to measure individual transit times.
        """

        # Prepare phase or time array
        if not isphase:
            t0 = self.orbital_parameters.T0 + dt0
            # Compute phase
            ph = (t - t0)/self.orbital_parameters.P % 1.0
        else:
            ph = t

        ph = ph.astype('d')
        ###
        # PREPARE INPUT PARAMETERS FOR JKTEBOP
        ###

        # Prepare input parameters array
        v = n.zeros(67, 'd')

        # Fix quadratic limbdarkening and get coefficients
        ldtype = n.empty((2,), dtype='i')
        ldtype[0] = int(4)  # quad LD law type for star A
        ldtype[1] = int(4)  # quad LD law type for star B

        # Prepare output array
        nph = int(len(ph))
        y = n.zeros(ph.shape, 'd')

        # Get wavelength dependent parameters for primary and secondary star
        ldc1 = self.get_LD(photband, 'primary')
        gd1 = self.get_GD(photband, 'primary')
        albedo1 = self.get_albedo(photband, 'primary')
        ldc2 = self.get_LD(photband, 'secondary')
        gd2 = self.get_GD(photband, 'secondary')
        albedo2 = self.get_albedo(photband, 'secondary')
        sbr = self.get_sbr(photband)
        # Useful orbital parameters
        incl = self.orbital_parameters.incl
        ecc = self.orbital_parameters.ecc
        omega = self.orbital_parameters.omega  # in radians
        # Ps = self.orbital_parameters.P * 86400  # in seconds

        # Get fractional fluxes
        kll = sbr * self.kr ** 2 * (1 - ldc2[0] / 3.0 - ldc2[1] / 6.0) / (
              1 - ldc1[0] / 3.0 - ldc1[1] / 6.0)
        ll1 = 1.0 / (1.0 + kll)
        ll2 = 1.0 - ll1

        # Get radii
        r1 = self.sumr / (1 + self.kr)
        r2 = r1 * self.kr

        v[1 - 1] = n.double(sbr)  # Surface brightness ratio
        v[2 - 1] = n.double(self.sumr)  # Sum of the radii (normalised to sma)
        v[3 - 1] = n.double(self.kr)  # Ratio of the radii
        v[4 - 1] = n.double(ldc1[0])  # LD star A (linear coeff)
        v[5 - 1] = n.double(ldc2[0])  # LD star B (linear coeff)
        v[6 - 1] = n.double(incl * 180.0 / pi)  # Orbital inclination (deg)
        v[7 - 1] = n.double(ecc * n.cos(omega))  # e OR ecos(omega)
        v[8 - 1] = n.double(ecc * n.sin(omega))  # omega OR e sin(omega)
        v[9 - 1] = n.double(gd1)  # Gravity darkening (star A)
        v[10 - 1] = n.double(gd2)  # Grav darkening (star B)
        # Reflection effect star A
        v[11 - 1] = n.double(0.5 * albedo1 * ll2 * r1 ** 2.0)
        # Reflection effect star B
        v[12 - 1] = n.double(0.5 * albedo2 * ll1 * r2 ** 2.0)
        v[13 - 1] = n.double(self.q)  # Mass ratio of system
        v[15 - 1] = n.double(0.0)  # Amount of third light
        v[16 - 1] = n.double(0.0)  # Phase shift of primary min
        v[17 - 1] = n.double(0.0)  # Light scale factor (mag)
        v[18 - 1] = n.double(1.0)  # Integ. ring size (deg)
        v[21 - 1] = n.double(ldc1[1])  # LD star A (nonlin coeff)
        v[22 - 1] = n.double(ldc2[1])  # LD star B (nonlin coeff)

        # get light curve
        y = run_EBOP(v, ldtype, ph, nph, Nmax=Nmax, components=True)

        # include spots
        y1 = y[0] * self.get_spots(t, photband, component='primary')
        y2 = y[1] * self.get_spots(t, photband, component='secondary')

        if beaming:
            # Include Beaming effect
            # y = [[LP],[LS],[ECL],[REFL]]

            # get componet radial velocity in m/s
            VR1 = (self.get_RV(t, component='primary') - self.v0) * 1e3
            VR2 = (self.get_RV(t, component='secondary') - self.v0) * 1e3

            B1 = self.get_beamingB(photband, component='primary')
            B2 = self.get_beamingB(photband, component='secondary')

            # neglecting beaming effect in reflection
            y = y1 * (1. - B1 * VR1 / c) + y2 * (1. - B2 * VR2 / c) - y[2] + y[
                3]

        else:
            # No Beaming
            y = y1 + y2 - y[2] + y[3]

        # Dilute EBOP output using third light
        y = y * (1 - self.f3) + self.f3

        return y


class IsoBinary(FitBinary):
    """
    Star + Star
    """

    def __init__(self, orbital_parameters, Star1, Star2):

        self.star1 = Star1
        self.star2 = Star2
        self.mact = self.star1.mact + self.star2.mact

        # Compute a/R1 for FitObsBinary
        ar = (G / (4.0 * pi ** 2) * (Star1.mact + Star2.mact) * Msun *
              (orbital_parameters.P * 24 * 3600.0) ** 2) ** (1. / 3.) / (
              Star1.R * Rsun)

        # Get parameters from Star1 and Star2 to use with FitObsBinary
        params = {'kr': Star2.R / Star1.R,
                  'ar': ar,
                  'q': Star2.mact / Star1.mact,
                  'v0': Star1.v0,
                  'drift1': Star1.drift1,
                  'drift2': Star1.drift2,
                  'drift3': Star1.drift3,
                  'vsini1': Star1.vsini,
                  'vsini2': Star2.vsini,
                  'orbital_parameters': orbital_parameters,
                  'b': orbital_parameters.b,
                  'ua1': Star1.ua,
                  'ub1': Star1.ub,
                  'ua2': Star2.ua,
                  'ub2': Star2.ub,
                  'gd1': Star1.gd,
                  'gd2': Star2.gd,
                  'albedo1': Star1.albedo,
                  'albedo2': Star2.albedo,
                  'B1': Star1.B,
                  'B2': Star2.B
                  }

        FitBinary.__init__(self, **params)

        self.dist = Star1.dist
        self.z = Star1.z
        self.logage = Star1.logage
        self.v0 = Star1.v0
        self.drift1 = Star1.drift1
        self.drift2 = Star1.drift2
        self.drift3 = Star1.drift3
        self.ebmv = Star1.ebmv

        Star1._parent = self
        Star2._parent = self
        self._parent = None

    def get_LD(self, photband, component='primary'):
        if component == 'primary':
            return self.star1.get_LD(photband)
        elif component == 'secondary':
            return self.star2.get_LD(photband)

    def get_GD(self, photband, component='primary'):
        if component == 'primary':
            return self.star1.get_GD(photband)
        elif component == 'secondary':
            return self.star2.get_GD(photband)

    def get_albedo(self, photband, component='primary'):
        if component == 'primary':
            return self.star1.get_albedo(photband)
        elif component == 'secondary':
            return self.star2.get_albedo(photband)

    def get_beamingB(self, photband, component='primary'):
        if component == 'primary':
            return self.star1.get_beamingB(photband)
        elif component == 'secondary':
            return self.star2.get_beamingB(photband)

    def get_sbr(self, photband):

        sp1 = self.star1.get_spectrum()
        sp2 = self.star2.get_spectrum()

        f1 = phot.get_flux(sp1, photband)
        f2 = phot.get_flux(sp2, photband)
        ll1 = f1 / (f1 + f2)
        ll2 = 1. - ll1

        ldc1 = self.star1.get_LD(photband)
        ldc2 = self.star2.get_LD(photband)

        return n.double(ll2 / ll1 * (self.star1.R / self.star2.R) ** 2.0 *
                        (1. - ldc1[0] / 3. - ldc1[1] / 6. /
                        (1. - ldc2[0] / 3. - ldc2[1] / 6.)))

    def get_spots(self, t, photband, component='primary'):
        # IsoBinary
        if component == 'primary':
            return self.star1.get_spots(t, photband)
        elif component == 'secondary':
            return self.star2.get_spots(t, photband)

    def get_mag(self, photband):
        """
        Get combined magnitude of binary system
        """
        m1 = self.star1.get_mag(photband)
        m2 = self.star2.get_mag(photband)

        return m1 - 2.5 * log10(1 + 10 ** (-0.4 * (m2 - m1)))


class qBinary(IsoBinary):
    """
    Star + mass ratio
    """

    def __init__(self, orbital_parameters, primary, q, **kwargs):
        # Mass of the secondary (solar mass)
        M2 = q * primary.mact

        # Create secondary Blend object
        secparams = {'minit': M2,
                     'logage': primary.logage,
                     'z': primary.z,
                     'dist': primary.dist,
                     'v0': primary.v0,
                     'ebmv': primary.ebmv,
                     'vsini': kwargs.pop('vsini2', None),
                     'albedo': kwargs.pop('albedo2', None),
                     'ua': kwargs.pop('ua2', None),
                     'ub': kwargs.pop('ub2', None),
                     'gd': kwargs.pop('gd2', None),
                     'B': kwargs.pop('B', None),
                     'alphaS': kwargs.pop('alphaS2', 0.0)
                     }

        secondary = Blend(**secparams)

        # Initialize IsoBinary
        IsoBinary.__init__(self, orbital_parameters, primary, secondary)


class FitPlanet(FitBinary):
    def __init__(self, **kwargs):
        FitBinary.__init__(self, **kwargs)
        self._parent = None
        # parameter needed for Rossiter-McLaughlin effect
        self.BmV = kwargs.pop('BmV', None)

        # Fixed parameters for planet
        self.sbr = kwargs.pop('sbr', 0.0)
        self.gd2 = kwargs.pop('gd2', 0.0)
        self.B2 = kwargs.pop('B2', 0.0)
        self.ub2 = kwargs.pop('ub2', 0.0)
        self.ua2 = kwargs.pop('ua2', 0.0)
        self.albedo1 = kwargs.pop('albedo1', 0.4)

        # By default, set self.q to zero
        self.q = kwargs.pop('q', 0.0)


class orbital_parameters(object):
    """
    Sub class for orbital parameters
    ChangeLog :
    2012-03-28 : omega changed to radian
    """

    def __init__(self, **kwargs):

        self.P = kwargs.pop('P', None)

        # # Eccentricity and omega
        self.ecc = kwargs.pop('ecc', None)
        self.omega = kwargs.pop('omega', None)
        if self.omega is not None:
            self.omega = self.omega * pi / 180.0

        self.ecos = kwargs.pop('ecos', None)
        self.esin = kwargs.pop('esin', None)

        self.secos = kwargs.pop('secos', None)
        self.sesin = kwargs.pop('sesin', None)

        # Inclination
        self.incl = kwargs.pop('incl', None)
        if self.incl is not None:
            self.incl = self.incl * pi / 180.0

        self.b = kwargs.pop('b', None)
        self.bprime = kwargs.pop('bprime', None)

        # Times of passage
        self.T0 = kwargs.pop('T0', None)
        self.Tp = kwargs.pop('Tp', None)

        # In case an epoch and some angles at that epoch are given.
        self.epoch = kwargs.pop('epoch', None)
        self.M0 = kwargs.pop('M0', None)  # Mean anomaly at epoch
        if self.M0 is not None:
            self.M0 = self.M0 * pi / 180.0
        self.L0 = kwargs.pop('L0', None)  # Mean longitude at epoch
        if self.L0 is not None:
            self.L0 = self.L0 * pi / 180.0

        # Argument of ascending node; longitude of periapsis
        self.Omega = kwargs.pop('Omega', 0.0)
        if self.Omega is not None:
            self.Omega = self.Omega * pi / 180.0
        self.pomega = kwargs.pop('pomega', None)
        if self.pomega is not None:
            self.pomega = self.pomega * pi / 180.0

        # Obliquity
        self.spinorbit = kwargs.pop('spinorbit', None)
        if self.spinorbit is not None:
            self.spinorbit *= pi / 180.0

        if self.omega is not None and self.Omega is not None:
            self.pomega = self.omega + self.Omega

        # If omega is not given, check if longitude of ascending node and
        # longitude of periapsis are given, and compute omega from these.
        elif self.pomega is not None and self.Omega is not None:
            self.omega = self.pomega - self.Omega

        # Convert sqrt(ecos)(omega) and sqrt(e)sin(omega) to ecc and omega
        if self.secos is not None and self.sesin is not None:

            if self.ecc is not None or self.omega is not None:
                # Both ecos/esin and ecc/omega given; print warning
                print('WARNING! Langrangean orbital elements and eccentricity '
                      'and/or omega given. Will keep consistency with '
                      'Lagrangean elements.')

            self.ecc = self.secos ** 2 + self.sesin ** 2
            self.omega = n.arctan2(self.sesin, self.secos)

        # Convert ecos(omega) and esin(omega) to ecc and omega
        if self.ecos is not None and self.esin is not None:

            if self.ecc is not None or self.omega is not None:
                # Both ecos/esin and ecc/omega given; print warning
                print('WARNING! Langrangean orbital elements and eccentricity '
                      'and/or omega given. Will keep consistency with '
                      'Lagrangean elements.')

            self.ecc = n.sqrt(self.ecos ** 2 + self.esin ** 2)
            self.omega = n.arctan2(self.esin, self.ecos)

        # If ecc > 1, return error
        # if self.ecc >= 1:
        if n.any(self.ecc >= 1):
            raise ValueError('Eccentricity larger than 1.')

        # Compute T0 from Tp and viceversa; compute M0, L0 and epoch
        if self.Tp is not None:

            # Set epoch and mean anomaly at epoch
            self.epoch = self.Tp
            self.M0 = 0.0

            # Set mean longitude
            if self.Omega is not None:
                self.L0 = self.M0 + self.Omega + self.omega

            if self.T0 is not None:
                # Both Tp and T0 given; print warning
                print('WARNING! Both Tp and T0 given. Will keep consistency '
                      'with Tp.')
            try:
                self.T0 = self.get_T0()
            except AttributeError:
                print('Orbital parameters missing!')

        elif self.T0 is not None:
            try:
                self.Tp = self.get_Tp()
            except AttributeError:
                print('Orbital parameters missing!')

            # Set epoch and mean anomaly at epoch
            self.epoch = self.Tp
            self.M0 = 0.0

            # Set mean longitude
            if self.Omega is not None:
                self.L0 = self.M0 + self.Omega + self.omega

        # If nor T0 not Tp are given, then an epoch must be given
        elif self.epoch is not None and \
                (self.M0 is not None or self.L0 is not None):

            if self.L0 is not None and self.Omega is not None:

                if self.M0 is not None:
                    print('WARNING! Both mean anomaly at epoch and mean '
                          'longitude at epoch are given. Will keep '
                          'consistency with L0.')

                # Compute mean anomaly at eopch
                self.M0 = self.L0 - self.Omega - self.omega

            # Compute time of passage by periastron and
            # central time of transit
            self.Tp = self.epoch - self.M0 * self.P / (2 * pi)
            self.T0 = self.get_T0()

            if self.L0 is None and self.Omega is not None:
                # Compute mean longitude at epoch, if longitude of
                # ascending node is given
                self.L0 = self.M0 + self.Omega + self.omega
        else:
            print('Orbital parameters missing!')
        return

    def get_Tp(self):
        """
        Return the epoch of periastron from other orbital parameters
        """
        E_0 = self.get_E0()
        Tp = self.T0 - self.P / (2. * n.pi) * (E_0 - self.ecc * n.sin(E_0))
        return Tp

    def get_T0(self):
        """
        Return the transit epoch from other orbital parameters
        """
        E_0 = self.get_E0()
        T0 = self.Tp + self.P / (2. * n.pi) * (E_0 - self.ecc * n.sin(E_0))
        return T0

    def get_E0(self):
        return n.arctan2(n.sqrt(1. - self.ecc ** 2) * n.cos(self.omega),
                         n.sin(self.omega) + self.ecc)


class Drift(object):
    def __init__(self, **kwargs):
        self.rv0 = kwargs.pop('rv0', 0.0)
        self.lin = kwargs.pop('lin', 0.0)
        self.quad = kwargs.pop('quad', 0.0)
        self.cub = kwargs.pop('cub', 0.0)
        self.unitconstant = kwargs.pop('unitconstant', 1.0)

        # Define reference time. If absent from dict, use value in module.
        self.TrefRV = kwargs.pop('TrefRV', TrefRV)

    def get_RV(self, t):

        t = n.atleast_1d(t)
        # FUTURE: perform computation for list of params
        # t = n.atleast_1d(t)[:, n.newaxis]

        return self.rv0 + self.unitconstant * (
            self.lin * (t - self.TrefRV)/365.25 +
            self.quad * (t - self.TrefRV)**2 / 365.25**2 +
            self.cub * (t - self.TrefRV)**3 / 365.25**3)


__all__ = ['Star', 'Target', 'Blend', 'PlanetHost', 'WhiteDwarf', 'Triple',
           'Planet', 'PlanSys', 'FitBinary', 'IsoBinary', 'qBinary',
           'FitPlanet', 'Drift']
