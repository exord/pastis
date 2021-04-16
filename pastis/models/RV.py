"""
New PASTIS_RV dev

# autopep8 used on this file on 2020-07-22; not tested yet.
"""
from math import *
import numpy as n
from scipy import optimize

# Intra-package imports
from ..paths import libpath
from . import RVgaussianFitError
from .. import AstroClasses as ac
from ..velocimetry import *
from . import Pyarome

# New imports and definitions for Fast Gauss version
from ctypes import cdll, c_double, c_int, POINTER
import os

# lib = cdll.LoadLibrary(os.path.dirname(os.path.realpath(__file__)) +
# "/fitgauss_c.so")
lib = cdll.LoadLibrary(os.path.join(libpath, 'cpp', 'fitgauss_c.so'))
# lib.gaussian_res.restype    = c_void
lib.gaussian_res.argtypes = [POINTER(c_double), POINTER(c_double),
                             POINTER(c_double), c_int, POINTER(c_double)]
# lib.gaussian_res_J.restype  = c_void
lib.gaussian_res_J.argtypes = [POINTER(c_double), POINTER(c_double),
                               POINTER(c_double), c_int, POINTER(c_double)]
##


def CCF_interpolate(RV, CCF, params):
    # Bissector calculus using Chris-like method of squared interpolation
    # The calculus itself is based on Chris variables and code

    k, v0, fwhm, c = params

    sigma = fwhm / 2. / n.sqrt(2. * n.log(2.))
    norm_CCF = -c / k * (1. - CCF / c)
    nstep = 100
    margin = 5
    depth = n.arange(nstep - 2 * margin + 1, dtype=float) / \
        nstep + float(margin) / nstep

    # mean RV for each segment of the CCF
    MeanRV = [(RV[i] + RV[i + 1]) / 2. for i in range(len(CCF) - 1)]

    # derivatives for each segment of the CCF
    ExpV = [n.exp(-(v - v0)**2 / 2 / sigma**2) / sigma**2 for v in MeanRV]
    dCCFdRV = [-(v - v0) * expV for v, expV in zip(MeanRV, ExpV)]
    d2CCFdRV2 = [
        ((v - v0)**2 / sigma**2 - 1) * expV for v,
        expV in zip(
            MeanRV,
            ExpV)]
    d2RVdCCF2 = n.array([-d1 / d2**3 for d1, d2 in zip(d2CCFdRV2, dCCFdRV)])

    # not-null range (a ver como simplificar)
    iRange = [ii for ii in range(len(CCF) - 1) if (max(norm_CCF[ii],
                                                       norm_CCF[ii + 1]) >= depth[0]) & (min(norm_CCF[ii],
                                                                                             norm_CCF[ii + 1]) <= depth[-1])]

    # parameters ?? for each segment of the CCF
    p = n.zeros([len(CCF), 3], 'd')
    p[iRange, 2] = n.array(d2RVdCCF2[iRange]) / 2.
    p[iRange, 1] = n.array([(RV[i +
                                1] -
                             RV[i] -
                             p[i, 2] *
                             (norm_CCF[i +
                                       1] ** 2 -
                              norm_CCF[i]**2)) /
                            (norm_CCF[i +
                                      1] -
                             norm_CCF[i]) for i in iRange])
    p[iRange, 0] = n.array(
        [RV[i] - p[i, 1] * norm_CCF[i] - p[i, 2] * norm_CCF[i]**2 for i in iRange])

    # Indexes where "norm_CCF > dd"
    Indexes = n.array([[n.where(norm_CCF > dd)[0][0] - 1,
                        n.where(norm_CCF > dd)[0][-1]] for dd in depth])
    IndexBlue, IndexRed = [ind[0]
                           for ind in Indexes], [ind[-1] for ind in Indexes]

    # Bisector definition
    bis_b = n.array([p[i_b, 0] + p[i_b, 1] * dd + p[i_b, 2] *
                     dd**2 for i_b, dd in zip(IndexBlue, depth)])
    bis_r = n.array([p[i_r, 0] + p[i_r, 1] * dd + p[i_r, 2] *
                     dd**2 for i_r, dd in zip(IndexRed, depth)])

    return depth, bis_b, bis_r, v0


def weight_bouchy(lambdas, flux):
    # calculates the weight as a function of the pixel using the formula 8 of
    # Bouchy, Pepe & Queloz (2001)
    dA_over_dlambda = (flux[1:] - flux[:-1]) / (lambdas[1:] - lambdas[:-1])
    # add one value in the end
    dA_over_dlambda = n.append(dA_over_dlambda, 0.0)

    weight = lambdas * lambdas * dA_over_dlambda * dA_over_dlambda / flux

    return weight


def Vasy(RV, CCF, Gauss_params):
    """ Function that calculates the vasy asymetry indicator """
    center_G = Gauss_params[1]

    depth, lambda_b, lambda_r, v0 = CCF_interpolate(
        RV - center_G, CCF, Gauss_params)

    # calculate weight of the left wing and the right one
    weight_r = weight_bouchy(lambda_r, depth)
    weight_b = weight_bouchy(lambda_b, depth)

    average_weight = (weight_r + weight_b) / 2.0

    Vasy_res = n.sum((weight_r - weight_b) * average_weight) / \
        n.sum(average_weight)

    return Vasy_res


def asym_gaussian(p, x):
    """
    p : [c, rv0, sig, Asym]
    """
    Gb = 1. - p[0] / 100. * n.exp(-(x[n.where(x < p[1])[0]] - p[1])
                                  ** 2 / (2. * (p[2] * (1. - p[3]))**2))
    Gr = 1. - p[0] / 100. * n.exp(-(x[n.where(x >= p[1])[0]] - p[1])
                                  ** 2 / (2. * (p[2] * (1. + p[3]))**2))
    return n.concatenate((Gb, Gr))


def gaussian(p, x):
    """
    p : [c, rv0, sig]
    """
    return 1. - p[0] / 100. * n.exp(-(x - p[1])**2 / (2. * p[2]**2))


def res(p, x, y):
    return y - gaussian(p, x)


def asym_res(p, x, y):
    return y - asym_gaussian(p, x)


def gaussian_res(p, x, y):
    rez = n.empty(x.size)
    lib.gaussian_res(
        p.ctypes.data_as(
            POINTER(c_double)), x.ctypes.data_as(
            POINTER(c_double)), y.ctypes.data_as(
                POINTER(c_double)), len(x), rez.ctypes.data_as(
                    POINTER(c_double)))
    return rez


def gaussian_res_J(p, x, y):
    rez = n.empty((4, x.size))
    lib.gaussian_res_J(
        p.ctypes.data_as(
            POINTER(c_double)), x.ctypes.data_as(
            POINTER(c_double)), y.ctypes.data_as(
                POINTER(c_double)), len(x), rez.ctypes.data_as(
                    POINTER(c_double)))
    return rez


def fitgauss_old(x, y, p, mode='normal'):
    """
    Fit a gaussian function to x,y data, with starting parameters p.

    Model: p[0]*1/sqrt(2*pi*p[2])*exp(-(x-p[1])**2/p[2]**2)

    mode : 'normal', 'asym'
    """
    from scipy import optimize
    if mode == 'normal':
        # precision required of 10cm/s
        p1, wf = optimize.leastsq(res, p, args=(
            x, y), full_output=False, xtol=1e-5)
    elif mode == 'asym':
        p1, wf = optimize.leastsq(asym_res, p, args=(
            x, y), full_output=False, xtol=1e-5)  # precision required of 10cm/s
    if wf not in (1, 2, 3, 4):
        print('Warning: problem with RV gaussian fit')
        return [999., 999., 999.], 0
    else:
        return p1, 1


def fitgauss(x, y, p, mode='normal'):
    """
    Fit a gaussian function to x,y data, with starting parameters p.

    Model: p[0]*1/sqrt(2*pi*p[2])*exp(-(x-p[1])**2/p[2]**2)

    mode : 'normal', 'asym'
    """
    if mode == 'normal':
        p = n.hstack([[1.], p])
        p[1] = -p[1] / 100
        p[3] = p[3] * n.sqrt(2.0)
        xtol = 1e-5  # precision required of 10cm/s
        p1, wf = optimize.leastsq(
            gaussian_res, p, Dfun=gaussian_res_J, args=(
                x, y), col_deriv=True, full_output=False, xtol=xtol)
        if wf not in (1, 2, 3, 4):
            print('Warning: problem with RV gaussian fit')
            return [999., 999., 999.], 0
        else:
            p1[1] = -p1[1] * 100
            # p1[3] could be negative! we should make abs
            p1[3] = abs(p1[3]) / n.sqrt(2.0)
            return p1[1:], 1
    elif mode == 'asym':
        p1, wf = optimize.leastsq(asym_res, p, args=(
            x, y), full_output=False, xtol=1e-5)  # precision required of 10cm/s
    if wf not in (1, 2, 3, 4):
        print('Warning: problem with RV gaussian fit')
        return [999., 999., 999.], 0
    else:
        return p1, 1


def vspan(x, y, p, mode='TB'):
    """
    Return the VSPAN bisector as defined in Boisse et al. (2009??)

    mode : TB (top and bottom) or RB (Red and Blue)
    """
    ll1 = p[1] - 3. * p[2]
    ll2 = p[1] - 1. * p[2]
    ll3 = p[1] + 1. * p[2]
    ll4 = p[1] + 3. * p[2]
    ll5 = p[1]

    stepccf = x[1] - x[0]
    if ll1 > x[0] and ll4 < x[-1]:

        ind1 = int((ll1 - x[0]) / stepccf)
        ind2 = int((ll2 - x[0]) / stepccf)
        ind3 = int((ll3 - x[0]) / stepccf)
        ind4 = int((ll4 - x[0]) / stepccf)
        ind5 = int((ll5 - x[0]) / stepccf)

        if mode == 'TB':
            x1 = n.concatenate((x[:ind2], x[ind3 + 1:]))
            y1 = n.concatenate((y[:ind2], y[ind3 + 1:]))
            x2 = n.concatenate((x[:ind1], x[ind2:ind3 + 1], x[ind4 + 1:]))
            y2 = n.concatenate((y[:ind1], y[ind2:ind3 + 1], y[ind4 + 1:]))
        elif mode == 'RB':
            x1 = n.concatenate((x[:ind5], x[ind4:]))
            y1 = n.concatenate((y[:ind5], y[ind4:]))
            x2 = n.concatenate((x[:ind1], x[ind5:]))
            y2 = n.concatenate((y[:ind1], y[ind5:]))

        try:
            p6, conv_6 = fitgauss(x2, y2, p)
            p5, conv_5 = fitgauss(x1, y1, p)
            CCF_RV_high = p5[1]
            CCF_RV_low = p6[1]

            if conv_5 and conv_6:
                vspan = CCF_RV_high - CCF_RV_low
            else:
                vspan = 999.
        except BaseException:
            vspan = 999.
    else:
        vspan = 999.
    return vspan


def fit_BiGauss(x, y, p):
    """
    p : [c, rv0, sig]
    """
    return p[1] - fitgauss(x, y, n.hstack([p, [0.0]]), mode='asym')[0][1]


def fit_BIS(x, y, p):
    """
    Return the Bisector Inverse Slope (BIS) as defined in Queloz et al. (2001)
    """
    seuil = n.linspace(0., 1., 101)
    y -= min(y)
    y /= max(y)
    left = n.where(n.logical_and(x <= p[1], x > p[1] - 5. * p[2]))[0]
    right = n.where(n.logical_and(x >= p[1], x < p[1] + 5. * p[2]))[0]
    Vleft = n.interp(seuil, 1. - y[left], x[left])
    Vright = n.interp(seuil, y[right], x[right])
    Vtop = 0.
    for j in range(10, 40):
        Vtop += (Vright[-(j + 1)] + Vleft[j]) / 2.
    Vtop /= 30.
    Vbottom = 0.
    for j in range(60, 90):
        Vbottom += (Vright[-(j + 1)] + Vleft[j]) / 2.
    Vbottom /= 30.
    return Vtop - Vbottom


def get_FWHM(vsini, BmV, spectrograph, output='FWHM'):
    """
    Return the FWHM of a CCF assuming a given vsini.

    Calibrated for HARPS, SOPHIE HE and SOPHIE HR
    """
    if spectrograph == 'CORALIE':
        sigma0 = 6.603 - 6.357 * BmV + 5.533 * BmV**2 - 1.454 * BmV**3
        A = 1.9
        FWHM = n.sqrt((vsini**2. + A**2. * sigma0**2.) / A **
                      2.) * (2. * n.sqrt(2. * n.log(2.)))
    if spectrograph == 'HARPS':
        sigma0 = 8.625 - 20.037 * BmV + 23.388 * \
            BmV**2 - 10.364 * BmV**3 + 1.273 * BmV**4
        A = 1.95
        FWHM = n.sqrt((vsini**2. + A**2. * sigma0**2.) / A **
                      2.) * (2. * n.sqrt(2. * n.log(2.)))
    if spectrograph == 'SOPHIE HE':
        sigma0 = 10.52 - 22.56 * BmV + 22.37 * BmV**2. - 6.95 * BmV**3.
        A = 1.64
        FWHM = n.sqrt((vsini**2. + A**2. * sigma0**2.) / A **
                      2.) * (2. * n.sqrt(2. * n.log(2.)))
    if spectrograph == 'SOPHIE HR':
        sigma0 = 9.90 - 22.56 * BmV + 22.37 * BmV**2. - 6.95 * BmV**3.
        A = 1.95
        FWHM = n.sqrt((vsini**2. + A**2. * sigma0**2.) / A **
                      2.) * (2. * n.sqrt(2. * n.log(2.)))

    if output == 'FWHM':
        return FWHM
    elif output == 'sigma0':
        return sigma0
    elif output == 'sigma':
        return FWHM / (2. * n.sqrt(2. * n.log(2.)))
    else:
        return


def get_contrast(FWHM, BmV, z, spectrograph, mask, output='ctrs'):
    """
    Return the contrast of a CCF assuming a given FWHM, (B-V) and metallicity z.

    Calibrated for SOPHIE and HARPS
    """
    if spectrograph == 'HARPS':

        if mask == 'K5':
            c = [-0.4232454896, 1.688452363, 0.2941150665, -
                 0.7595982552, -0.05445498228, -0.1592344046]
        if mask == 'G2':
            c = [-0.1295745075, 1.534696698, 0.2792906761, -
                 0.7796853781, -0.05642914772, -0.1540470123]
        logW = c[0] + c[1] * BmV + c[2] * z + c[3] * \
            BmV**2 + c[4] * z**2 + c[5] * BmV * z

    if spectrograph.find('SOPHIE') > -1:

        if mask == 'K5':
            logW = (z + 2.2606 * BmV - 0.2615) / 3.9553  # Boisse et al. 2010

        if mask == 'G2':
            logW = (z + 1.4992 * BmV + 0.9440) / 3.8807  # Boisse et al. 2010

    if spectrograph.find('CORALIE') > -1:

        if mask == 'K0':
            logW = (z - 2.573 + 8.142 * BmV - 5.583 * BmV**2) / \
                4.587  # Santos et al. 2002

    if output == 'ctrs':
        return 10**logW / n.sqrt(2. * n.pi) / \
            (FWHM / (2. * n.log(2.))) * 100.  # in %
    elif output == 'W':
        return 10**logW
    else:
        return


def get_Wrot(ctrs, Wccf, rv, vsini, epsilon, FWHM_instru):
    """
    Not working, do not use !
    """
    from scipy import integrate
    return Wccf - integrate.trapz(1. - profile_rot(rv,
                                                   0., vsini, epsilon, ctrs[0], FWHM_instru), rv)


def get_contrast_rot(rv, vsini, BmV, z, spectrograph, mask, epsilon):
    """
    Not working, do not use !
    """
    Wccf = get_contrast(1., BmV, z, spectrograph, mask, output='W')
    FWHM_instru = get_FWHM(0., BmV, spectrograph, output='sigma0')
    CRot, toto = optimize.leastsq(
        get_Wrot, [1.], args=(
            Wccf, rv, vsini, epsilon, FWHM_instru))
    return CRot


def profile_rot(x, x0, vsini, epsilon, ctrs, FWHM_instru=9.9):
    """
    Not working, do not use !
    """
    sig = 2. * n.sqrt(2. * n.log(2.))
    y1 = 1. / n.sqrt(2. * n.pi * (FWHM_instru / sig)**2) * \
        n.exp(-(x - x0)**2 / (2. * (FWHM_instru / sig)**2))
    # y1=n.exp(-(x-x0)**2/(2.*(FWHM_instru/sig)**2))*ctrs
    y2 = rot_profile(x, vsini, epsilon)
    # y2/=max(y2)
    Conv2 = n.convolve(y1[1:], y2[:-1], mode=1)
    Conv2 = n.array(list(Conv2) + [0.0])
    Conv2 /= max(Conv2)
    CCF = (1. - Conv2 * ctrs)
    # CCF -= min(CCF) - (1. - ctrs)
    # CCF /= max(CCF)
    return CCF


def rot_profile(v, vsini, epsilon):

    c1 = 2. * (1. - epsilon)
    c2 = 0.5 * n.pi * epsilon

    G = n.zeros(len(v), float)
    for i in range(len(v)):
        if abs(v[i]) <= vsini:
            Dv = (v[i] / vsini)**2.
            a = 1. - Dv
            b = n.sqrt(a)
            G[i] = (c1 * b) + (c2 * a)
            G[i] /= (n.pi * vsini * (1. - epsilon / 3.))
        else:
            G[i] = 0

    return G


def make_CCF(rv, rv0, FWHM, contrast, flux):
    """
    Make a synthetic gaussian CCF for a given radial velocity, FWHM, contrast and flux.
    """
    sigma = FWHM / (2. * n.sqrt(2. * n.log(2)))
    return gaussian([contrast, rv0, sigma], rv) * flux


def CCF_prop(t_rv, spectro, mask, *args):
    # count number of stars
    nb_star = 0
    for obj in args:
        if isinstance(obj, ac.Star):
            nb_star += 1
        if isinstance(obj, ac.IsoBinary):
            nb_star += 2
        if isinstance(obj, ac.Triple):
            for comp in (obj.object1, obj.object2):
                if isinstance(comp, ac.Star):
                    nb_star += 1
                elif isinstance(comp, ac.IsoBinary):
                    nb_star += 2
                elif isinstance(comp, ac.PlanSys):
                    nb_star += 1
        if isinstance(obj, ac.PlanSys):
            nb_star += 1
    FWHM = n.zeros(nb_star, float)
    contrast = n.zeros(nb_star, float)
    rv0 = n.zeros([nb_star, len(t_rv)], float)
    flux = n.zeros(nb_star, float)
    # define CCF proprieties of each stars
    count = 0
    for obj in args:
        if isinstance(obj, ac.Star):
            obj.BmV = obj.get_BmV()
            FWHM[count] = get_FWHM(obj.vsini, obj.BmV, spectro)
            contrast[count] = get_contrast(
                FWHM[count], obj.BmV, obj.z, spectro, mask)
            flux[count] = obj.get_flux('Johnson-V')
            rv0[count] = obj.v0 * n.ones(len(t_rv), float)
            count += 1
        if isinstance(obj, ac.IsoBinary):
            obj.star1.BmV = obj.star1.get_BmV()
            FWHM[count] = get_FWHM(obj.star1.vsini, obj.star1.BmV, spectro)
            contrast[count] = get_contrast(
                FWHM[count], obj.star1.BmV, obj.star1.z, spectro, mask)
            flux[count] = obj.star1.get_flux('Johnson-V')
            rv0[count] = obj.get_RV(t_rv, component='primary')
            count += 1
            obj.star2.BmV = obj.star2.get_BmV()
            FWHM[count] = get_FWHM(obj.star2.vsini, obj.star2.BmV, spectro)
            contrast[count] = get_contrast(
                FWHM[count], obj.star2.BmV, obj.star2.z, spectro, mask)
            flux[count] = obj.star2.get_flux('Johnson-V')
            rv0[count] = obj.get_RV(t_rv, component='secondary')
            count += 1

        if isinstance(obj, ac.Triple):
            for component, labelcomp in zip(
                    (obj.object1, obj.object2), ('object1', 'object2')):
                if isinstance(component, ac.Star):
                    component.BmV = component.get_BmV()
                    FWHM[count] = get_FWHM(
                        component.vsini, component.BmV, spectro)
                    contrast[count] = get_contrast(
                        FWHM[count], component.BmV, component.z, spectro, mask)
                    flux[count] = component.get_flux('Johnson-V')
                    rv0[count] = obj.get_RV(t_rv, component=labelcomp)
                    count += 1
                elif isinstance(component, ac.IsoBinary):
                    component.star1.BmV = component.star1.get_BmV()
                    FWHM[count] = get_FWHM(
                        component.star1.vsini, component.star1.BmV, spectro)
                    contrast[count] = get_contrast(
                        FWHM[count], component.star1.BmV, component.star1.z, spectro, mask)
                    flux[count] = component.star1.get_flux('Johnson-V')
                    rv0[count] = obj.get_RV(t_rv,
                                            component=labelcomp) + component.get_RV(t_rv,
                                                                                    component='primary') - component.v0
                    count += 1
                    component.star2.BmV = component.star2.get_BmV()
                    FWHM[count] = get_FWHM(
                        component.star2.vsini, component.star2.BmV, spectro)
                    contrast[count] = get_contrast(
                        FWHM[count], component.star2.BmV, component.star2.z, spectro, mask)
                    flux[count] = component.star2.get_flux('Johnson-V')
                    rv0[count] = obj.get_RV(t_rv,
                                            component=labelcomp) + component.get_RV(t_rv,
                                                                                    component='secondary') - component.v0
                    count += 1
                elif isinstance(component, ac.PlanSys):
                    component.star.BmV = component.star.get_BmV()
                    FWHM[count] = get_FWHM(
                        component.star.vsini, component.star.BmV, spectro)
                    contrast[count] = get_contrast(
                        FWHM[count], component.star.BmV, component.star.z, spectro, mask)
                    flux[count] = component.star.get_flux('Johnson-V')
                    rv0[count] = obj.get_RV(
                        t_rv, component=labelcomp) + component.get_RV(t_rv) - component.v0
                    count += 1
        if isinstance(obj, ac.PlanSys):
            obj.star.BmV = obj.star.get_BmV()
            FWHM[count] = get_FWHM(obj.star.vsini, obj.star.BmV, spectro)
            contrast[count] = get_contrast(
                FWHM[count], obj.star.BmV, obj.star.z, spectro, mask)
            flux[count] = obj.star.get_flux('Johnson-V')
            rv0[count] = obj.get_RV(t_rv)
            count += 1

    return nb_star, FWHM, contrast, flux, rv0


def PASTIS_RV(t_rv, RVdatadict, *args):
    """
    return models for RV, bisector, FWHM and contrast

    :param t_rv:
    :param RVdatadict:
    :param args:
    :return:
    """
    spectro = RVdatadict['spectro']
    mask = RVdatadict['mask']

    output_dict = {}

    observables = ['RV', 'CTRS', 'FWHM', 'BIS', 'Vspan', 'Wspan', 'BiGauss',
                   'Vasy']

    # print(RVdatadict)
    # If any of the observables is in RVdatadict, initialise arrays
    if any([obs in RVdatadict for obs in observables]):
        rv_simu = n.zeros(len(t_rv), float)
        fwhm_simu = n.zeros(len(t_rv), float)
        contrast_simu = n.zeros(len(t_rv), float)

    if 'BIS' in RVdatadict:
        bis_simu = n.zeros(len(t_rv), float)
    if 'Vspan' in RVdatadict:
        vspan_simu = n.zeros(len(t_rv), float)
    if 'Wspan' in RVdatadict:
        wspan_simu = n.zeros(len(t_rv), float)
    if 'BiGauss' in RVdatadict:
        bigauss_simu = n.zeros(len(t_rv), float)
    if 'Vasy' in RVdatadict:
        vasy_simu = n.zeros(len(t_rv), float)

    # condFitBin = map(isinstance, args, [ac.FitBinary]*len(args))
    condFitBin = [isinstance(x, ac.FitBinary) for x in args]
    # condIsoBin = map(isinstance, args, [ac.IsoBinary]*len(args))
    condIsoBin = [isinstance(x, ac.IsoBinary) for x in args]
    # condFitPlanet = map(isinstance, args, [ac.FitPlanet]*len(args))
    condFitPlanet = [isinstance(x, ac.FitPlanet) for x in args]
    # cond_drift = map(isinstance, args, [ac.Drift]*len(args))
    cond_drift = [isinstance(x, ac.Drift) for x in args]

    # If sole object is a Planetary System
    if len(args) == 1 and isinstance(args[0], ac.PlanSys):

        # Get stellar B-V, FWHM (and sigma0), and contrast
        args[0].star.BmV = args[0].star.get_BmV()
        FWHM = get_FWHM(args[0].star.vsini, args[0].star.BmV, spectro)
        sigma0 = get_FWHM(args[0].star.vsini, args[0].star.BmV, spectro,
                          output='sigma0')
        contrast = get_contrast(FWHM, args[0].star.BmV, args[0].star.z,
                                spectro, mask)
        CCF_width = FWHM / (2. * n.sqrt(2. * n.log(2.)))

        if 'RV' in RVdatadict:
            output_dict['RV'] = args[0].get_RV(
                t_rv) - RVdatadict['RV']['offset']
            for p in args[0].planets:
                if p.orbital_parameters.spinorbit is not None:
                    if isinstance(p, ac.Planet):
                        output_dict['RV'] += Pyarome.arome(
                            p.get_true_lat(t_rv),
                            p.ar *
                            (
                                1. -
                                p.orbital_parameters.ecc**2) /
                            (
                                1. +
                                p.orbital_parameters.ecc *
                                n.sin(
                                    p.orbital_parameters.omega)),
                            p.orbital_parameters.incl,
                            p.orbital_parameters.spinorbit,
                            [
                                args[0].star.ua,
                                args[0].star.ub],
                            beta0=sigma0,
                            Vsini=args[0].star.vsini,
                            sigma0=CCF_width,
                            zeta=args[0].star.zeta,
                            Rp=p.Rp /
                            p.star.R,
                            units='degree')
                    elif isinstance(p, ac.FitPlanet):
                        output_dict['RV'] += Pyarome.arome(
                            p.get_true_lat(t_rv),
                            p.ar *
                            (
                                1. -
                                p.orbital_parameters.ecc**2) /
                            (
                                1. +
                                p.orbital_parameters.ecc *
                                n.sin(
                                    p.orbital_parameters.omega)),
                            p.orbital_parameters.incl,
                            p.orbital_parameters.spinorbit,
                            [
                                args[0].star.ua,
                                args[0].star.ub],
                            beta0=sigma0,
                            Vsini=args[0].star.vsini,
                            sigma0=CCF_width,
                            zeta=args[0].star.zeta,
                            Rp=p.kr,
                            units='radian')
        if 'CTRS' in RVdatadict:
            output_dict['CTRS'] = (n.ones(len(t_rv), float) * contrast -
                                   RVdatadict['CTRS']['offset'])
        if 'FWHM' in RVdatadict:
            output_dict['FWHM'] = (n.ones(len(t_rv), float) * FWHM -
                                   RVdatadict['FWHM']['offset'])

        # Set all diagnosis to zero + offset except for RV.
        for obs in observables[3:]:
            if obs in RVdatadict:
                output_dict[obs] = (n.zeros_like(t_rv, dtype=float) -
                                    RVdatadict[obs]['offset'])

    elif all(condFitBin or cond_drift) and \
            not any(condIsoBin) and not any(condFitPlanet):
        # If all objects are FitBinary and none is an IsoBinary, and none is a FitPlanet
        # just return sum of RVs (in this way qBinaries will go through
        # the entire procedure).

        for obj in args:
            rv_simu += obj.get_RV(t_rv)

        if 'RV' in RVdatadict:
            output_dict['RV'] = rv_simu - RVdatadict['RV']['offset']

        # Set all diagnosis to zero + offset except for RV.
        for obs in observables[1:]:
            if obs in RVdatadict:
                output_dict[obs] = (n.zeros_like(t_rv, dtype=float) -
                                    RVdatadict[obs]['offset'])

    # TODO: merge with previos condition (yes, we can!)
    elif all(n.array(condFitPlanet) | n.array(cond_drift)):
        for obj in args:
            rv_simu += obj.get_RV(t_rv)

            # Add RM-effect, only if FitPlanet instance
            if (isinstance(obj, ac.FitPlanet) and
                    obj.orbital_parameters.spinorbit is not None):
                sigma0 = get_FWHM(obj.vsini1, obj.BmV, spectro,
                                  output='sigma0')
                CCF_width = get_FWHM(obj.vsini1, obj.BmV, spectro,
                                     output='sigma')
                rv_simu += Pyarome.arome(obj.get_true_lat(t_rv),
                                         obj.ar * (1. - obj.orbital_parameters.ecc**2) / (1. + obj.orbital_parameters.ecc * n.sin(obj.orbital_parameters.omega)),
                                         obj.orbital_parameters.incl,
                                         obj.orbital_parameters.spinorbit,
                                         [obj.ua1,
                                          obj.ub1],
                                         beta0=sigma0,
                                         Vsini=obj.vsini1,
                                         sigma0=CCF_width,
                                         zeta=obj.zeta1,
                                         Rp=obj.kr,
                                         units='radians')

        if 'RV' in RVdatadict:
            output_dict['RV'] = rv_simu - RVdatadict['RV']['offset']
        if 'CTRS' in RVdatadict:
            output_dict['CTRS'] = (n.ones(len(t_rv), float) * contrast -
                                   RVdatadict['CTRS']['offset'])
        if 'FWHM' in RVdatadict:
            output_dict['FWHM'] = (n.ones(len(t_rv), float) * FWHM -
                                   RVdatadict['FWHM']['offset'])

        # Set all diagnosis to zero + offset except for RV, CTRS, and FWHM.
        for obs in observables[3:]:
            if obs in RVdatadict:
                output_dict[obs] = (n.zeros_like(t_rv, dtype=float) -
                                    RVdatadict[obs]['offset'])

    else:
        nb_star, FWHM, contrast, flux, rv0 = CCF_prop(
            t_rv, spectro, mask, *args)

        # CCF parameters
        cond_RV_CCF = n.logical_and(
            RV_CCF_ALL > n.min(rv0) - 5. * max(FWHM) / 2.3548,
            RV_CCF_ALL < n.max(rv0) + 5. * max(FWHM) / 2.3548)
        RV_CCF = RV_CCF_ALL[n.where(cond_RV_CCF)[0]]
        CCF = n.zeros([nb_star, len(RV_CCF)], float)

        # Make fixed CCFs
        count = 0
        CCF_sum = n.zeros(len(RV_CCF), float)
        while isinstance(args[count], Star):
            CCF[count] = make_CCF(RV_CCF, rv0[count, 0],
                                  FWHM[count], contrast[count], flux[count])
            CCF_sum += CCF[count]
            count += 1
        CCF_star = CCF_sum.copy()

        for i in range(len(t_rv)):
            CCF_sum = CCF_star.copy()
            # Make non-fixed CCFs
            for o in range(count, nb_star):
                CCF[o] = make_CCF(RV_CCF, rv0[o, i], FWHM[o],
                                  contrast[o], flux[o])
                CCF_sum += CCF[o]

            normalized_CCF = n.array(CCF_sum / max(CCF_sum))

            # Fit the blended CCF

            if n.any([oo in RVdatadict for oo in observables]):
                p0 = [contrast[n.argmax(flux *
                                        contrast /
                                        100.)], rv0[n.argmax(flux *
                                                             contrast /
                                                             100.), i], FWHM[n.argmax(flux *
                                                                                      contrast /
                                                                                      100.)] /
                      2.3548]
                flag = 0
                try:
                    CCF_res, flag = fitgauss(RV_CCF, normalized_CCF, p0)
                except ValueError:
                    raise RVgaussianFitError(
                        "Problem with CCF fitting in PASTIS_RV.py")

                if flag == 0:
                    raise RVgaussianFitError(
                        "Problem with CCF fitting in PASTIS_RV.py")
                else:
                    rv_simu[i] = CCF_res[1]
                    fwhm_simu[i] = CCF_res[2] * 2.3548
                    contrast_simu[i] = CCF_res[0]

            if 'BIS' in RVdatadict:
                bis_simu[i] = fit_BIS(RV_CCF, normalized_CCF, CCF_res)
            if 'Vspan' in RVdatadict:
                vspan_simu[i] = vspan(
                    RV_CCF, normalized_CCF, CCF_res, mode='TB')
            if 'Wspan' in RVdatadict:
                wspan_simu[i] = vspan(
                    RV_CCF, normalized_CCF, CCF_res, mode='RB')
            if 'BiGauss' in RVdatadict:
                bigauss_simu[i] = fit_BiGauss(RV_CCF, normalized_CCF, CCF_res)
            if 'Vasy' in RVdatadict:
                vasy_simu[i] = weight_interp_depth(RV_CCF - CCF_res[1],
                                                   normalized_CCF,
                                                   [CCF_res[0] / 100., 0.,
                                                   CCF_res[2] * 2.3548])

        if 'RV' in RVdatadict:
            output_dict['RV'] = rv_simu - RVdatadict['RV']['offset']
        if 'CTRS' in RVdatadict:
            output_dict['CTRS'] = contrast_simu - RVdatadict['CTRS']['offset']
        if 'FWHM' in RVdatadict:
            output_dict['FWHM'] = fwhm_simu - RVdatadict['FWHM']['offset']
        if 'BIS' in RVdatadict:
            output_dict['BIS'] = bis_simu - RVdatadict['BIS']['offset']
        if 'Vspan' in RVdatadict:
            output_dict['Vspan'] = vspan_simu - RVdatadict['Vspan']['offset']
        if 'Wspan' in RVdatadict:
            output_dict['Wspan'] = wspan_simu - RVdatadict['Wspan']['offset']
        if 'BiGauss' in RVdatadict:
            output_dict['BiGauss'] = bigauss_simu - \
                RVdatadict['BiGauss']['offset']
        if 'Vasy' in RVdatadict:
            output_dict['Vasy'] = vasy_simu - RVdatadict['Vasy']['offset']

    return output_dict
