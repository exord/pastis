#!/usr/bin/env python
"""
Module with functions to build objects from dictionary of parameters.

The main function, ObjectBuilder, recursively explores the input dictionary,
constructs the individual instances of the classes it requires, and then
puts them together in systems as required.

05-04-2012: free extinction included in the star call
"""

from . import TypeList
from .photometry import Allpbands
from . import AstroClasses as ac


def make_triple_component(objname, configdict, imposeobj=None):

    occurrence = 0

    ## Test if Object is a binary
    for binary in TypeList['binary']:
        occurrence = occurrence + objname.count(binary)

    if occurrence == 1:
        #print('Found component of triple: %s, of type binary'%objname)
        component = make_binary(objname, inputdict[objname], imposeobj)

    else:
        ## Test if Object is a planetary system
        for plansys in TypeList['plansys']:
            occurrence = occurrence + objname.count(plansys)

        if occurrence == 1:
            #print('Found component of triple: %s, of type plansys'%objname)
            component = make_plansys(objname, inputdict[objname], imposeobj)

        else:
            ## Test if Object is a star
            for star in TypeList['star']:
                occurrence = occurrence + objname.count(star)

            if occurrence == 1:
                #print('Found component of triple: %s, of type star'%objname)
                ddstar = {}
                for key in inputdict[objname].keys():
                    ddstar[key] = inputdict[objname][key][0]

                component = make_star(objname, ddstar, imposeobj)

            else:
                raise Exception('Object %s does not belong to a class that can be associated in a Triple.'%objname)

    """
    try:
        ## Remove component from list; inside try to use as standalone
        inputdict.pop(objname)
    except NameError:
        pass
    """
    return component


def make_binary(objname, dd, imposeobj = None):

    # If no orbital parameters are given, construct them from input 
    # dictionary.
    if not 'orbital_parameters' in dd.keys():
        ddop = {}
        for key in dd.keys():
            if not isinstance(dd[key], dict):
                ddop[key] = dd[key][0]
        #dd['orbital_parameters'] = orbital_parameters(**ddop)
        orbit = ac.orbital_parameters(**ddop)

    if 'FitBinary' in objname or 'FitPlanet' in objname:
        ddbin = {}
        for key in inputdict[objname].keys():
            if not isinstance(inputdict[objname][key], dict):
                ddbin[key] = inputdict[objname][key][0]
            else:
                ddbin[key] = dict((photband, 
                                   inputdict[objname][key][photband][0]) for
                                   photband in inputdict[objname][key])
        ddbin['orbital_parameters'] = orbit

        if 'FitBinary' in objname:
            binary = ac.FitBinary(**ddbin)

        elif 'FitPlanet' in objname:
            binary = ac.FitPlanet(**ddbin)

    elif 'IsoBinary' in objname or 'qBinary' in objname:
        ## Construct primary star
        dd1 = {}

        starname = inputdict[objname]['star1'] ## CHECK NAME OF STAR IN QBIN
        for key in inputdict[starname].keys():
            dd1[key] = inputdict[starname][key][0]

        ## If impositions are to be made (e.g. if binary belongs to a triple system)
        if imposeobj != None:
            dd1['logage'] = imposeobj.logage
            dd1['z'] = imposeobj.z
            dd1['v0'] = imposeobj.v0
            dd1['dist'] = imposeobj.dist
            dd1['ebmv'] = imposeobj.ebmv

        star1 = make_star(starname, dd1)

        try:
            ## Remove star from list; inside try to use as standalone
            inputdict.pop(starname)
        except NameError:
            pass

        if 'qBinary' in objname:
            ddbin = {}
            for key in inputdict[objname].keys():
                if key == 'star1': continue
                ddbin[key] = inputdict[objname][key][0]
            ddbin['orbital_parameters'] = orbit

            binary = ac.qBinary(primary=star1, **ddbin)

        elif 'IsoBinary' in objname:
            ## Construct secondary imposing primary age, z, dist, v0, and E(B-V)
            secondaryname = inputdict[objname]['star2']

            dd2 = {}
            for key in inputdict[secondaryname].keys():
                dd2[key] = inputdict[secondaryname][key][0]

            dd2['logage'] = star1.logage
            dd2['z'] = star1.z
            dd2['dist'] = star1.dist
            dd2['v0'] = star1.v0
            dd2['ebmv'] = star1.ebmv

            star2 = make_star(secondaryname, dd2)

            try:
                ## Remove star from list; inside try to use as standalone
                inputdict.pop(secondaryname)
            except NameError:
                pass

            binary = ac.IsoBinary(orbit, star1, star2)

    return binary


def make_plansys(objname, dd, imposeobj = None):

    # Make central star
    dd1 = {}

    # Get name of central star in planetary system
    starname = inputdict[objname]['star1'] ## CHECK NAME OF STAR IN QBIN
    
    for key in inputdict[starname].keys():
        if not isinstance(inputdict[starname][key], dict):
            dd1[key] = inputdict[starname][key][0]
        else:
            dd1[key] = dict((photband, 
                             inputdict[starname][key][photband][0]) for
                            photband in inputdict[starname][key])
        
    # Impose conditions as necessary
    # (e.g. if plansys belongs to a triple system)
    if imposeobj != None:
        dd1['logage'] = imposeobj.logage
        dd1['z'] = imposeobj.z
        dd1['v0'] = imposeobj.v0
        dd1['dist'] = imposeobj.dist
        dd1['ebmv'] = imposeobj.ebmv

    # Build central star instance
    star = make_star(starname, dd1)

    try:
        # Remove star from list; inside try to use as standalone
        inputdict.pop(starname)
    except NameError:
        pass

    # Construct planet instances
    planets = []
    for planet in dd.keys():

        if planet == 'star1':
            continue

        ddp = {}
        
        # Get planet name plansys dictionary
        planetname = dd[planet]
                
        for key in inputdict[planetname].keys():

            if not isinstance(inputdict[planetname][key], dict):
                ddp[key] = inputdict[planetname][key][0]
            else:
                ddp[key] = dict((photband, 
                                 inputdict[planetname][key][photband][0]) for
                                 photband in inputdict[planetname][key])

        if not 'orbital_parameters' in ddp.keys():
            ddp['orbital_parameters'] = ac.orbital_parameters(**ddp)

        ### BE CAREFUL; orbital_parameters DOES NOT REMOVE ELEMENTS FROM
        ### DICTIONARY. AS A CONSECUENCE, PLANET EXPLODES.

        ## Construct Classical Planets
        if 'FitPlanet' in dd[planet]:
            planets.append(ac.FitPlanet(**ddp))
            inputdict.pop(dd[planet])
        elif 'Planet' in dd[planet]:
            planets.append(ac.Planet(**ddp))
            inputdict.pop(dd[planet])

    ## Construct PlanSys
    plansys = ac.PlanSys(star, *planets)

    return plansys


def make_star(objname, dd, imposeobj = None):

    if imposeobj != None:
        dd['logage'] = imposeobj.logage
        dd['z'] = imposeobj.z
        dd['v0'] = imposeobj.v0
        dd['dist'] = imposeobj.dist
        dd['ebmv'] = imposeobj.ebmv

    if 'Target' in objname:
        star = ac.Target(**dd)

    elif 'PlanetHost' in objname:
        star = ac.PlanetHost(**dd)

    elif 'Blend' in objname:
        star = ac.Blend(**dd)

    elif 'Star' in objname:
        star = ac.Star(**dd)

    elif 'WhiteDwarf' in objname:
        star = ac.WhiteDwarf(**dd)

    return star


def ObjectBuilder(dictionary) :
    """Construct objects from scenario dictionary."""
    import numpy as n
    objects = {}

    global inputdict
    inputdict = dictionary.copy()

    ###
    ### Get names of all objects to be constructed.
    ###
    toconstruct = inputdict.keys()

    # Build all objects
    while len(inputdict.keys()) > 0:
        
        ## Start construction of Triples
        indTriple = n.array([a.find('Triple') for a in inputdict])

        triples = n.array(list(inputdict.keys()))[indTriple != -1]

        for obj in triples:

            ## Get config dict of triple
            dd = dictionary[obj]

            ## Construct first object of triple
            obj1 = make_triple_component(dd['object1'],
                                         dictionary[dd['object1']])

            ## Remove component from list to construct.
            inputdict.pop(dd['object1'])

            ## Construct second object of triple
            obj2 = make_triple_component(dd['object2'],
                                         dictionary[dd['object2']],
                                         imposeobj = obj1)

            ## Remove component from list to construct.
            inputdict.pop(dd['object2'])

            ## Make orbital parameters
            if not 'orbital_parameters' in dd.keys():
                ddop = {}
                for key in dd.keys():
                    ddop[key] = dd[key][0]
                dd['orbital_parameters'] = ac.orbital_parameters(**ddop)

            ## Build Triple
            objects[obj] = ac.Triple(dd['orbital_parameters'], obj1, obj2)


            ## Remove Triple object from list to construct.
            inputdict.pop(obj)

        ###
        ## Start construction of planetary systems
        ###
        plansystems = []
        for plansys in TypeList['plansys']:
            indPlanSys = n.array([a.find(plansys) for a in inputdict])

            plansystems.extend(n.array(list(inputdict.keys())
                              )[indPlanSys != -1])

        for obj in plansystems:

            ## Construct
            objects[obj] = make_plansys(obj, inputdict[obj])

            ## Remove Planetary System from list
            inputdict.pop(obj)

        ###
        ## Start construction of binaries
        ###
        binaries = []
        for binary in TypeList['binary']:
            indBinary = n.array([a.find(binary) for a in inputdict])

            binaries.extend(n.array(list(inputdict.keys())
                                    )[indBinary != -1])

        for obj in binaries:

            ## Construct
            objects[obj] = make_binary(obj, inputdict[obj])

            ## Remove created binaries from list
            inputdict.pop(obj)

        ###
        ## Start construction of stars
        ###
        stars = []

        for star in TypeList['star']:
            
            indStars = n.array([a.find(star) for a in inputdict])

            stars.extend(n.array(list(inputdict.keys()))[indStars != -1])

        for obj in stars:
            dds = {}
            for key in inputdict[obj].keys():
                dds[key] = inputdict[obj][key][0]

            objects[obj] = make_star(obj, dds)

            ## Remove created stars from list
            inputdict.pop(obj)

        ####
        # Construct Drifts
        ####
        ind_drifts = n.array([a.find('Drift') for a in inputdict])
        
        drifts = n.array(list(inputdict.keys()))[ind_drifts != -1]

        for obj in drifts:

            dd = {}
            for key in inputdict[obj].keys():
                dd[key] = inputdict[obj][key][0]

            objects[obj] = ac.Drift(**dd)

        ####
        # Remove data-related objects
        ####
        pbands = []

        for pband in Allpbands:
            indPbands = n.array([a.find(pband) for a in inputdict])

            pbands.extend(n.array(list(inputdict.keys()))[indPbands != -1])

        for pband in pbands:
            ## Remove photometric objects from list
            inputdict.pop(pband)

        ####
        # Remove all remaining keys, if they're not recognized
        ####
        ### BE CAREFUL! THIS IS DANGEROUS!
        ### IT MAY LEAD TO ERRORS
        allobjs = []
        for objtype in TypeList.keys():
            allobjs.extend(TypeList[objtype])

        for obj in allobjs:
            ind = n.array([a.find(obj) for a in inputdict])

        for remain in n.array(list(inputdict.keys()))[ind == -1]:
            inputdict.pop(remain)

    return list(objects.values())
