#!/usr/bin/env python

'''
05-04-2012: free extinction included in the star call

'''
import string
from .photometry import Allpbands
from .AstroClasses import *


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

    ## If no orbital parameters are
    ## specified, construct them from input dictionary.
    if not 'orbital_parameters' in dd.keys():
        ddop = {}
        for key in dd.keys():
            ddop[key] = dd[key][0]
        #dd['orbital_parameters'] = orbital_parameters(**ddop)
        orbit = orbital_parameters(**ddop)
    
    if 'FitBinary' in objname or 'FitPlanet' in objname:
        ddbin = {}
        for key in inputdict[objname].keys():  
            ddbin[key] = inputdict[objname][key][0]
        ddbin['orbital_parameters'] = orbit

        if 'FitBinary' in objname:
            binary = FitBinary(**ddbin)
    
        elif 'FitPlanet' in objname:
            binary = FitPlanet(**ddbin)

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

            binary = qBinary(primary = star1, **ddbin)

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

            binary = IsoBinary(orbit, star1, star2)

    return binary


def make_plansys(objname, dd, imposeobj = None):

    ## Make central star
    dd1 = {}
    
    starname = inputdict[objname]['star1'] ## CHECK NAME OF STAR IN QBIN
    for key in inputdict[starname].keys():  
        dd1[key] = inputdict[starname][key][0]
        
    ## If impositions are to be made (e.g. if plansys belongs to a triple system)
    if imposeobj != None:
        dd1['logage'] = imposeobj.logage
        dd1['z'] = imposeobj.z
        dd1['v0'] = imposeobj.v0
        dd1['dist'] = imposeobj.dist
        dd1['ebmv'] = imposeobj.ebmv

    star = make_star(starname, dd1)

    try:
        ## Remove star from list; inside try to use as standalone
        inputdict.pop(starname)
    except NameError:
        pass

    
    ## Construct planet objects
    planets = []
    for planet in dd.keys():

        if planet == 'star1' : continue

        ddp = {}
        for key in inputdict[dd[planet]].keys():
            ddp[key] = inputdict[dd[planet]][key][0]

        if not 'orbital_parameters' in ddp.keys():
            ddp['orbital_parameters'] = orbital_parameters(**ddp)
            
        ### BE CAREFUL; orbital_parameters DOES NOT REMOVE ELEMENTS FROM 
        ### DICTIONARY. AS A CONSECUENCE, PLANET EXPLODES.
        
        ## Construct Classical Planets
        if 'FitPlanet' in dd[planet]:
            planets.append(FitPlanet(**ddp))
            inputdict.pop(dd[planet])
        elif 'Planet' in dd[planet]:
            planets.append(Planet(**ddp))
            inputdict.pop(dd[planet])

    ## Construct PlanSys
    plansys = PlanSys(star, *planets)

    return plansys
    

def make_star(objname, dd, imposeobj = None):

    if imposeobj != None:
        dd['logage'] = imposeobj.logage
        dd['z'] = imposeobj.z
        dd['v0'] = imposeobj.v0
        dd['dist'] = imposeobj.dist
        dd['ebmv'] = imposeobj.ebmv

    if 'Target' in objname:
        star = Target(**dd)

    elif 'PlanetHost' in objname:
        star = PlanetHost(**dd)

    elif 'Blend' in objname:
        star = Blend(**dd)

    elif 'Star' in objname:
        star = Star(**dd)

    elif 'WhiteDwarf' in objname:
        star = WhiteDwarf(**dd)

    return star


def ObjectBuilder(dictionary) :
    """
    Construct objects from scenario dictionary.
    """
    import numpy as n
    objects = {}

    global inputdict
    inputdict = dictionary.copy()
    
    ###
    ### Get names of all objects to be constructed.
    ###
    toconstruct = inputdict.keys()

    ## Build all objects
    while len(inputdict.keys()) > 0:

        ## Start construction of Triples
        indTriple = n.array(map(string.find, inputdict.keys(),
                                ['Triple']*len(inputdict.keys())
                                ))

        triples = n.array(inputdict.keys())[indTriple != -1]
        
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
                dd['orbital_parameters'] = orbital_parameters(**ddop)

            ## Build Triple
            objects[obj] = Triple(dd['orbital_parameters'], obj1, obj2)
            

            ## Remove Triple object from list to construct.
            inputdict.pop(obj)

        ###
        ## Start construction of planetary systems
        ###
        plansystems = []
        for plansys in TypeList['plansys']:
            indPlanSys = n.array(map(string.find, inputdict.keys(),
                                     [plansys]*len(inputdict.keys())
                                     )
                                 )

            plansystems.extend(n.array(inputdict.keys())[indPlanSys != -1])

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
            indBinary = n.array(map(string.find, inputdict.keys(),
                                    [binary]*len(inputdict.keys())
                                    )
                                )

            binaries.extend(n.array(inputdict.keys())[indBinary != -1])

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
            indStars = n.array(map(string.find, inputdict.keys(),
                                   [star]*len(inputdict.keys())
                                   )
                               )

            stars.extend(n.array(inputdict.keys())[indStars != -1])
            
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
        ind_drifts = n.array(map(string.find, inputdict.keys(),
                                 ['Drift']*len(inputdict.keys())))

        drifts = n.array(inputdict.keys())[ind_drifts != -1]

        for obj in drifts:

            dd = {}
            for key in inputdict[obj].keys():
                dd[key] = inputdict[obj][key][0]

            objects[obj] = Drift(**dd)

        ####
        # Remove data-related objects
        ####
        pbands = []
        
        for pband in Allpbands:
            indPbands = n.array(map(string.find, inputdict.keys(),
                                    [pband]*len(inputdict.keys())
                                    )
                                )

            pbands.extend(n.array(inputdict.keys())[indPbands != -1])

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
            ind = n.array(map(string.find, inputdict.keys(),
                              [obj]*len(inputdict.keys())
                              )
                          )

        for remain in n.array(inputdict.keys())[ind == -1]:
            inputdict.pop(remain)

    return objects.values()
