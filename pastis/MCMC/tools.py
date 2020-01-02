# Functions to transform input dictionary to list of parameters
# and to update this dictionary to feed it back to the Object Builder

import numpy as n

# Intra-package imports
from . import Objects_MCMC as objMCMC


def state_constructor(input_dict):
    """
    Construct the state of a Markov Chain starting from a input dictionary
    with a given structure (see PASTIS documentation?).
    """

    theta = []
    #global labeldict
    labeldict = {}

    # Iteration over all objects
    for objkey in input_dict.keys():

        # Iteration over all parameters of a given object
        for parkey in input_dict[objkey]:
            
            family = [objkey, parkey]

            if parkey == 'object':
                continue

            parlist = input_dict[objkey][parkey]

            # Special treatment for limb darkening
            if isinstance(parlist, dict):
                
                # Iterate over all photbands
                for photband in parlist:
                    
                    # Set previous level as parent
                    familytree = family.copy()
                    familytree.append(photband)
                    
                    # Define name
                    parname = '_'.join(familytree)
                                    
                    par = create_parameter(parname, familytree, 
                                           parlist[photband])

                    theta.append(par)
                    labeldict[par.label] = par
                    
            # If parameter has None value, continue
            elif parlist[0] is None:
                continue

            elif isinstance(parlist, list):
                if len(parlist) < 2:
                    continue
                ##
                # Construct parameter instance with information on dictionary
                parname = '_'.join(family)
                #parname = objkey+'_'+parkey
                par = create_parameter(parname, family, parlist)
                theta.append(par)
                labeldict[par.label] = par

    return n.array(theta), labeldict


def create_parameter(name, family, parlist):
    # Search jump size proposed by user
    try:
        if 'ebmv' not in name:
            jumpsize = parlist[8]
        else:
            raise IndexError
    except IndexError:
        # If not present, get prior information to do a reasonable
        # guess on the jump size
        if parlist[2] in ('Uniform', 'Jeffreys', 'Sine'):
            priorsize = 0.683*(parlist[4] - parlist[3])

        elif parlist[2] in ('Normal', 'LogNormal', 'TruncatedUNormal'):
            priorsize = parlist[4]*2.0

        elif parlist[2] == 'Binormal':
            priorsize = 2.0*max(parlist[4], parlist[6])

        elif parlist[2] == 'AssymetricNormal':
            priorsize = 2.0*max(parlist[4], parlist[5])

        elif parlist[2] == 'PowerLaw':
            priorsize = 0.683*(parlist[5] - parlist[4])

        elif parlist[0] != 0:
            priorsize = 2.0*abs(parlist[0]*0.1)

        else:
            priorsize = 1.0

        jumpsize = priorsize*0.5

        return objMCMC.Parameter(parlist[0], None, jump=parlist[1], label=name,
                                 family=family, proposal_scale=jumpsize)
    
    
    
""""
def state_deconstructor(state, input_dict):
    output_dict = {}
    for key in input_dict.keys():
        output_dict[key] = input_dict[key].copy()

    # Iterate over all parameters
    for par in state:

        objkey = par.label.split('_')[0]
        parkey = par.label.split('_')[1]

        if not (objkey in output_dict):
            print('Warning! Dictionary does not have key \"{}\"'.format(objkey))
            continue

        if not (parkey in output_dict[objkey]):
            raise KeyError('Error! Dictionary of object {0} does not have key '
                           '\"{1}\"'.format(objkey, parkey))

        output_dict[objkey][parkey] = []
        output_dict[objkey][parkey].append(par.get_value())
        for i in range(1, len(input_dict[objkey][parkey])):
            output_dict[objkey][parkey].append(input_dict[objkey][parkey][i])

    return output_dict
"""

"""
# WARNING! CHECK WHICH VERSION OF state_deconstructor IS BETTER IN TERMS
# OF OVERWRITING THE INPUT STATE
def state_deconstructor(state, output_dict):
    #
    # Iterate over all parameters
    for par in state:

        objkey = par.label.split('_')[0]
        parkey = par.label.split('_')[1]

        if not (objkey in output_dict):
            print('Warning! Dictionary does not have key \"{}\"'.format(objkey))
            continue

        if not (parkey in output_dict[objkey]):
            raise KeyError('Error! Dictionary of object {0} does not have key '
                           '\"{1}\"'.format(objkey, parkey))

        output_dict[objkey][parkey][0] = par.get_value()


    return output_dict
"""

def state_deconstructor(state, output_dict):
    #
    # Iterate over all parameters
    for par in state:

        objkey = par.family[0]
        parkey = par.family[1]

        if not (objkey in output_dict):
            print('Warning! Dictionary does not have key \"{}\"'.format(objkey))
            continue

        if not (parkey in output_dict[objkey]):
            raise KeyError('Error! Dictionary of object {0} does not have key '
                           '\"{1}\"'.format(objkey, parkey))
            
        # Get object dict for this parameter
        adict = output_dict[par.family[0]]
        
        # Descend on the dictionary until we reach the parameter level
        for i in range(1, len(par.family)):
            adict = adict[par.family[i]]
                
        # Set new value for parameter.
        adict[0] = par.get_value()

    return output_dict



def get_jitter(data, instrument, paramdict, observable=None):
    """
    Compute jitter given an element of a datadict (data).
    :return:
    """

    # # Construct jitter key (to match new keys in RV diagnostics).
    if observable is None:
        jitterkey = instrument + '_jitter'
    else:
        jitterkey = instrument + '_' + observable + 'jitter'

    if 'jittermodel' not in data or data['jittermodel'] == 'constant':
        try:
            return paramdict[jitterkey].get_value()
        except KeyError:
            return 0.0

    elif data['jittermodel'] == 'linear_rhk':
        minjitterkey = jitterkey.replace('jitter', 'minjitter')
        alphajitterkey = jitterkey.replace('jitter', 'alphajitter')

        alpha = paramdict[alphajitterkey].get_value()
        minjitter = paramdict[minjitterkey].get_value()

        return alpha * (data['data']['rhk'] + 5.0) + minjitter

    elif data['jittermodel'] == 'linear_halpha':
        minjitterkey = jitterkey.replace('jitter', 'minjitter')
        alphajitterkey = jitterkey.replace('jitter', 'alphajitter')

        try:
            alpha = paramdict[alphajitterkey].get_value()
        except KeyError:
            # Use same alpha for SOPHIEm and SOPHIEp
            alpha = paramdict[alphajitterkey.replace('SOPHIEm',
                                                     'SOPHIEp')].get_value()
        minjitter = paramdict[minjitterkey].get_value()

        return alpha * (data['data']['halpha'] -
                        data['data']['halpha'].min()) + minjitter


def chain2inputdict(vddict, index=None):
    """
    Convert a chain dictionary to an input_dict appropiate to construct a model.

    The function returns an input dict with nonesense prior information that
    can be passed to the object builder to construct objects

    :param dict vddict: a dictionary instance with the parameter names and the
     chain traces.

    :param int index: if not None, use this element of chain to build input dict.
    """

    vddict.pop('logL', None)
    vddict.pop('posterior', None)

    # Prepare input dict skeleton
    import string
    input_dict = dict((string.split(s, '_')[0], {}) for s in vddict)

    for i, p in enumerate(vddict.keys()):
        # split object and param name
        obj, par = p.split('_')
        input_dict[obj][par] = [vddict[p][index], 0, 'Uniform', 0.0, 0.0,
                                0.0, 0.0, '']

    return input_dict
