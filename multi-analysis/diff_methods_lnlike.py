#!/usr/bin/env python

################################################################################
# mcmc_lnlike.py

# purpose: Contains the definition of the likelihood for the mcmc.
# Author: Nicolas Baillot d'Etivaux

################################################################################

# Imported packages:
import os
import sys
import numpy as np
from distutils.dir_util import copy_tree
from shutil import copyfile
from fnmatch import fnmatch
import netCDF4 as nc

from diff_methods_functions import *


################################################################################
def ln_23(p, parameters_dict, parameters_to_sample, exec_command,
            directories, verb):
    """Computes the log prob and monitor it:
    """
    res = ln_prob(p, parameters_dict, parameters_to_sample, exec_command,
            directories, verb)
    if not np.isfinite(res):
        print('ln_prob = ', res)
        return -np.inf
    else:
        print('ln_prob = ', res)
        return res
################################################################################
################################################################################
def ln_prob(p, parameters_dict, parameters_to_sample, exec_command,
            directories, verb):
    """Compute the log probability.
    """
    #----------------------------------------------------------------------------
    # Set the parameters:
    parameters = dict(parameters_dict)
    i = 0
    print('p=',p)
    for k, key in enumerate(parameters):
        if key in parameters_to_sample:
            # Check the boundaries:
            par_min = parameters[key]['min']
            par_max = parameters[key]['max']
            if p[i] < par_min or p[i] > par_max:
                print(f'param {i} out of [{par_min}:{par_max}] with value: {p[i]}')
                return -np.inf
            else:
                # Put the sampled value in the parameters:
                par_type = parameters[key]['type']
                if par_type == 'float':
                    parameters[key]['val'] = p[i]
                    i += 1
                elif par_type == 'int':
                    # /!\ keep float to avoid conflicts:
                    parameters[key]['val'] = p[i]
                    i+=1
                elif par_type=='geom':
                    print("do the geometry parameters later")
                    return -np.inf
                else:
                    print("type error in ln_prob")
                    return -np.inf
    #--------------------------------------------------------------------------
    # Run the two methods:

    # Store the results and namelists:
    
    namelist_and_output = namelist_and_output_write(p, parameters, parameters_to_sample, directories, verb)
    namelist = namelist_and_output['namelist']
    output = namelist_and_output['output']
    #except:
    #    print("Error in namelist_write -- fix later: the origin seems to be in parallelization of writing tasks")
    #    return -np.inf
    os.chdir(directories['multi'])
    path_to_namelist = os.path.join(directories["namelists"] + namelist)
    command_line = 'nohup ' + exec_command + path_to_namelist + ' > test.log'
    #command_line = exec_command + path_to_namelist
    if verb:
        print(command_line)
    os.system(command_line)
    
    os.remove(os.path.join(directories["namelists"] + namelist))
    
    # Get the results of the code and store them:
    #if True:
    try:    
        ds = netcdf_extract(output)
        os.remove(output)
        os.chdir(os.path.join(directories["analysis"]))
        #-----------------------------------------------------------------------
        # Compute the difference:
        algo = 'lanczos'
        method = 'alternative'
        cost_func_id = 'j'
        diff = diff_compute_cost_function(ds, algo, cost_func_id)
    except:     
       print("Error in diff_compute with the following parameters:")
       print('par=',parameters)
       return -np.inf
    #-----------------------------------------------------------------------
    return diff
    
################################################################################
################################################################################
def diff_compute_cost_function(ds, algo, cost_func_id):
    """Compute the difference between the costs functions.
    """
    # Get the cost functions:
    cost_function={}
    for m, met in enumerate(ds.groups):
        cost_function[met] = []
        for io in ds[met].groups:
            cost_function[met].append(np.array(ds[met][io][algo][cost_func_id][:]))
        cost_function[met]=np.reshape(cost_function[met],-1)
        
    # Compute the difference and return it:
    diff=0.
    for i in range(len(cost_function['theoretical'])):
        diff += abs(cost_function['alternative'][i] - cost_function['theoretical'][i])
    return diff
################################################################################
