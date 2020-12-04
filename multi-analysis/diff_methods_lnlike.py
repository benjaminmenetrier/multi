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



################################################################################
def ln_23(p):
    namelist=open("../namelists_tmp/"+str(p[0]),"w")
    print("../namelists_tmp/"+str(p[0]))
    namelist.close()
    #os.remove("../namelists_tmp/"+str(p))
    return 0
################################################################################
################################################################################
def ln_prob(p,parameters,parameters_to_sample,methods_list,exec_command,directories,verb):
    """Compute the log probability.
    """
    #----------------------------------------------------------------------------
    # Set the parameters:
    # copy the dict and factorize this function:
    i=0
    for k,key in enumerate(parameters):
        if key in parameters_to_sample:
            # Check the boundaries:
            par_min=parameters[key]['min']
            par_max=parameters[key]['max']
            if p[i]<par_min or p[i]>par_max:
                print('parameter {} sampled out of its range [{}:{}] with value: {}'.format(i, par_min, par_max, p[i]))
                return -np.inf
            else:
                # Put the sampled value in the parameters:
                par_type=parameters[key]['type']
                if par_type=='float':
                    parameters[key]['val']=p[i]
                    i+=1
                elif par_type=='int':
                    # /!\ keep float to avoid conflicts:
                    parameters[key]['val']=p[i]
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
    results={}
    namelists=[]
    outputs=[]
    for m,method in enumerate(methods_list):
        parameters['method']['val']=method
        # write the namelist file:
        #try:
        if True:
            namelist_and_output=namelist_and_output_write(p,parameters,parameters_to_sample,directories,verb)
            namelists.append(namelist_and_output['namelist']) 
            outputs.append(namelist_and_output['output'])
        #except:
        #    print("Error in namelist_write -- fix later: the origin seems to be in parallelization of writing tasks")
        #    return -np.inf
        os.chdir(directories['multi'])
        path_to_namelist=os.path.join(directories["namelists"]+namelists[m])
        print(path_to_namelist)
        command_line=os.path.join('echo "'+path_to_namelist+'" | '+exec_command)
        print(command_line)
        os.system(command_line)
        
        os.remove(os.path.join(directories["namelists"]+namelists[m]))
        
        # Get the results of the code and store them:
        #try:
        if True:
            ds=netcdf_extract(outputs[m])
            results[method]=ds
            os.remove(outputs[m])
        #except:
        #    print("Error with output of the following command line:")
        #    print(command_line)
        #    return -np.inf
        os.chdir(os.path.join(directories["analysis"]))
    #-----------------------------------------------------------------------
    # Compute the difference:
    algo='lanczos'
    cost_func_id='j_nl'
    try:
        diff=diff_compute_cost_function(results,algo,methods_list,cost_func_id)
    except:
       print("Error in diff_compute with the following parameters:")
       print('par=',parameters)
       return -np.inf
    #-----------------------------------------------------------------------
    return diff
################################################################################
################################################################################
def diff_compute_cost_function(results,algo,methods_list,cost_func_id):
    """Compute the difference between the costs functions.
    """
    # Get the cost functions:
    if len(methods_list)!=2:
        print("Error in diff_compute: too much methods to compare")
    #cost_function=[]
    cost_function={}
    for m,method in enumerate(results):
        #cost_function.append([])
        cost_function[method]=[]
        ds=results[method]
        for io in ds.groups:
            #cost_function[m].append(np.array(ds[io][algo][cost_func_id][:]))
            cost_function[method].append(np.array(ds[io][algo][cost_func_id][:]))
        #cost_function[m]=np.reshape(cost_function[m],-1)
        cost_function[method]=np.reshape(cost_function[method],-1)
        
    # Compute the difference and return it:
    diff=0
    for i in range(len(cost_function[methods_list[0]])):
        diff+=abs(cost_function[methods_list[0]][i]-cost_function[methods_list[1]][i])
    return diff
################################################################################
