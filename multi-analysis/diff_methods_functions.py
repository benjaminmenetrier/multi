#!/usr/bin/env python

################################################################################
# multi-analysis.py

# purpose: Contains functions for the methods comparision analysis.
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
################################################################################
def netcdf_extract(output):
    """Extract the netCDF data from output file:
    """
    # option 'r' or 'rb' can cause errors: (fix later)
    ds=nc.Dataset(output,"r")
    return ds
################################################################################
################################################################################
def walkers_create(nwalkers,parameters,parameters_to_sample,verb):
    """Generates the walkers in the parameter space to sample:
    """
    all_walkers=[]
    for w in range(nwalkers):
        walker_values=[]
        for key in parameters:
            if key in parameters_to_sample:
                # Sample the parameter in its range of variation:
                par_min=parameters[key]['min']
                par_max=parameters[key]['max']
                par_type=parameters[key]['type']
                if par_type=='float':
                    p=np.random.uniform(par_min, par_max)
                elif par_type=='int':
                    p=int(np.random.uniform(par_min, par_max))
                elif par_type=='geom':
                    print("do the geometry parameters later")
                else:
                    print("type error in walkers_create")
                    break
                walker_values.append(p)
        all_walkers.append(walker_values)
    return all_walkers
################################################################################
################################################################################
# Make it later:
# def walkers_continue(state):
#     """Take an existing state of the chain and reuse it as first state for the run.
#     """
################################################################################
################################################################################
def namelist_and_output_write(p,parameters,parameters_to_sample,directories,verb):
    """Write the namelist file according to the position of the walker p in the 
       parameter space.
    """
    namelist_id="namelist_"
    output_id=os.path.join(directories['outputs']+'output_')
    i=0
    for key in parameters:
        #namelist_id+=key+"_"
        #output_id+=key+"_"
        if key in parameters_to_sample:
            par_type=parameters[key]['type']
            if par_type=='float':
                namelist_id+=str(p[i])+"_"
                output_id+=str(p[i])+"_"
            elif par_type=='int':
                # /!\ keep floats for the name to avoid conflicts:
                namelist_id+=str(p[i])+"_"
                output_id+=str(p[i])+"_"
            elif par_type=='geom':
                print("do the geometry parameters later")
            else:
                print("type error in walkers_create")
                break
            i+=1
        elif key=='filename':
            pass
        else:
            namelist_id+=str(parameters[key]['val'])+"_"
            output_id+=str(parameters[key]['val'])+"_"
    namelist_id=namelist_id.replace('"','')
    output_id=output_id.replace('"','')
    output_id+='.nc'
    
    # Write the parameters in namelist file:
    namelist=open(os.path.join(directories["namelists"]+namelist_id),"w")
    # Solver:
    namelist.write("&solver\n")
    namelist.write("no = {}\n".format(parameters['no']['val']))
    namelist.write("ni = {}\n".format(parameters['ni']['val']))
    namelist.write("lmp_mode = {}\n".format(parameters['lmp_mode']['val']))
    namelist.write("test_ortho = {}\n".format(parameters['test_ortho']['val']))
    namelist.write("shutoff_type = {}\n".format(parameters['shutoff_type']['val']))
    namelist.write("shutoff_value = {}\n".format(parameters['shutoff_value']['val']))
    namelist.write("method = {}\n".format(parameters['method']['val']))
    namelist.write("transitive_interp = {}\n".format(parameters['transitive_interp']['val']))
    namelist.write("projective_Bmatrix = {}\n".format(parameters['projective_Bmatrix']['val']))
    namelist.write("/\n\n")
    # Resolutions:
    namelist.write("&resolutions\n")
    namelist.write("nx = {}\n".format(parameters['nx']['val']))
    namelist.write("ny = {}\n".format(parameters['ny']['val']))
    namelist.write("/\n\n")
    # Observations:
    namelist.write("&observations\n")
    namelist.write("nobs = {}\n".format(int(parameters['nobs']['val'])))
    namelist.write("sigma_obs = {}\n".format(parameters['sigma_obs']['val']))
    namelist.write("/\n\n")
    # Background:
    namelist.write("&background\n")
    namelist.write("sigmabvar = {}\n".format(parameters['sigmabvar']['val']))
    namelist.write("Lb = {}\n".format(parameters['Lb']['val']))
    namelist.write("spvarmin = {}\n".format(parameters['spvarmin']['val']))
    namelist.write("/\n\n")
    # Miscellanous:
    namelist.write("&miscellanous\n")
    namelist.write("new_seed = {}\n".format(parameters['new_seed']['val']))
    
    # The same name format is used for the corresponding ouput file:
    namelist.write("filename = {}\n".format('"'+output_id+'"'))
    namelist.write("/\n")
    namelist.close()

    namelist_and_output={'namelist':namelist_id, 'output':output_id}
    return namelist_and_output
################################################################################ 
################################################################################
