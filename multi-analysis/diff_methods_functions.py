#!/usr/bin/env python

################################################################################
# multi-analysis.py

# purpose: This code plots the results of the code multi.
# Author: Nicolas Baillot d'Etivaux

################################################################################
# Imported packages:
import os
import sys
import numpy as np
from distutils.dir_util import copy_tree
from shutil import copyfile
from multi_analysis_2D import *
from fnmatch import fnmatch
import netCDF4 as nc

# Allow the use of tex format for the labels:
plt.rcParams.update({"text.usetex": True, "font.size" : 15})

################################################################################
# usefull or not for the next steps ? is there options that we want to manage ?
def netcdf_extract(output):
    """ Extract the netcdf data from the result directory res_dir:
        Args: 
            res_dir: result directory in which is stored the netcdf output file.
        Return: 
            ds: netCDF4 dataset containing the results.
    """
    ds=nc.Dataset(output,"r")
    return ds
################################################################################
def ln_23(p):
    namelist=open("../namelists_tmp/"+str(p[0]),"w")
    print("../namelists_tmp/"+str(p[0]))
    namelist.close()
    #os.remove("../namelists_tmp/"+str(p))
    return 0
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
            namelist_and_output=namelist_write(p,parameters,parameters_to_sample,directories,verb)
            namelists.append(namelist_and_output['namelist']) 
            outputs.append(namelist_and_output['output'])
            print(namelist_and_output)
        #except:
        #    print("Error in namelist_write -- fix later: the origin seems to be in parallelization of writing tasks")
        #    return -np.inf
        os.chdir(directories["code"])
        command_line=os.path.join('echo "'+directories["namelists"]+'/'+namelists[m]+'" | '+exec_command)
        print(command_line)
        os.system(command_line)
        
        #os.remove(os.path.join(directories["namelists"]+namelists[m]))
        
        # Get the results of the code and store them:
        try:
            print(outputs[m])
            ds=netcdf_extract(outputs[m])
            results[method]=ds
            os.remove(outputs[m])
        except:
            print("Error with output of the following command line:")
            print(command_line)
            return -np.inf
        os.chdir(os.path.join(directories["analysis"]))
    #-----------------------------------------------------------------------
    # Compute the difference:
    algo='lanczos'
    cost_func_id='j_nl'
    #try:
    #if True:
        #diff=diff_compute_cost_function(results,algo,methods_list,cost_func_id)
    #except:
        #print("Error in diff_compute with the following parameters:")
        #print('par=',parameters)
        #return -np.inf
    #-----------------------------------------------------------------------
    return -parameters['nobs']['val']
    #return diff
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
# def walkers_create(nwalkers,parameters,parameters_to_sample,verb):
#     """Generates the walkers in the parameter space to sample:
#     """
#     all_walkers=[]
#     for w in range(nwalkers):
#         walker_values=[]
#         for key in parameters:
#             if key in parameters_to_sample:
#                 # Sample the parameter in its range of variation:
#                 par_min=parameters[key]['min']
#                 par_max=parameters[key]['max']
#                 par_type=parameters[key]['type']
#                 if par_type=='float':
#                     p=np.random.uniform(par_min, par_max)
#                 elif par_type=='int':
#                     p=int(np.random.uniform(par_min, par_max))
#                 elif par_type=='geom':
#                     print("do the geometry parameters later")
#                 else:
#                     print("type error in walkers_create")
#                     break
#                 walker_values.append(p)
#             elif key=="file_name":
#                 # /!\ check here if the name of output file is unic enough.
#                 walker_values.append('"outputs_tmp/output_walker{}.nc"'.format(w))    
#             else:
#                 # Set the parameter to its default value:
#                 walker_values.append(parameters[key]['val'])
#         all_walkers.append(walker_values)
#     return all_walkers    
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
def namelist_write(p,parameters,parameters_to_sample,directories,verb):
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
