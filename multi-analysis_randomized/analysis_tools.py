#!/usr/bin/env python

################################################################################
# analysis_tools.py

# purpose: Contains usefull tools.

# Author: Nicolas Baillot d'Etivaux

################################################################################
# Imported packages:
import os
import netCDF4 as nc
import numpy as np
################################################################################

################################################################################
################################################################################
def init_state_gen(full_res, nobs_max, directory):
    """Generates random vectors to produce the truth and the background initial states at full resolution
    """

    # Truth and background:
    npts = full_res*full_res
    xt_full = np.random.normal(0, 1, npts)
    xb_full = np.random.normal(0, 1, npts)
    
    #print('xt_full shape:', np.shape(xt_full))
    #print('xb_full shape:', np.shape(xb_full))
    print(nobs_max)
    
    xt_full_file = open(os.path.join(directory,"xt_full.dat"), "w")
    xb_full_file = open(os.path.join(directory,"xb_full.dat"), "w")
    
    for xt_val in xt_full:
        xt_full_file.write(f'{xt_val} ')
    xt_full_file.write('\n')    
    xt_full_file.close()

    for xb_val in xb_full:
        xb_full_file.write(f'{xb_val} ')    
    xb_full_file.write('\n')
    xb_full_file.close()

    # Observations:
    x_obs = []
    y_obs = []
    obs_err = []
    for obs in range(nobs_max):
        x_obs.append(np.random.random())
        y_obs.append(np.random.random())
        
    obs_err = np.random.normal(0, 1, nobs_max)
        
    x_obs_file = open(os.path.join(directory,"x_obs.dat"), "w")
    y_obs_file = open(os.path.join(directory,"y_obs.dat"), "w")
    obs_err_file = open(os.path.join(directory,"obs_err.dat"), "w")

    for x in x_obs:
        x_obs_file.write(f'{x} ')
    x_obs_file.write('\n')
    x_obs_file.close()    

    for y in y_obs:
        y_obs_file.write(f'{y} ')
    y_obs_file.write('\n')
    y_obs_file.close()

    for obs in obs_err:
        obs_err_file.write(f'{obs} ')
    obs_err_file.write('\n')
    obs_err_file.close()

    #-------------------------------------------------------------------------------
    # Old trial with netcdf files:
    

    # Read the netcdf file containing xb, xt and observations:
    #file = sys.argv[1]
    
    #ds=nc.Dataset("init_state.nc","r")
    
    # # Get the truth:
    # xt_full = np.array(ds['xt'])
    
    # # Get the observations:
    # x_obs = np.array(ds['x_obs'])
    # y_obs = np.array(ds['y_obs'])
    # obs_val = np.array(ds['obs_val'])
    
    # # Get the last iteration (full resolution):
    # for io in ds['theoretical'].groups:
    #     no = io
    
    # # Get the background:    
    # xb_full = np.array(ds['theoretical'][no]['xb'])
    #--------------------------------------------------------------------------------
    
################################################################################
################################################################################
def netcdf_extract(res_dir):
    """ Extract the netcdf data from the result directory res_dir:
        Args: 
            res_dir: result directory in which is stored the netcdf output file.
        Return: 
            ds: netCDF4 dataset containing the results.
    """
    ds=nc.Dataset(res_dir+"/output.nc","r")
    return ds
################################################################################
################################################################################
# Write the parameters in namelist file:
def namelist_write(parameters, directory, rand_seed):
    """Write the namelist file according to the given parameters.
    """
    
    parameters_name={}
    parameters_name['solver'] = ['nm', 'method', 'na', 'algorithm', 'no', 'ni', 'lmp_mode', 'test_ortho',
                                 'shutoff_type', 'shutoff_value',
                                 'interp_method', 'projective_Bmatrix',]
    parameters_name['resolution'] = ['nx', 'ny']
    parameters_name['obs'] = ['nobs', 'sigma_obs', 'Hnl_coeff']
    parameters_name['background'] = ['sigmabvar', 'Lb', 'spvarmin']
    parameters_name['miscellanous'] = ['new_seed', 'filename']

    namelist = open(os.path.join(directory + f"/namelist_{rand_seed}"),"w")

    # Solver:
    namelist.write("&solver\n")
    for p,par in enumerate(parameters['solver']):
        namelist.write(parameters_name['solver'][p] + " = {}\n".format(par))
        print(parameters_name['solver'][p] + " = {}\n".format(par))
    namelist.write("/\n\n")
    # Resolutions:
    namelist.write("&resolutions\n")
    for p,par in enumerate(parameters['resolution']):
        namelist.write(parameters_name['resolution'][p] + " = {}\n".format(par))
    namelist.write("/\n\n")
    # Observations:
    namelist.write("&observations\n")
    for p,par in enumerate(parameters['obs']):
        namelist.write(parameters_name['obs'][p] + " = {}\n".format(par))
    namelist.write("/\n\n")
    # Background:
    namelist.write("&background\n")
    for p,par in enumerate(parameters['background']):
        namelist.write(parameters_name['background'][p] + " = {}\n".format(par))
    namelist.write("/\n\n")
    # Miscellanous:
    namelist.write("&miscellanous\n")
    for p,par in enumerate(parameters['miscellanous']):
        namelist.write(parameters_name['miscellanous'][p] + " = {}\n".format(par))
    namelist.write("/\n")
    
    namelist.close()
################################################################################
