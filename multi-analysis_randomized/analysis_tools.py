#!/usr/bin/env python

################################################################################
# analysis_tools.py

# purpose: Contains usefull tools.

# Author: Nicolas Baillot d'Etivaux

################################################################################
# Imported packages:
import os
import netCDF4 as nc
################################################################################

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
def namelist_write(parameters,directory):
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

    namelist = open(os.path.join(directory + "/namelist"),"w")

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
