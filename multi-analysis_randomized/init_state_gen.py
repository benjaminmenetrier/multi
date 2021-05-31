#!/usr/bin/env python

################################################################################
# xb_xt_obs_gen.py

# purpose: Generate full resolution background, truth and obs state.

# Author: Nicolas Baillot d'Etivaux

################################################################################
# Imported packages:
import numpy as np
import sys
import netCDF4 as nc

print('Exiting: file already exists, are you sure to create a new file ?')
sys.exit()

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

def init_state_gen(full_res):
    """Generates random vectors to produce the truth and the background initial states at full resolution
    """
    xt_full = np.random.normal(0, 1, full_res)
    xb_full = np.random.normal(0, 1, full_res)
    
    #print('xt_full shape:', np.shape(xt_full))
    #print('xb_full shape:', np.shape(xb_full))
    
    xt_full_file = open("../xt_full.dat", "w")
    xb_full_file = open("../xb_full.dat", "w")
    
    for xt_val in xt_full:
        xt_full_file.write(f'{xt_val} ')
    xt_full_file.close()

    for xb_val in xb_full:
        xb_full_file.write(f'{xb_val} ')    
    xb_full_file.close()
