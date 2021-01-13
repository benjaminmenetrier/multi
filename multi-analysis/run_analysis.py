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
from fnmatch import fnmatch
from analysis_tools import *
from analysis_1D_plot import *
from matrix_plot import *

#---------------------------
# Default values for the parameters in namelist:

# Solver:
nm = 3
method = '"theoretical", "standard", "alternative"'
na = 2
algo = '"lanczos", "planczosif"'
no = 2
ni = 4
lmp_mode = '"none"'
test_ortho = "F"
shutoff_type = 0
shutoff_value = 1e-2
interp_method = '"spectral"'
projective_Bmatrix = "T"
# Resolution:
nx = "11,31,51,101"
ny = "11,31,51,101"
# Observations:
nobs = 100
sigma_obs = 0.1
# Background:
sigmabvar = 0.0
Lb = 0.12
spvarmin = 1.0e-5
# Miscellanous:
new_seed = "F"
filename = '"output.nc"'

# if True, the results already existing for the same set of parameters
# will be rewritten:
rewrite_results = True
#---------------------------

directories = {}

# Get the current working directory:
cwd = os.getcwd()
directories['run_analysis'] = cwd

os.chdir('./analysis_results')
cwd = os.getcwd()
directories['analysis_results'] = cwd

os.chdir('../..')
cwd=cwd=os.getcwd()
directories['multi'] = cwd

# Compile the code:
directories['build'] = os.path.join(directories['multi'] + '/build')
os.chdir(directories['build'])
os.system('make')

# Path to the executable from multi:
exec_command = os.path.join(directories['multi'] + '/build/bin/multi ')

# Location of the raw output of the code:
code_output = os.path.join(directories['multi'] + '/output.nc')
#---------------------------

# Analysis results paths:
out_dirs = []

# Root directory of the results of the analysis:
results_dir_root = os.path.join(directories['analysis_results'] + '/cubic_H_test_fct')
out_dirs.append(results_dir_root)

# Raw results of the analysis: 
res_dir_raw = os.path.join(results_dir_root + '/raw_results')
out_dirs.append(res_dir_raw)

# Results for the comparision between LMP modes:
# compare_methods_dir = os.path.join(results_dir_root + '/compare_methods/')
# out_dirs.append(compare_methods_dir)

# Creates the results directories:
for dir in out_dirs:
    if not os.path.exists(dir):
        os.mkdir(dir)

os.chdir(directories['run_analysis'])    

################################################################################
# Results directories:
res_dir_list = []

# Store the outer iteraions for plotting:
outer_iterations_list = []

# Loop over the parametrizations and run the code:
i_no = [3]
i_ni = [6]

i_nx = ['21,51,101']

i_nobs = [1000]
i_sigma_obs = [0.01]

i_sigmabvar = [0.]
i_Lb = [0.1]

i_interp_method = ['"spectral"']#['"bilinear"','"spectral"','"nearest"']
i_project_B = ["T"]
i_test_ortho = ["T"]

iter_params = itertools.product(i_no, i_ni, i_nx, i_nobs, i_sigma_obs,
                                i_sigmabvar, i_Lb, i_interp_method,
                                i_project_B, i_test_ortho, repeat=1)

for no, ni, nx, nobs, sigma_obs, sigmabvar, Lb, interp_method, projective_Bmatrix, test_ortho in iter_params:
    # square grid:
    ny=nx
    
    # outer iteraions for plotting:
    outer_iterations = []
    for io in range(no+1):
        outer_iterations.append((ni)*io)
    outer_iterations_list.append(outer_iterations)
        
    # Create the results directory:
    name_string = f'res_no{no}_ni{ni}'
    name_string += '_lmp-{}'.format(lmp_mode.replace('"',''))
    #name_string += '_met-{}'.format(method.replace('"',''))
    name_string += f'_nx{nx}_ny{ny}_n-obs{nobs}_sigmaobs{sigma_obs}'
    name_string += f'_sigbvar{sigmabvar}_Lb{Lb}'
    name_string += f'_orth{test_ortho}_shut_{shutoff_type}-{shutoff_value}'
    name_string += '_interp-{}'.format(interp_method.replace('"',''))
    name_string += f'_proj{projective_Bmatrix}'
    name_string = name_string.replace(',','-').replace(' ','')
    
    res_dir = os.path.join(res_dir_raw + '/' +  name_string)
    res_dir_list.append(res_dir)
    
    # Rewrite the results or not:
    if not rewrite_results and os.path.exists(res_dir):
        continue
    
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)
        
    parameters={}
    parameters['solver'] = [nm, method, na, algo, no, ni, lmp_mode, test_ortho,
                            shutoff_type, shutoff_value,
                            interp_method, projective_Bmatrix]
    parameters['resolution'] = [nx, ny]
    parameters['obs'] = [nobs, sigma_obs]
    parameters['background'] = [sigmabvar, Lb, spvarmin]
    parameters['miscellanous'] = [new_seed, filename]
    
    namelist_write(parameters,directories['multi'])
        
    # Run the code:
    os.chdir(directories['multi'])
    #os.system('echo "namelist" | '+exec_command)
    os.system(exec_command + " namelist")
    copyfile("namelist", os.path.join(res_dir + "/namelist"))                                    
    os.chdir(directories['run_analysis'])
    # Copy the results:
    copyfile(code_output, os.path.join(res_dir + "/output.nc"))
                                    
################################################################################
# Analysis:
print('Starting analysis:')

# # Loop over the results directories produced:
for r, res_dir in enumerate(res_dir_list):
    
    # Get the data:
    ds = netcdf_extract(res_dir)
    
    # Plots in observation space:
    obs_plot(ds, res_dir)
    hxg_plot(ds, res_dir)
    innovation_plot(ds, res_dir)
    x_true = np.array(ds['xt'][:])
    out_name = os.path.join(res_dir + '/xt.png')
    field_plot(x_true, out_name)
    
    # # Plots comparision between lanczos and planczosif:
    lanczos_vs_planczosif_plot(ds, res_dir, outer_iterations_list[r])
    lanczos_vs_planczosif_2D_outer(ds, res_dir)
    # Plots in model space:
    # At outer loop level:
    for met in ds.groups:
        met_dir = os.path.join(res_dir + f'/{met}')
        if not os.path.exists(met_dir):
            os.mkdir(met_dir)
            
        for io in ds[met].groups:
            
            x_coord = np.array(ds[met][io]['x_coord'])
            y_coord = np.array(ds[met][io]['y_coord'])

            background_fields = ['sigmab', 'dirac_cov', 'dirac_cor', 'dirac_cov_bis',
                      'dirac_cor_bis', 'xb']
            
            background_dir = os.path.join(met_dir + '/background')
            if not os.path.exists(background_dir):
                os.mkdir(background_dir)
            
            for field in background_fields:
                try:
                    matrix = np.array(ds[met][io][field][:])
                    out_name = os.path.join(background_dir + f'/{field}_{io}')
                    field_plot(matrix, out_name)
                except:
                    print('Cannot plot matrix for ', met, io, field)

                # At algorithms level:
                for algo in ds[met][io].groups:
                    algo_dir = os.path.join(met_dir + f'/{algo}')
                    if not os.path.exists(algo_dir):
                        os.mkdir(algo_dir)
                
                    algo_fields = ['xg']
                    for field in algo_fields :
                        try:
                            matrix = np.array(ds[met][io][algo][field][:])
                            out_name = os.path.join(algo_dir + f'/{algo}_{field}_{io}')
                            field_plot(matrix, out_name)
                        except:
                            print('Cannot plot matrix for ', met, io, field)
                            
                        # At inner loop level:
                        inner_fields = ['dx']
                        
                        inner_dir = os.path.join(algo_dir + f'/inner_loops')
                        if not os.path.exists(inner_dir):
                            os.mkdir(inner_dir)
                
                        for field in inner_fields:
                            inner_field_dir = os.path.join(inner_dir + f'/{field}')
                            if not os.path.exists(inner_field_dir):
                                os.mkdir(inner_field_dir)
                    
                            for ii in range(ds[met][io][algo].dimensions['nimax'].size):
                                try:
                                    dx = np.array(ds[met][io][algo][field][ii][:])
                                    out_name = os.path.join(inner_field_dir + f'/{algo}_{field}_{io}_inner_{ii}')
                                    field_plot(dx, out_name)
                                except:
                                    print('Cannot plot matrix for ', met, io, field)    
                                                
################################################################################
# Comparision between the different methods:

for r, res_dir in enumerate(res_dir_list):
    #try:
    if True:                        
        ds = netcdf_extract(res_dir)
        # Comparision for 1D variables (cost functions, rho, beta ...):
        compare_methods_plot(ds, outer_iterations_list[r], res_dir)

        compare_methods_plot2(ds, outer_iterations_list[r], res_dir)

        # Comparision for 2D variables at outer loop level:
        compare_methods_2D_outer(ds, res_dir)
    #except:
    #    print("Cannot compare methods: the following file does not exist:\n",res)
print('analysis completed')    
################################################################################
