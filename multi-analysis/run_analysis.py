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
no = 2
ni = 4

lmp_mode = '"none"'
test_ortho = "F"
shutoff_type = 0
shutoff_value = 1e-2
method = '"theoretical"'
transitive_interp = "T"
projective_Bmatrix = "T"

nx = "11,31,51,101"
ny = "11,31,51,101"

nobs = 100
sigma_obs = 0.1

sigmabvar = 0.0
Lb = 0.12
spvarmin = 1.0e-5

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

# Path to the executable from multi:
exec_command = os.path.join(directories['multi'] + '/build/bin/multi ')

# Location of the raw output of the code:
code_output = os.path.join(directories['multi'] + '/output.nc')
#---------------------------

# Analysis results paths:
out_dirs = []

# Root directory of the results of the analysis:
results_dir_root = os.path.join(directories['analysis_results'] + '/analysis_results_no_resolution_changing')
out_dirs.append(results_dir_root)

# Raw results of the analysis: 
res_dir_raw = os.path.join(results_dir_root + '/raw_results')
out_dirs.append(res_dir_raw)

# Results for the comparision between LMP modes:
compare_methods_dir = os.path.join(results_dir_root + '/compare_methods/')
out_dirs.append(compare_methods_dir)

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
# (future improvement: use itertools)
for no in [1,2]:
    for ni in [4,8]:
        for lmp_mode in ['"none"']:
            for method in ['"theoretical"','"standard"','"alternative"']:
                for nx in ['101']:
                    for nobs in [2000]:
                        for sigma_obs in [0.01]:
                            for sigmabvar in [0.1]:
                                for Lb in [0.1]:
                                    # square grid:
                                    ny=nx
                                    #ny='1'
                                    
                                    # outer iteraions for plotting:
                                    outer_iterations = []
                                    for io in range(no+1):
                                        outer_iterations.append((ni)*io)
                                    outer_iterations_list.append(outer_iterations)
                                        
                                    # Create the results directory:
                                    # use f'' syntaxe here!
                                    # name_string1 = 'res_no{}_ni{}_lmp-{}_met-{}'
                                    # name_string1 = name_string1.format(no,ni,lmp_mode.replace('"',''),method.replace('"',''))
                                    # name_string2 = '_nx{}_ny{}_n-obs{}_sigmaobs{}'
                                    # name_string2 = name_string2.format(nx,ny,nobs,sigma_obs)
                                    # name_string3 = '_sigbvar{}_Lb{}'
                                    # name_string3 = name_string3.format(sigmabvar,Lb)
                                    # name_string4 = '_orth{}_shut_{}-{}_trans{}_proj{}'w
                                    # name_string4 = name_string4.format(test_ortho,shutoff_type,shutoff_value,transitive_interp,projective_Bmatrix)
                                    name_string = f'res_no{no}_ni{ni}'
                                    name_string += '_lmp-{}'.format(lmp_mode.replace('"',''))
                                    name_string += '_met-{}'.format(method.replace('"',''))
                                    name_string += f'_nx{nx}_ny{ny}_n-obs{nobs}_sigmaobs{sigma_obs}'
                                    name_string += f'_sigbvar{sigmabvar}_Lb{Lb}'
                                    name_string += f'_orth{test_ortho}_shut_{shutoff_type}-{shutoff_value}'
                                    name_string += f'_trans{transitive_interp}_proj{projective_Bmatrix}'

                                    res_dir = os.path.join(res_dir_raw + '/' +  name_string)
                                    #name_string1+name_string2+name_string3+name_string4)
                                    res_dir_list.append(res_dir)
                                    # Rewrite the results or not:
                                    if not rewrite_results and os.path.exists(res_dir):
                                        continue
                                    if not os.path.exists(res_dir):
                                        os.mkdir(res_dir)
                                        
                                    parameters = [no, ni, lmp_mode, test_ortho, shutoff_type, shutoff_value, method,
                                                  transitive_interp, projective_Bmatrix, nx, ny, nobs,
                                                  sigma_obs, sigmabvar, Lb, spvarmin, new_seed, filename]

                                    namelist_write(parameters,directories['multi'])

                                    # Compile and run the code:
                                    #os.chdir("../build")
                                    # ecbuild not necessary:
                                    #os.system("ecbuild ..")
                                    #os.system("make")
                                    #os.chdir(cwd)

                                    # Run the code:
                                    os.chdir(directories['multi'])
                                    os.system('echo "namelist" | '+exec_command)
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
    obs_plot(ds,res_dir)
    hxg_plot(ds,res_dir)
    innovation_plot(ds,res_dir)
    
    # Plots comparision between lanczos and planczosif:
    lanczos_vs_planczosif_plot(ds,res_dir,outer_iterations_list[r])
    
    # Plots in model space:
    # At outer loop level:
    for io in ds.groups:

        x_coord = np.array(ds[io]['x_coord'])
        y_coord = np.array(ds[io]['y_coord'])

        background_fields = ['sigmab', 'dirac_cov', 'dirac_cor', 'dirac_cov_bis',
                      'dirac_cor_bis', 'xb']

        background_dir = os.path.join(res_dir + '/background')
        if not os.path.exists(background_dir):
            os.mkdir(background_dir)
            
        for field in background_fields:
            matrix = np.array(ds[io][field][:])
            out_name = os.path.join(background_dir + f'/{field}_{io}')
            field_plot(matrix,out_name)

        # At algorithms level:
        for algo in ds[io].groups:
            algo_dir = os.path.join(res_dir + f'/{algo}')
            if not os.path.exists(algo_dir):
                os.mkdir(algo_dir)
                
            algo_fields = ['xg']
            for field in algo_fields :
                matrix = np.array(ds[io][algo][field][:])
                out_name = os.path.join(algo_dir + f'/{algo}_{field}_{io}')
                field_plot(matrix, out_name)

                # At inner loop level:
                inner_fields = ['dx']

                inner_dir = os.path.join(algo_dir + f'/inner_loops')
                if not os.path.exists(inner_dir):
                    os.mkdir(inner_dir)
                    
                for field in inner_fields:
                    inner_field_dir = os.path.join(inner_dir + f'/{field}')
                    if not os.path.exists(inner_field_dir):
                        os.mkdir(inner_field_dir)
                        
                    for ii in range(ds[io][algo].dimensions['nimax'].size):
                        dx = np.array(ds[io][algo][field][ii][:])
                        out_name = os.path.join(inner_field_dir + f'/{algo}_{field}_{io}_{ii}')
                        field_plot(dx, out_name)
                                                
################################################################################
# Comparision between the different methods:

methods_list = ['theoretical', 'standard', 'alternative']

for r, res_dir in enumerate(res_dir_list):
    #try:
    if True:    
        if methods_list[0] in res_dir:
            res_tmp = res_dir.split('theoretical')
            res_tmp1, res_tmp2 = res_tmp[0], res_tmp[1]
            # Store the output files names
            compare_methods_out = res_tmp1 + 'compare' + res_tmp2
            compare_methods_out = compare_methods_out.replace(res_dir_raw, compare_methods_dir)
            if not os.path.exists(compare_methods_out):
                os.mkdir(compare_methods_out)
            for met in methods_list:
                compare_methods_2D_dir = os.path.join(compare_methods_out + "/theoretical_vs_"+met)
                if not os.path.exists(compare_methods_2D_dir):
                    os.mkdir(compare_methods_2D_dir)
                    
            compare_methods=[]
            for method in methods_list:
                compare_methods.append(res_tmp1 + method + res_tmp2)

            compare_methods_data = []
            for res in compare_methods:
                ds = netcdf_extract(res)
                compare_methods_data.append(ds)

            # Comparision for 1D variables (cost functions, rho, beta ...):
            compare_methods_plot(compare_methods_data, methods_list,
                                 outer_iterations_list[r], compare_methods_out)

            compare_methods_plot2(compare_methods_data, methods_list,
                                  outer_iterations_list[r], compare_methods_out)

            # Comparision for 2D variables at outer loop level:
            compare_methods_2D_outer(compare_methods_data, methods_list,
                                    compare_methods_out)
    #except:
    #    print("Cannot compare methods: the following file does not exist:\n",res)
################################################################################
