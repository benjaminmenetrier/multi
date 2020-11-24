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

#---------------------------
# Default values for the parameters in namelist:
no=2
ni=4

lmp_mode='"none"'
test_ortho=".false."
shutoff_type=0
shutoff_value=0.0
method='"theoretical"'
transitive_interp=".true."
projective_Bmatrix=".true."

nx="101,121,151"
ny="101,121,151"

nobs=3000
sigma_obs=0.1

sigmabvar=0.0
Lb=0.12
spvarmin=1.0e-5

new_seed=".false."
filename='"output"'


# if True, the results already existing for the same set of parameters
# will be rewritten:
rewrite_results=False
#---------------------------

# Path to the executable from multi:
exec_command='./build/bin/multi '

# Location of the raw outputs of the code:
code_output='../output.nc'

# Analysis results paths:
out_dirs=[]

# Root of the results of the analysis:
results_dir_root='./analysis_results_dev/'
out_dirs.append(results_dir_root)

# Raw results of the analysis: 
res_dir_raw=results_dir_root+'raw_results/'
out_dirs.append(res_dir_raw)

# Results for the comparision between LMP modes:
compare_methods_dir=results_dir_root+'compare_methods/'
out_dirs.append(compare_methods_dir)

# Creates the results directories:
for dir in out_dirs:
    if not os.path.exists(dir):
        os.mkdir(dir)
    
# Get the current working directory:
cwd=os.getcwd()

################################################################################
# Results directories:
res_dir_list=[]

# Store the outer iteraions for plotting:
outer_iterations_list=[]

# Loop over the parametrizations and run the code:
# (future improvement: use itertools)
for no in [4]:
    for ni in [6]:
        for lmp_mode in ['"none"']:#['"none"','"ritz"','"spectral"']:
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
                                    outer_iterations=[]
                                    for io in range(no+1):
                                        outer_iterations.append((ni)*io)
                                        outer_iterations_list.append(outer_iterations)
                                        
                                    # Create the results directory:
                                    name_string1='res_no{}_ni{}_lmp-{}_met-{}'
                                    name_string1=name_string1.format(no,ni,lmp_mode.replace('"',''),method.replace('"',''))
                                    name_string2='_nx{}_ny{}_n-obs{}_sigmaobs{}'
                                    name_string2=name_string2.format(nx,ny,nobs,sigma_obs)
                                    name_string3='_sigbvar{}_Lb{}'
                                    name_string3=name_string3.format(sigmabvar,Lb)
                                    
                                    res_dir=res_dir_raw+name_string1+name_string2+name_string3
                                    res_dir_list.append(res_dir)
                                    
                                    # Rewrite the results or not:
                                    if not rewrite_results and os.path.exists(res_dir):
                                        continue
                                    if not os.path.exists(res_dir):
                                        os.mkdir(res_dir)
                                        
                                    # Write the parameters in namelist file:
                                    namelist=open("../namelist","w")
                                    # Solver:
                                    namelist.write("&solver\n")
                                    namelist.write("no = {}\n".format(no))
                                    namelist.write("ni = {}\n".format(ni))
                                    namelist.write("lmp_mode = {}\n".format(lmp_mode))
                                    namelist.write("test_ortho = {}\n".format(test_ortho))
                                    namelist.write("shutoff_type = {}\n".format(shutoff_type))
                                    namelist.write("shutoff_value = {}\n".format(shutoff_value))
                                    namelist.write("method = {}\n".format(method))
                                    namelist.write("transitive_interp = {}\n".format(transitive_interp))
                                    namelist.write("projective_Bmatrix = {}\n".format(projective_Bmatrix))
                                    namelist.write("/\n\n")
                                    # Resolutions:
                                    namelist.write("&resolutions\n")
                                    namelist.write("nx = {}\n".format(nx))
                                    namelist.write("ny = {}\n".format(ny))
                                    namelist.write("/\n\n")
                                    # Observations:
                                    namelist.write("&observations\n")
                                    namelist.write("nobs = {}\n".format(nobs))
                                    namelist.write("sigma_obs = {}\n".format(sigma_obs))
                                    namelist.write("/\n\n")
                                    # Background:
                                    namelist.write("&background\n")
                                    namelist.write("sigmabvar = {}\n".format(sigmabvar))
                                    namelist.write("Lb = {}\n".format(Lb))
                                    namelist.write("spvarmin = {}\n".format(spvarmin))
                                    namelist.write("/\n\n")
                                    # Miscellanous:
                                    namelist.write("&miscellanous\n")
                                    namelist.write("new_seed = {}\n".format(new_seed))
                                    namelist.write("filename = {}\n".format(filename))
                                    namelist.write("/\n")
                                    namelist.close()
                                    
                                    # Compile and run the code:
                                    os.chdir("../build")
                                    os.system("ecbuild ..")
                                    os.system("make")
                                    os.chdir(cwd)
                                    os.chdir('..')
                                    os.system(exec_command)
                                    os.chdir(cwd)
                                    # Copy the results:
                                    copyfile(code_output,res_dir+"/output.nc")
                                    
################################################################################
# # Plots:

# # Loop over the results directories produced:
# for r,res_dir in enumerate(res_dir_list):
    
#     # Get the data:
#     ds=netcdf_extract(res_dir)
    
#     # Plots in observation space:
#     obs_plot(ds,res_dir)
#     hxg_plot(ds,res_dir)
#     innovation_plot(ds,res_dir)
    
#     # Plots comparision between lanczos and planczosif:
#     lanczos_vs_planczosif_plot(ds,res_dir,outer_iterations_list[r])
    
#     # Plots in model space:
#     # At outer loop level:
#     for io in ds.groups:
#         x_coord=np.array(ds[io]['x_coord'])
#         y_coord=np.array(ds[io]['y_coord'])
#         for field in ['sigmab','dirac_cov','dirac_cor','xb']:
#             matrix=np.array(ds[io][field][:])
#             out_name=res_dir+'/'+field+'_'+io
#             field_plot(matrix,out_name)
#         for algo in ds[io].groups:
#             for field in ['xg']:
#                 matrix=np.array(ds[io][algo][field][:])
#                 out_name=res_dir+'/'+algo+'_'+field+'_'+io
#                 field_plot(matrix,out_name)
#                 # At inner loop level:
#                 for ii in range(ds[io][algo].dimensions['nimax'].size):
#                     dx=np.array(ds[io][algo]['dx'][ii][:])
#                     out_name=res_dir+'/'+algo+'_dx_'+io+'_'+'inner_'+str(ii)
#                     field_plot(dx,out_name)
# ################################################################################

# Comparision between the different methods:
methods_list=['theoretical','standard','alternative']
for r,res_dir in enumerate(res_dir_list):
    #try:
    aaa=True
    if aaa==True:
        if methods_list[0] in res_dir:
            res_tmp=res_dir.split('theoretical')
            res_tmp1,res_tmp2=res_tmp[0],res_tmp[1]
            # Store the output files names
            compare_methods_out=res_tmp1+'compare'+res_tmp2
            compare_methods_out=compare_methods_out.replace(res_dir_raw,compare_methods_dir)
            if not os.path.exists(compare_methods_out):
                os.mkdir(compare_methods_out)
            compare_methods=[]
            for method in methods_list:
                compare_methods.append(res_tmp1+method+res_tmp2)
            compare_methods_data=[]
            for res in compare_methods:
                ds=netcdf_extract(res)
                compare_methods_data.append(ds)
            # Comparision for 1D variables (cost functions, rho, beta ...):
            #compare_methods_plot(compare_methods_data,methods_list,outer_iterations_list[r],compare_methods_out)
            # Comparision for 2D variables:
            compare_methods_2D_outer(compare_methods_data,methods_list,compare_methods_out)
    #except:
    aaa=False
        #print("Cannot compare methods: the following file does not exist:\n",res)
################################################################################
