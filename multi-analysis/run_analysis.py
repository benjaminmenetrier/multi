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
# from multi_analysis import multi_plot
# from multi_analysis import lmp_compare
# #from multi_analysis import diff_plot
from multi_analysis import *
from fnmatch import fnmatch

# Default values for the parameters:
nres=128
no=2
ni=4
obsdist=4
lmp_mode='ritz'
sigma_obs=0.01
sigmabvar=0.1
Lb=0.005
full_res='F'
new_seed='F'
gp_from_sp='T' # rajouter cette option
shutoff_type=10 # stop criterion: 1-Jb, 2-beta, (else: no stop criterion)
shutoff_value=0.9 # stop criterion value: if Jb: close to 1, if beta: close to 0

# if True, the results already existing will be rewritten
rewrite_results=True

# Paths and executable file:
path_to_code='..'
exec_code='./run/multi '

# Output paths:
out_dirs=[]

# Where are written the raw outputs of the code:
code_output='../results'
out_dirs.append(code_output)

# Roots of the results of the analysis:
results_dir_root='./analysis_results_dev/'
out_dirs.append(results_dir_root)

# Raw results of the analysis: 
res_dir_raw=results_dir_root+'raw_results/'
out_dirs.append(res_dir_raw)

# Results for the comparision between LMP modes:
res_dir_lmp_compare_J=results_dir_root+'lmp_compare_J/'
out_dirs.append(res_dir_lmp_compare_J)

# Results for the check of second-level lmp:
res_dir_lmp_check=results_dir_root+'check_second_level_lmp/'
out_dirs.append(res_dir_lmp_check)

# Results for rho:
res_dir_lmp_compare_rho=results_dir_root+'lmp_compare_rho/'
out_dirs.append(res_dir_lmp_compare_rho)

# Results for beta:
res_dir_lmp_compare_beta=results_dir_root+'lmp_compare_beta/'
out_dirs.append(res_dir_lmp_compare_beta)

# Results for the innovation:
res_dir_outer_vectors=results_dir_root+'outer_vectors/'
out_dirs.append(res_dir_outer_vectors)

# Results for the evolution of the difference in J vs nres
res_dir_diff_vs_nres=results_dir_root+'diff_vs_nres/'
out_dirs.append(res_dir_diff_vs_nres)

# Creates the results directories:
for dir in out_dirs:
    if not os.path.exists(dir):
        os.mkdir(dir)
    
# Get the current working directory:
cwd=os.getcwd()

# Compile the code:
os.chdir(path_to_code)
os.system('make')
os.chdir(cwd)

################################################################################
# Run the code, and store the results:

# Store the results directories produced in the following loop:
res_dir_list=[]
# Outer iteraions for plotting:
outer_iterations_list=[]

# Loop over the parameters and run the code:
for lmp_mode in ['ritz','spectral','none']:
    for nres in [128]:
        for no in [3]:
            for ni in [4]:
                for obsdist in [4]:
                    for sigma_obs in [0.1]:
                        for sigmabvar in [0.1]:
                            for Lb in [1]:

                                # Outer iteraions for plotting:
                                outer_iterations=[]
                                for io in range(no):
                                    outer_iterations.append((ni+1)*io)
                                outer_iterations_list.append(outer_iterations)

                                # parameters of the code:
                                parameters=[nres, no, ni, obsdist, lmp_mode, sigma_obs,
                                            sigmabvar, Lb, full_res, new_seed,
                                            shutoff_type,shutoff_value]

                                # Create the results directory:
                                name_string='res_reso{}_no{}_ni{}_obsdist{}_lmp-{}'
                                name_string=name_string.format(nres,no,ni,obsdist,lmp_mode)
                                name_string2='_sigmao{}_sigmab{}_Lb{}_reso-{}/'
                                name_string2=name_string2.format(sigma_obs,sigmabvar,Lb,full_res)
                                res_dir=res_dir_raw+name_string+name_string2
                                res_dir_list.append(res_dir)

                                if not rewrite_results and os.path.exists(res_dir):
                                    continue
                                if not os.path.exists(res_dir):
                                    os.mkdir(res_dir)

                                # define the command line to run the code:
                                arguments=''
                                for par in parameters:
                                    arguments+=' {} '.format(par)
                                exec_command ='echo '+ arguments + ' | ' + exec_code
                                print("\n",exec_command,"\n")

                                # Run the code and save the results:
                                os.chdir(path_to_code)
                                os.system(exec_command)
                                os.chdir(cwd)
                                copy_tree(code_output,res_dir)
################################################################################



################################################################################
# Plot the results of the code:
multi_plot(res_dir_list,outer_iterations_list)
################################################################################



# ##############################################################################
# Plots the outer vectors for each outer loop:

# Plots the innovation:
#for res_dir in res_dir_list:
#   yo_vs_hxg_plot(res_dir)

# Plots vectors on the grid points:
res_file_names=['lanczos_control_space_outer_grid.dat','PlanczosIF_model_space_outer_grid.dat']
for f,res_file_name in enumerate(res_file_names):

    coord_column=2
    # Plots the guess:
    out_name='_guess.png'
    label=r'$x^g'
    column_of_interest=-1
    for res_dir in res_dir_list:
        results_file=res_dir+res_file_name
        out_file_name=results_file[:-4]+out_name
        print('out_file',out_file_name)
        vec_plot(results_file,column_of_interest,coord_column,label,out_file_name)

    # Plots the background:
    out_name='_background.png'
    label=r'$x^b'
    column_of_interest=-2
    for res_dir in res_dir_list:
        results_file=res_dir+'/'+res_file_name
        out_file_name=results_file[:-4]+out_name
        vec_plot(results_file,column_of_interest,coord_column,label,out_file_name)

    # Plots the increment:
    out_name='_increment.png'
    label=r'$\delta x^b'
    column_of_interest=-3
    for res_dir in res_dir_list:
        results_file=res_dir+'/'+res_file_name
        out_file_name=results_file[:-4]+out_name
        vec_plot(results_file,column_of_interest,coord_column,label,out_file_name)

    # Plots the preconditionned vectors:
    column_of_interest=-4
    if f==0:
        out_name='_dvb.png'
        label=r'$\delta v^b'
    elif f==1:
        out_name='_dxbbar.png'
        label=r'$\delta \bar{x^b}'
    for res_dir in res_dir_list:
        results_file=res_dir+'/'+res_file_name
        out_file_name=results_file[:-4]+out_name
        vec_plot(results_file,column_of_interest,coord_column,label,out_file_name)


    column_of_interest=-5
    if f==0:
        out_name='_dva_interp.png'
        label=r'$\Pi \delta v^{a}'
    elif f==1:
        out_name='_dxabar_interp.png'
        label=r'$\Pi \delta \bar{x^{a}}'
    for res_dir in res_dir_list:
        results_file=res_dir+'/'+res_file_name
        out_file_name=results_file[:-4]+out_name
        vec_plot(results_file,column_of_interest,coord_column,label,out_file_name)


# Plots vectors in the obs space:
res_file_names=['lanczos_control_space_outer_obs.dat','PlanczosIF_model_space_outer_obs.dat']
for f,res_file_name in enumerate(res_file_names):

    # Plots Hxg:
    out_name='_hxg.png'
    label=r'$H x^g'
    column_of_interest=-3
    for res_dir in res_dir_list:
        results_file=res_dir+'/'+res_file_name
        out_file_name=results_file[:-4]+out_name
        vec_plot(results_file,column_of_interest,coord_column,label,out_file_name)

    # Plots the innovation:
    out_name='_innovation.png'
    #label=r'$y^o_{}-H x^g'
    label=r'$d'
    column_of_interest=-1
    for res_dir in res_dir_list:
        results_file=res_dir+'/'+res_file_name
        out_file_name=results_file[:-4]+out_name
        vec_plot(results_file,column_of_interest,coord_column,label,out_file_name)

    # Plots the obs:
    out_name='_obs.png'
    label=r'$y^o'
    column_of_interest=-2
    for res_dir in res_dir_list:
        results_file=res_dir+res_file_name
        out_file_name=results_file[:-4]+out_name
        vec_plot(results_file,column_of_interest,coord_column,label,out_file_name)

################################################################################
# Comparision of LMP methods:

# Output filenames:
out_names_J=[]
out_names_check=[]
out_names_rho=[]
out_names_beta=[]
# Store the outer_itertaions:
outer_iterations_list_tmp=[]
# Store the results files to compare:
lmp_to_compare=[]
check_second_level_lmp_dirs=[]

for r,res_dir in enumerate(res_dir_list):
    if 'ritz' in res_dir:
        res_tmp=res_dir.split('ritz')
        res_tmp1,res_tmp2=res_tmp[0],res_tmp[1]
        # Store the output files names
        out_name=res_tmp1+'compare'+res_tmp2
        out_name=out_name.split(res_dir_raw)[1]

        out_name_J=res_dir_lmp_compare_J+out_name[:-1]+'.png'
        out_name_check=res_dir_lmp_check+out_name[:-1]+'.png'
        out_name_rho=res_dir_lmp_compare_rho+out_name[:-1]+'.png'
        out_name_beta=res_dir_lmp_compare_beta+out_name[:-1]+'.png'
        out_names_J.append(out_name_J)
        out_names_check.append(out_name_check)
        out_names_rho.append(out_name_rho)
        out_names_beta.append(out_name_beta)
        
        # Store the outer iterations:
        outer_iterations_list_tmp.append(outer_iterations_list[r])

        # Store the results files to compare:
        lmp_to_compare_tmp=[]
        check_second_level_lmp_tmp=[]
        for lmp in ['ritz','spectral','none']:
            lmp_to_compare_tmp.append(res_tmp1+lmp+res_tmp2)
            if not lmp=='none':
                check_second_level_lmp_tmp.append(res_tmp1+lmp+res_tmp2)
        lmp_to_compare.append(lmp_to_compare_tmp)
        check_second_level_lmp_dirs.append(check_second_level_lmp_tmp)
        
# Plots the comparision of LMP methods according to J:
column_of_interest=3
ylabel1=r'$J=J_o+J_b$'
ylabel2=r'$J_{B^{1/2}}-J_{B}$'
lmp_compare(out_names_J,lmp_to_compare,column_of_interest,ylabel1,ylabel2,outer_iterations_list)
# Plots the comparision of LMP methods according to rho:
column_of_interest=6
ylabel1=r'$\rho$'
ylabel2=r'$\rho_{B^{1/2}}-\rho_{B}$'
lmp_compare(out_names_rho,lmp_to_compare,column_of_interest,ylabel1,ylabel2,outer_iterations_list)
# Plots the comparision of LMP methods according to the beta:
column_of_interest=7
ylabel1=r'$\beta$'
ylabel2=r'$\beta_{B^{1/2}}-\rho_{B}$'
lmp_compare(out_names_beta,lmp_to_compare,column_of_interest,ylabel1,ylabel2,outer_iterations_list)
################################################################################

# Check that the difference between the second level preconditionners:
#check_second_level_lmp(out_names_check,check_second_level_lmp_dirs,outer_iterations_list)

################################################################################
# matrix representation:
for res_dir in res_dir_list:
    # For the B matrix:
    results_file='Bdelta_test.dat'
    out_file_name='B_matrix'
    matrix_monitoring(res_dir,results_file,out_file_name)
    # For the H matrix:
    results_file='Hdelta_test.dat'
    out_file_name='H_matrix'
    matrix_monitoring(res_dir,results_file,out_file_name)

################################################################################



# See later the diff_plots:
################################################################################
# Plot the differences as a function of the resolution, see later.
# # Output filenames:
# out_names=[]
# # Store the outer_itertaions:
# outer_iterations_list_tmp=[]
# # Store the results files to compare:
# diff_vs_nres=[]

# for r,res_dir in enumerate(res_dir_list):
#         res_tmp=res_dir.split('nres')
#         res_tmp1,res_tmp2=res_tmp[0],res_tmp[1]
#         # Store the output files names
#         out_name=res_tmp1+'compare'+res_tmp2
#         out_name=out_name.split(res_dir_raw)[1]
#         out_name=res_dir_diff_vs_nres+out_name[:-1]+'.png'
#         out_names.append(out_name)
#         # Store the outer iterations:
#         outer_iterations_list_tmp.append(outer_iterations_list[r])
#         # Store the results files to compare:
#         diff_vs_nres_tmp=[]
#         for nres in [128,2048]:
#             diff_vs_nres_tmp.append(res_tmp1+str(nres)+res_tmp2)
#         diff_vs_nres.append(diff_vs_nres_tmp)
#         #print(diff_vs_nres)
# # Plots the comparision of LMP methods:
# #diff_plot(out_names,diff_vs_nres,outer_iterations_list)        
# ################################################################################


