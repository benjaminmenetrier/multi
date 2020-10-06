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
no=4 # error when no>4 ?
ni=5
lmp_mode='ritz'
sigma_obs=0.1
sigmabvar=0.0
Lb=0.005
full_res='F'
new_seed='F'
gp_from_sp='T' # rajouter cette option

# not applied yet:
nobs=nres
gp_from_sp='T'
res_changing=2



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
results_dir_root='./analysis_results_res2_more_complex/'
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


# # Loop over the parameters and run the code:
# for lmp_mode in ['ritz','spectral','none']:
#     for nres in [128,2048]:
#         for no in [4,6]:
#             for ni in [2,6]:
#                 for sigma_obs in [0.01,0.1]:
#                     for sigmabvar in [0.01,0.1]:
#                         for Lb in [0.001,0.1]:

# Loop over the parameters and run the code:
for lmp_mode in ['ritz','spectral','none']:
    for nres in [2048]:
        for no in [4]:
            for ni in [1,2,6,8]:
                for sigma_obs in [0.01]:
                    for sigmabvar in [0.1]:
                        for Lb in [0.001]:

                            # Outer iteraions for plotting:
                            outer_iterations=[]
                            for io in range(no):
                                outer_iterations.append((ni+1)*io)
                            outer_iterations_list.append(outer_iterations)

                            # parameters of the code:
                            parameters=[nres, no, ni, lmp_mode, sigma_obs, sigmabvar, Lb, full_res, new_seed]

                            # Create the results directory:
                            res_dir=res_dir_raw+'res_nres{}_no{}_ni{}_lmp-{}_sigmao{}_sigmab{}_Lb{}_reso-{}/'.format(nres,no,ni,lmp_mode,sigma_obs,sigmabvar,Lb,full_res)
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



# ################################################################################
# # Comparision of LMP methods:

# # Output filenames:
# out_names=[]
# # Store the outer_itertaions:
# outer_iterations_list_tmp=[]
# # Store the results files to compare:
# lmp_to_compare=[]

# for r,res_dir in enumerate(res_dir_list):
#     if 'ritz' in res_dir:
#         res_tmp=res_dir.split('ritz')
#         res_tmp1,res_tmp2=res_tmp[0],res_tmp[1]
#         # Store the output files names
#         out_name=res_tmp1+'compare'+res_tmp2
#         out_name=out_name.split(res_dir_raw)[1]
#         out_name=res_dir_lmp_compare+out_name[:-1]+'.png'
#         out_names.append(out_name)
#         # Store the outer iterations:
#         outer_iterations_list_tmp.append(outer_iterations_list[r])
#         # Store the results files to compare:
#         lmp_to_compare_tmp=[]
#         for lmp in ['ritz','spectral','none']:
#             lmp_to_compare_tmp.append(res_tmp1+lmp+res_tmp2)
#         lmp_to_compare.append(lmp_to_compare_tmp)

# # Plots the comparision of LMP methods:
# lmp_compare(out_names,lmp_to_compare,outer_iterations_list)
# ################################################################################


################################################################################
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


