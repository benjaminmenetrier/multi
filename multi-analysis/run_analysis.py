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
from multi_analysis import multi_plot
from multi_analysis import lmp_compare
from multi_analysis import diff_plot
from fnmatch import fnmatch

# Default values for the parameters:
nobs=128
no=4 # error when no>4 ?
ni=5
lmp_mode='ritz'
full_res='T'
new_seed='F'
sigma_obs=0.1
sigmabvar=0.0
Lb=0.005

# if True, the results already existing will be rewritten
rewrite_results=False

# Paths and executable file:
path_to_code='..'
exec_code='./run/multi '

# Output paths:
out_dirs=[]

# Where are written the raw outputs of the code:
code_output='../results'
out_dirs.append(code_output)

# Global results of the analysis:
results_dir_root='./analysis_results/'
out_dirs.append(results_dir_root)

# Raw results of the analysis: 
res_dir_raw=results_dir_root+'raw_results/'
out_dirs.append(res_dir_raw)

# Results for the comparision between LMP modes:
res_dir_lmp_compare=results_dir_root+'lmp_compare/'
out_dirs.append(res_dir_lmp_compare)
    
# Results for the evolution of the difference in J vs nobs
res_dir_diff_vs_nobs=results_dir_root+'diff_vs_nobs/'
out_dirs.append(res_dir_diff_vs_nobs)

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
    for nobs in [128,2048]:
        for no in [3]:
            for ni in [4]:
                for sigma_obs in [0.1]:
                    for sigmabvar in [0]:
                        for Lb in [0.005]:
                            # Outer iteraions for plotting:
                            outer_iterations=[]
                            for io in range(no):
                                outer_iterations.append((ni+1)*io)
                            outer_iterations_list.append(outer_iterations)

                            # parameters of the code:
                            parameters=[nobs, no, ni, lmp_mode, full_res, new_seed, sigma_obs, sigmabvar, Lb]

                            # Create the results directory:
                            res_dir=res_dir_raw+'res_nobs{}_no{}_ni{}_lmp-{}_sigmao{}_sigmab{}_Lb{}_reso-{}/'.format(nobs,no,ni,lmp_mode,sigma_obs,sigmabvar,Lb,full_res)
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
#multi_plot(res_dir_list,outer_iterations_list)
################################################################################



################################################################################
# Comparision of LMP methods:

# Output filenames:
out_names=[]
# Store the outer_itertaions:
outer_iterations_list_tmp=[]
# Store the results files to compare:
lmp_to_compare=[]

for r,res_dir in enumerate(res_dir_list):
    if 'ritz' in res_dir:
        res_tmp=res_dir.split('ritz')
        res_tmp1,res_tmp2=res_tmp[0],res_tmp[1]
        # Store the output files names
        out_name=res_tmp1+'compare'+res_tmp2
        out_name=out_name.split(res_dir_raw)[1]
        out_name=res_dir_lmp_compare+out_name[:-1]+'.png'
        out_names.append(out_name)
        # Store the outer iterations:
        outer_iterations_list_tmp.append(outer_iterations_list[r])
        # Store the results files to compare:
        lmp_to_compare_tmp=[]
        for lmp in ['ritz','spectral','none']:
            lmp_to_compare_tmp.append(res_tmp1+lmp+res_tmp2)
        lmp_to_compare.append(lmp_to_compare_tmp)

# Plots the comparision of LMP methods:
lmp_compare(out_names,lmp_to_compare,outer_iterations_list)
################################################################################


################################################################################
# Comparision of LMP methods:

# Output filenames:
out_names=[]
# Store the outer_itertaions:
outer_iterations_list_tmp=[]
# Store the results files to compare:
lmp_to_compare=[]

for r,res_dir in enumerate(res_dir_list):
    if 'ritz' in res_dir:
        res_tmp=res_dir.split('ritz')
        res_tmp1,res_tmp2=res_tmp[0],res_tmp[1]
        # Store the output files names
        out_name=res_tmp1+'compare'+res_tmp2
        out_name=out_name.split(res_dir_raw)[1]
        out_name=res_dir_lmp_compare+out_name[:-1]+'.png'
        out_names.append(out_name)
        # Store the outer iterations:
        outer_iterations_list_tmp.append(outer_iterations_list[r])
        # Store the results files to compare:
        lmp_to_compare_tmp=[]
        for lmp in ['ritz','spectral','none']:
            lmp_to_compare_tmp.append(res_tmp1+lmp+res_tmp2)
        lmp_to_compare.append(lmp_to_compare_tmp)

# Plots the comparision of LMP methods:
lmp_compare(out_names,lmp_to_compare,outer_iterations_list)
################################################################################


################################################################################
# Output filenames:
out_names=[]
# Store the outer_itertaions:
outer_iterations_list_tmp=[]
# Store the results files to compare:
diff_vs_nobs=[]

for r,res_dir in enumerate(res_dir_list):
        res_tmp=res_dir.split('nobs')
        res_tmp1,res_tmp2=res_tmp[0],res_tmp[1]
        # Store the output files names
        out_name=res_tmp1+'compare'+res_tmp2
        out_name=out_name.split(res_dir_raw)[1]
        out_name=res_dir_diff_vs_nobs+out_name[:-1]+'.png'
        out_names.append(out_name)
        # Store the outer iterations:
        outer_iterations_list_tmp.append(outer_iterations_list[r])
        # Store the results files to compare:
        diff_vs_nobs_tmp=[]
        for nobs in [128,2048]:
            diff_vs_nobs_tmp.append(res_tmp1+str(nobs)+res_tmp2)
        diff_vs_nobs.append(diff_vs_nobs_tmp)
        print(diff_vs_nobs)
# Plots the comparision of LMP methods:
diff_plot(out_names,diff_vs_nobs,outer_iterations_list)        
################################################################################


