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
from multi_analysis import compare_plot_N


# Default values for the parameters:
n=128
no=4 # error when no>4 ?
ni=5
lmp_mode='ritz'
full_res='F'
new_seed='T'
sigma_obs=0.1
sigmabvar=0.0
Lb=0.005

# Paths and executable file:
path_to_code='..'
exec_code='./run/multi '
code_output='../results'
results_dir_root='./analysis_results/'

# Get the current working directory:
cwd=os.getcwd()

# Compile the code:
os.chdir(path_to_code)
os.system('make')
os.chdir(cwd)

# Create the results directory:
if not os.path.exists(results_dir_root):
    os.mkdir(results_dir_root)

# outer iteraions for plotting
outer_iterations=[]
for io in range(no):
    outer_iterations.append((ni+1)*io)

# Loop over a given parameter:
for sigmabvar in [0.0]:
    
    parameters=[n, no, ni, lmp_mode, full_res, new_seed, sigma_obs, sigmabvar, Lb]
    print("\n ============== running with sigmabvar={} ============== \n".format(sigmabvar))
    
    # Create the results directory:
    res_dir=results_dir_root+'res_n{}_no{}_ni{}_lmp-{}_sigmao{}_sigmab{}_Lb{}_reso-{}/'.format(n,no,ni,lmp_mode,sigma_obs,sigmabvar,Lb,full_res)
    
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)
        
    # define the command line to run the code:
    arguments=''
    for par in parameters:
        arguments+=' {} '.format(par)
        #print(arguments)
    exec_command ='echo '+ arguments + ' | ' + exec_code
    print(exec_command)
    
    # Run the code and save the results:
    os.chdir(path_to_code)
    os.system(exec_command)
    os.chdir(cwd)
    copy_tree(code_output,res_dir)
    
    # Run the analysis over the results:
    multi_plot(res_dir,outer_iterations)
