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
from multi_analysis import compare_plots_2N


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

# Outer iteraions for plotting:
outer_iterations=[]
for io in range(no):
    outer_iterations.append((ni+1)*io)

# Store the results directories produced in the following loop:
res_dir_list=[]

# Loop over the parameters and run the code:
for lmp_mode in ['ritz','spectral']:
    
    parameters=[n, no, ni, lmp_mode, full_res, new_seed, sigma_obs, sigmabvar, Lb]
    
    # Create the results directory:
    res_dir=results_dir_root+'res_n{}_no{}_ni{}_lmp-{}_sigmao{}_sigmab{}_Lb{}_reso-{}/'.format(n,no,ni,lmp_mode,sigma_obs,sigmabvar,Lb,full_res)
    res_dir_list.append(res_dir)
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

# Run the analysis over the results:
for res_dir in res_dir_list:
    multi_plot(res_dir,outer_iterations)




# Compare spectral and ritz lmp modes:
results_directory_test=[]
labels=[]
legend=[]

for lmp_mode in ['ritz','spectral']:
    labels.append(lmp_mode)
    legend.append([lmp_mode+'-model',lmp_mode+'-control'])
    results_directory_test.append(results_dir_root+'res_n{}_no{}_ni{}_lmp-{}_sigmao{}_sigmab{}_Lb{}_reso-{}'.format(n,no,ni,lmp_mode,sigma_obs,sigmabvar,Lb,full_res))
    

diff_list=[]
obj_list=[]

for res_dir in results_directory_test:
    res1=np.genfromtxt(res_dir+'/lanczos_control_vs_PlanczosIF_model.dat', comments='#')
    res2=np.genfromtxt(res_dir+'/PlanczosIF_model_space.dat', comments='#')    
    res3=np.genfromtxt(res_dir+'/lanczos_control_space.dat', comments='#')
    diff_list.append(res1[:,3])
    obj_list.append([res2[:,3],res3[:,3]])
    

# maybe dirtyish but...
itot=list(range(len(obj_list[0][0])))

print(legend)

ylabel1=r'$J=J_o+J_b$'
ylabel2=r'$\Delta J$'
x=itot
xlabel='iterations'
io=outer_iterations
out_name=results_dir_root+"/test.png"

compare_plots_2N(out_name,obj_list,ylabel1,diff_list,ylabel2,x,xlabel,io,legend)

