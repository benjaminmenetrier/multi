#!/usr/bin/env python

################################################################################
# multi-analysis.py

# purpose: This code plots the results of the code multi.

# Author: Nicolas Baillot d'Etivaux

################################################################################
# Imported packages:
import os
import sys
import itertools
import multiprocessing
import concurrent.futures
# Personal packages:
from run_plots import *
from run_multi import *
from ensemble_analysis_plot import *
################################################################################
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

# Root directory of the results of the analysis:
results_dir_root = os.path.join(directories['analysis_results'] + '/test_randomize_mean')
if not os.path.exists(results_dir_root):
    os.mkdir(results_dir_root)

os.chdir(directories['run_analysis'])

################################################################################
# Parameters configurations:
i_no_ni = [[7,3],[3,2]]
i_nx = ['51,61,71,101']
i_nobs = [2000]
i_sigma_obs = [0.1]
i_Hnl_coeff = [0.0]
i_sigmabvar = [0.0]
i_Lb = [0.1]
i_interp_method = ['"spectral"','"bilinear"','"nearest"']
i_project_B = ["T"]
i_test_ortho = ["T"]

# Size of the ensemble tu run:
ensemble_size=5
i_rand_seed=list(range(ensemble_size))

# iterations to run:
iter_params = list(itertools.product(i_no_ni, i_nx, i_nobs, i_sigma_obs,
                                i_Hnl_coeff, i_sigmabvar, i_Lb,
                                i_interp_method, i_project_B,
                                i_test_ortho,repeat=1))

# Creates the results directories and outer_iterations_list:
res_dir_dict = {}
outer_iterations_dict = {}
for rand_seed in i_rand_seed:
    res_dir_list, outer_iterations_list = create_res_dirs(rand_seed, iter_params, results_dir_root)
    res_dir_dict[rand_seed] = res_dir_list
    outer_iterations_dict[rand_seed] = outer_iterations_list

print('Running the code multi')
################################################################################
#-------------------------------------------------------------------------------
# Analysis for each seed independently:
def run_element_analysis(rand_seed, iter_params, directories, res_dir_list, outer_iterations_list):
    """Runs the code multi and then produces the individual plots for every elements
    of the ensemble.
    """
    # If True, gives more plots of all the variables:
    extra_monitoring=False
    
    run_multi(rand_seed, iter_params, directories, res_dir_list)
    #run_plots(res_dir_list, outer_iterations_list, extra_monitoring)
#-------------------------------------------------------------------------------

print('Starting analysis:')

processes = []
with concurrent.futures.ProcessPoolExecutor() as executor:
    for rand_seed in i_rand_seed:
        res_dir_list = res_dir_dict[rand_seed]
        outer_iterations_list = outer_iterations_dict[rand_seed]
        args = (rand_seed, iter_params, directories, res_dir_list, outer_iterations_list)
        processes.append(executor.submit(run_element_analysis, *args))
        
    #for process in concurrent.futures.as_completed(processes):
    #    ensemble_results.append(process.result())
        
################################################################################
print('Starting ensemble analysis')

ensemble_compare_methods_plot(res_dir_dict, outer_iterations_dict)

print('analysis completed')   
################################################################################
