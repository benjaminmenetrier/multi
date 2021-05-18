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
import time
# Personal packages:
from run_plots import *
from run_multi import *
from ensemble_analysis_plot import *
from analysis_tools import init_state_gen
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
# /!\ the word "seed" is forbidden in this string because of the use of a split
# method in the function that builds the results files.

#results_dir_root = os.path.join(directories['analysis_results'] + '/')
results_dir_root = '/media/bayow/HDD2/multi-results/read_test_w_obs_xb_and_xt'

if not os.path.exists(results_dir_root):
    os.mkdir(results_dir_root)

os.chdir(directories['run_analysis'])

start = time.perf_counter()
################################################################################
# Parameters configurations:

i_no_ni = [[4,6]]
i_nx = ['31,41,61,101']
i_nobs = [2000]
i_sigma_obs = [0.1]
i_Hnl_coeff = [0.]
i_sigmabvar = [0.0]
i_Lb = [0.1]
i_interp_method = ['"spectral"','"bilinear"','"nearest"']
i_project_B = ["T"]
i_test_ortho = ["T"]

# Size of the ensemble tu run:
ensemble_size=1
i_rand_seed=list(range(ensemble_size))

# Determine the maximum full resolution:
# (random vectors will be truncated in the src/main code when the resolution is too high)
full_res = 0
for ii_nx in i_nx:
    if int(ii_nx.split(',')[-1]) > full_res:
        full_res = int(ii_nx.split(',')[-1])
print(full_res)

# Determine the maximum number of observations:
nobs_max = max(i_nobs)

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

# Creates the initial state (background, truth and observations):
init_state_gen(full_res, nobs_max, directories['multi'])

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

results_obj = build_results_object(res_dir_dict, outer_iterations_dict)

ensemble_compare_methods_plot(res_dir_dict, outer_iterations_dict, results_obj)

#linearization_check(res_dir_dict, outer_iterations_dict, results_obj)

stop = time.perf_counter()

print(f'analysis completed in {stop - start} seconds')   
################################################################################
