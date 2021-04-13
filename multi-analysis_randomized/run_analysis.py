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
from multiprocessing import Pool
# Personal packages:
from run_plots import *
from run_multi import *
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
results_dir_root = os.path.join(directories['analysis_results'] + '/test_randomize_T')
if not os.path.exists(results_dir_root):
    os.mkdir(results_dir_root)

os.chdir(directories['run_analysis'])

################################################################################
# Parameters configurations:
i_no_ni = [[4,6]]
i_nx = ['51,61,71,101']
i_nobs = [2000]
i_sigma_obs = [0.1]
i_Hnl_coeff = [0.0]
i_sigmabvar = [0.0]
i_Lb = [0.1]
#i_interp_method = ['"spectral"','"bilinear"','"nearest"']
i_interp_method = ['"spectral"']
i_project_B = ["T"]
i_test_ortho = ["T"]

# Size of the ensemble tu run:
ensemble_size=2
i_rand_seed=list(range(ensemble_size))

# If True, gives more plots of all the variables:
extra_monitoring=False

# iterations to run:
iter_params = list(itertools.product(i_no_ni, i_nx, i_nobs, i_sigma_obs,
                                i_Hnl_coeff, i_sigmabvar, i_Lb,
                                i_interp_method, i_project_B,
                                i_test_ortho,repeat=1))

print('Running the code multi')
################################################################################
for rand_seed in i_rand_seed:
   res_dir_list, outer_iterations_list = run_multi(rand_seed, iter_params, directories, results_dir_root)

#with Pool(4) as p:
#    p.map(run_multi, [(12, iter_params, directories, results_dir_root),(12, iter_params, directories, results_dir_root)])

print('Starting analysis:')
################################################################################

#run_plots(res_dir_list, outer_iterations_list, extra_monitoring)

print('analysis completed')   
################################################################################
