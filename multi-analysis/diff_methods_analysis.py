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
import emcee
import pickle
from multiprocessing import Pool

from diff_methods_functions import *
from diff_methods_lnlike import *

#---------------------------
# Verbose mode:
verb = True

# Analysis results paths:
out_dirs = []

# Get the current working directory and define other directories:
directories = {}
cwd = os.getcwd()
directories['analysis'] = cwd

# Path to the code multi:
os.chdir('..')
cwd = os.getcwd()
directories['multi'] = cwd

# Path to the executable from multi:
exec_command = os.path.join(directories['multi'] + '/build/bin/multi ')

# Compile the code (penser a le copier ds un autre repo):
os.chdir('./build')
os.system('ecbuild ..')
os.system('make')
os.chdir(directories['analysis'])

# Root of the results of the analysis:
results_dir_root = os.path.join(directories['analysis'] + '/diff_methods_results_dev/')
out_dirs.append(results_dir_root)

# Results directory:
diff_methods_dir = os.path.join(results_dir_root + '/diff_methods_analysis_results/')
out_dirs.append(diff_methods_dir)

# Location of the raw outputs of the code (temporary stored):
outputs_tmp_dir = os.path.join(directories['multi'] + '/outputs_tmp/')
directories["outputs"] = outputs_tmp_dir
out_dirs.append(outputs_tmp_dir)


# Location of the raw outputs of the code (temporary stored):
namelists_tmp_dir = os.path.join(directories['multi'] + '/namelists_tmp/')
directories["namelists"] = namelists_tmp_dir
out_dirs.append(namelists_tmp_dir)

# Creates the results directories:
# improvement: use os.mkdirs methods
                                 
for dir in out_dirs:
    if not os.path.exists(dir):
        os.mkdir(dir)

################################################################################
# Default values for the parameters and their range of variantion:
parameters = {}

# Improvement: use class for the parameters.

# Solver:
parameters['nm']                 = {'min':1 ,'max':3 ,'type':'int', 'val':3}
parameters['method']             = {'val':'"theoretical", "standard", "alternative"'}
parameters['na']                 = {'min':1 ,'max':2 ,'type':'int', 'val':2}
parameters['algorithm']               = {'val':'"lanczos", "planczosif"'}

parameters['no']                 = {'min':1 ,'max':10 ,'type':'int', 'val':3}
parameters['ni']                 = {'min':2 ,'max':10 ,'type':'int', 'val':6}

parameters['lmp_mode']           = {'val':'"none"'}
parameters['test_ortho']         = {'val':"F"}
parameters['shutoff_type']       = {'val':0}
parameters['shutoff_value']      = {'min':1e-4 ,'max':0.1 ,'type':'float' ,'val':1e-2}
parameters['transitive_interp']  = {'val':"T"}
parameters['projective_Bmatrix'] = {'val':"T"}

# Resolutions:
parameters['nx']                 = {'min':1 ,'max':1001 ,'type':'geom', 'val':'51,101,201'}
parameters['ny']                 = {'min':1 ,'max':1001 ,'type':'geom', 'val':'51,101,201'}

# Observations:
parameters['nobs']               = {'min':10 ,'max':1000 ,'type':'int', 'val':100}
parameters['sigma_obs']          = {'min':0. ,'max':1. ,'type':'float', 'val':0.01}

# Background:
parameters['sigmabvar']          = {'min':0. ,'max':1. ,'type':'float', 'val':0.}
parameters['Lb']                 = {'min':1e-15 ,'max':10. ,'type':'float', 'val':0.1}
parameters['spvarmin']           = {'min':0. ,'max':1. ,'type':'float', 'val':1e-5}

# Miscellanous:
parameters['new_seed']           = {'val':"F"}
parameters['filename']           = {'val':'"output.nc"'}

if verb:
    print('initial parameters: \n', parameters, '\n')

# Define the parameters space to sample:
parameters_to_sample = ['nobs', 'Lb']

#-------------------------------------------------------------------------------
# set the seed:
np.random.seed(42)

# Number of walkers, steps, dimensions and threads:
ndim = len(parameters_to_sample)
nwalkers = 80

nsteps = 100
nruns = int(nsteps/10.)
if nruns < 1:
    nruns = 1
    
# Scale factor of the stretch-move algorithm (see the doc: ...pdf)
scale_factor = 0.2

# Run the analysis:
p0 = walkers_create(nwalkers, parameters, parameters_to_sample, verb)
if verb:
    print('Generated walker example: \n', p0[0])
    print('shape of p0:', np.shape(p0),'\n')

# Test the behavior of the ln_prob:
#ln_prob(p0[0], parameters, parameters_to_sample, exec_command, verb)

#-------------------------------------------------------------------------------
for run in range(nruns):
    if not run == 0:
        p0 = state
    # Define the sampler object:
    ln_prob_args = (parameters, parameters_to_sample, exec_command,
                    directories,verb)
    
    print("---------- Starting run {} of the MCMC ---------- \n".format(run))
    
    # Parallelization of the code using pool:
    #with Pool() as pool:
    if True:
        #sampler = emcee.EnsembleSampler(nwalkers,ndim,ln_23,args=ln_prob_args,
 #                                       a=scale_factor,live_dangerously=True,
#                                        pool=pool)
        # Run the MCMC:
        #state = sampler.run_mcmc(p0, nsteps)
        
        sampler = emcee.EnsembleSampler(nwalkers,ndim,ln_23,args=ln_prob_args,
                                        a=scale_factor,live_dangerously=True)
        # Run the MCMC:
        state = sampler.run_mcmc(p0, nsteps)

        #-------------------------------------------------------------------------------
        # Save the results:
        print("save the results as pickles")
        # Save the chain containing the position (chain) and associated lnprobability (lnprob):
        results = {}
        results['chain'] = sampler.chain
        results['ln_prob'] = (-1)*sampler.lnprobability
        results_file = os.path.join(diff_methods_dir + f'results_run_{run}.py')
        print(results, results_file)
        with open(results_file, 'wb') as handle:
            pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print('results have been saved in:', results_file)
        
################################################################################
