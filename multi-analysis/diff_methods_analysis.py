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
import emcee
from diff_methods_functions import *
import pickle
from multiprocessing import Pool

#---------------------------
# if True, the results already existing for the same set of parameters
# will be rewritten:
rewrite_results=True

# Verbose mode:
verb=True
#---------------------------

# Path to the executable from multi:
exec_command='./build/bin/multi '

# Analysis results paths:
out_dirs=[]

# Get the current working directory and define other directories:
directories={}
cwd=os.getcwd()
directories["analysis"]=cwd

os.chdir('..')
cwd=os.getcwd()
directories["code"]=cwd

# Compile the code (penser a le copier ds un autre repo):
os.chdir('./build')
os.system('ecbuild ..')
os.system('make')
os.chdir('..')

# cwd=os.getcwd()
# directories["build"]=cwd

# Root of the results of the analysis:
results_dir_root=directories["analysis"]+'/diff_methods_results_dev/'
out_dirs.append(results_dir_root)

# Results directory:
diff_methods_dir=results_dir_root+'/diff_methods_analysis_results/'
out_dirs.append(diff_methods_dir)

# Location of the raw outputs of the code (temporary stored):
outputs_tmp_dir=directories["code"]+'/outputs_tmp/'
directories["outputs"]=outputs_tmp_dir
out_dirs.append(outputs_tmp_dir)


# Location of the raw outputs of the code (temporary stored):
namelists_tmp_dir=directories["code"]+'/namelists_tmp/'
directories["namelists"]=namelists_tmp_dir
out_dirs.append(namelists_tmp_dir)

# Creates the results directories:
for dir in out_dirs:
    if not os.path.exists(dir):
        os.mkdir(dir)

################################################################################
# Analysis design:

#-------------------------------------------------------------------------------
# Default values for the parameters and their range of variantion:
parameters={}

# Solver:
parameters['no']={'min':1 ,'max':10 ,'type':'int', 'val':4}
parameters['ni']={'min':2 ,'max':10 ,'type':'int', 'val':6}
parameters['lmp_mode']={'val':'"none"'} # lmps: none, ritz or spectral
parameters['test_ortho']={'val':".false."}
parameters['shutoff_type']={'val':0}
parameters['shutoff_value']={'min':1e-4 ,'max':0.1 ,'type':'float' ,'val':1e-2}
parameters['method']={'val':'"theoretical"'} # methods: theoretical, standard or alternative
parameters['transitive_interp']={'val':".true."}
parameters['projective_Bmatrix']={'val':".true."}
# Resolutions:
parameters['nx']={'min':1 ,'max':1001 ,'type':'geom', 'val':'11,31,51,101'}
parameters['ny']={'min':1 ,'max':1001 ,'type':'geom', 'val':'11,31,51,101'}
# Observations:
parameters['nobs']={'min':10 ,'max':1000 ,'type':'int', 'val':100}
parameters['sigma_obs']={'min':0. ,'max':1. ,'type':'float', 'val':0.01}
# Background:
parameters['sigmabvar']={'min':0. ,'max':1. ,'type':'float', 'val':0.1}
parameters['Lb']={'min':1e-15 ,'max':10. ,'type':'float', 'val':0.1}
parameters['spvarmin']={'min':0. ,'max':1. ,'type':'float', 'val':1e-5}
# Miscellanous:
parameters['new_seed']={'val':".false."}
parameters['filename']={'val':'"output.nc"'}

if verb:
    print('initial parameters: \n', parameters, '\n')

# Define the parameters space to sample:
parameters_to_sample=['nobs', 'sigma_obs','Lb']
#parameters_to_sample=['nx','ny']

# Define the methods to compare:
methods_list=['"standard"','"alternative"']

#-------------------------------------------------------------------------------
# set the seed:
np.random.seed(42)

# Number of walkers, steps, dimensions and threads:
ndim=len(parameters_to_sample)
nwalkers = 6

nsteps = 1
nruns = int(nsteps/10.)
if nruns <1:
    nruns=1
    
# Scale factor of the stretch-move algorithm (see the doc: ...pdf)
scale_factor = 0.2


# Run the analysis:
p0 = walkers_create(nwalkers,parameters,parameters_to_sample,verb)
if verb:
    print('Generated walker example: \n', p0[0])
    print('shape of p0:',np.shape(p0),'\n')
# Test the behavior of the ln_prob:
#ln_prob(p0[0],parameters,parameters_to_sample,methods_list,exec_command,verb)

#-------------------------------------------------------------------------------
for run in range(nruns):
    if not run==0:
        p0 = state

        
    print("---------- Starting the MCMC ---------- \n")
    # Parallelization of the code using pool:
    with Pool() as pool:
        # Define the sampler object:
        ln_prob_args=(parameters,parameters_to_sample,methods_list,exec_command,directories,verb)
        sampler = emcee.EnsembleSampler(nwalkers,ndim,ln_prob,args=ln_prob_args,
                                        a=scale_factor,live_dangerously=True,pool=pool)
        # Run the MCMC:
        state = sampler.run_mcmc(p0, nsteps)
        
        # sampler = emcee.EnsembleSampler(nwalkers,ndim,ln_prob,args=ln_prob_args,
        #                                 a=scale_factor,live_dangerously=True)
        # # Run the MCMC:
        # state = sampler.run_mcmc(p0, nsteps)


        #-------------------------------------------------------------------------------
        # Save the results:
        print("save the results as pickles")
        # Save the chain containing the position (chain) and associated lnprobability (lnprob):
        results={}
        results['chain'] = sampler.chain
        results['ln_prob'] = (-1)*sampler.lnprobability
        pickle.dump(results,open(diff_methods_dir+"results_run_{}.py".format(run),"wb"))
        print('results have been saved in:',diff_methods_dir+'results_run_{}.py'.format(run))
################################################################################
