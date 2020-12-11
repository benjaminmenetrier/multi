#!/usr/bin/env python

################################################################################
# diff_methods_plot

# purpose: Contains functions plotting comparision between the methods
# Author: Nicolas Baillot d'Etivaux

################################################################################

# imported pachages:
import numpy as np
import pickle
import corner
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from corner_plot import *

# plot style:
smoothing = 0.9

plt.rc('font', family='serif')
plt.rc('font', size=10)          # controls default text sizes
plt.rc('axes', labelsize=13)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels

# Allow the use of tex format for the labels:
plt.rcParams.update({"text.usetex": True, "font.size" : 15})

################################################################################
# Path of the results file:
res_file = './diff_methods_results_dev/diff_methods_analysis_results/all_results.py'

# Name of the output dir:
out_dir = './diff_methods_results_dev/diff_methods_analysis_results/'

# burnin to apply (default = 0):
burnin = 0.9

# Load the raw results data:
res = pickle.load(open(res_file,"rb"))

# MCMC Chains parameters:
chain = np.array(res["chain"])
print ("Shape of chain: ",np.shape(chain))

npars = len(chain[0][0])
nsteps = int(len(chain[0])*(1-burnin))
nwalkers = len(chain)

# MCMC parameters (walkers coordinates):
flat_chain = chain[:, :nsteps, :].reshape((-1, npars))
print ("Shape of flattened chain:",np.shape(flat_chain))

# MCMC log-likelyhood:
lnprob = np.array(res["ln_prob"])
print("Shape of ln_prob: ", np.shape(lnprob))
flat_lnprob = lnprob[:, :nsteps].reshape(-1)
print ("Shape of flattened lnprob: ",np.shape(flat_lnprob))

# Defines the quantiles to plot on the cornerplots:
quant = (0.02275, 0.16, 0.84, 0.97725)

################################################################################
# Corner plot for the chain:
print("plotting cornerplot for the chain:")

flat_lnprob = np.atleast_2d(flat_lnprob)

chain_vs_lnprob = np.append(flat_lnprob.T, flat_chain, axis = 1)
print('shape of chain_vs_lnprob', np.shape(chain_vs_lnprob))

labels = [ r'$diff$', r'$N_{obs}$', r'$L_b$']

out_name = out_dir + "/cornerplot_test.png"
figure = make_corner_plot(chain_vs_lnprob, labels, quant, smoothing)
figure.savefig(out_name)
figure.clf()
################################################################################
# Convergence of the MCMC:
# Raw plot with all the walkers:
out_name = out_dir + '/mean_difference.png'
ylabel = 'mean difference'
xlabel = 'steps'

print(lnprob[0][0])
print(chain[0][0][:])
mean_lnprob = []
for step in lnprob.T:
    mean_lnprob.append(sum(step) / (nwalkers * 18.))

fig = plt.figure()
plt.plot(range(len(mean_lnprob)), mean_lnprob[:])
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.savefig(out_name)
plt.close()
################################################################################
