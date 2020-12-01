#!/usr/bin/env python

################################################################################
#                        diff_methods_plot                                     #
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

################################################################################
# Path of the results file:
res_file = './diff_methods_results_dev/diff_methods_analysis_results/results.py'

# Name of the output dir:
out_dir = './diff_methods_results_dev/diff_methods_analysis_results/'

# burnin to apply (default = 0):
burnin = 0.9

# Load the raw results data:
res = pickle.load(open(res_file,"rb"))

# MCMC Chains parameters:
chain = np.array(res["chain"])
print ("shape of chain: ",np.shape(chain))

npars = len(chain[0][0])
nsteps = int(len(chain[0])*(1-burnin))
nwalkers = len(chain)

# MCMC parameters (walkers coordinates):
chain = chain[:,:nsteps,:].reshape((-1, npars))
print ("shape of flattened chain:",np.shape(chain))

# MCMC log-likelyhood:
lnprob = np.array(res["ln_prob"])
lnprob = lnprob[:,:nsteps].reshape(-1)
print ("shape of lnprob_data: ",np.shape(lnprob))

# Defines the quantiles to plot on the cornerplots:
quant=(0.02275,0.16, 0.84,0.97725)

################################################################################
labels=[r'$N_{obs}$',r'\sigma',r'L_b']
out_name=out_dir+"/cornerplot_test.png"
figure=make_corner_plot(chain,labels,quant,smoothing)
figure.savefig(out_name)
figure.clf()
