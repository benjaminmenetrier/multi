#!/usr/bin/env python

################################################################################
# analysis_tools.py

# purpose: Contains usefull tools for the analysis.
# Author: Nicolas Baillot d'Etivaux

################################################################################

# Imported packages:
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec 
import matplotlib.ticker as ticker
from matplotlib import cm
from matplotlib.offsetbox import AnchoredText
import numpy as np
import netCDF4 as nc

# Allow the use of tex format for the labels:
plt.rcParams.update({"text.usetex": True, "font.size" : 15})

################################################################################
################################################################################
def netcdf_extract(res_dir):
    """ Extract the netcdf data from the result directory res_dir:
        Args: 
            res_dir: result directory in which is stored the netcdf output file.
        Return: 
            ds: netCDF4 dataset containing the results.
    """
    ds=nc.Dataset(res_dir+"/output.nc","r")
    return ds
################################################################################
################################################################################
# def concatenate_pickles(res_dir,start,end)
# # Results directory
# res_dir=sys.argv[1]

# # First run and last run to consider (burnin):
# start_chain=int(sys.argv[3])
# end_chain=int(sys.argv[4])

# all_res=np.load(res_dir+"results_{}.py".format(start_chain))

# all_chains=all_res['chain']
# all_ln_probs=all_res['ln_prob']

# nwalkers=len(all_chains)

# for run in range(start_chain,end_chain):
#     print("Concatenating the run ", run)
#     res = np.load(res_dir+"results_{}.py".format(run))
#     all_chains=np.append(all_chains,res['chain'],axis=1)
#     all_ln_probs=np.append(all_ln_probs,res['ln_prob'],axis=1)

# print "final chain shape = ", np.shape(all_chains)
# print "final ln_prob shape = ", np.shape(all_ln_probs)

# # Save the results in a dict:
# Res={}
# Res["lnprobs"]=all_lnprobs
# Res["chain"]=all_chains

# pickle.dump(Res,open(res_dir+"all_results.py","wb"),protocol=pickle.HIGHEST_PROTOCOL)
################################################################################
################################################################################
