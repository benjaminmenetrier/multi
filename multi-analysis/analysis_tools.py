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
import os
import netCDF4 as nc
import itertools

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
# Write the parameters in namelist file:
def namelist_write(parameters,directory):
    """Write the namelist file according to the given parameters.
    """
    
    parameters_name={}
    parameters_name['solver'] = ['nm', 'method', 'na', 'algorithm', 'no', 'ni', 'lmp_mode', 'test_ortho',
                                 'shutoff_type', 'shutoff_value',
                                 'interp_method', 'projective_Bmatrix',]
    parameters_name['resolution'] = ['nx', 'ny']
    parameters_name['obs'] = ['nobs', 'sigma_obs']
    parameters_name['background'] = ['sigmabvar', 'Lb', 'spvarmin']
    parameters_name['miscellanous'] = ['new_seed', 'filename']

    namelist = open(os.path.join(directory + "/namelist"),"w")

    # Solver:
    namelist.write("&solver\n")
    for p,par in enumerate(parameters['solver']):
        namelist.write(parameters_name['solver'][p] + " = {}\n".format(par))
        print(parameters_name['solver'][p] + " = {}\n".format(par))
    namelist.write("/\n\n")
    # Resolutions:
    namelist.write("&resolutions\n")
    for p,par in enumerate(parameters['resolution']):
        namelist.write(parameters_name['resolution'][p] + " = {}\n".format(par))
    namelist.write("/\n\n")
    # Observations:
    namelist.write("&observations\n")
    for p,par in enumerate(parameters['obs']):
        namelist.write(parameters_name['obs'][p] + " = {}\n".format(par))
    namelist.write("/\n\n")
    # Background:
    namelist.write("&background\n")
    for p,par in enumerate(parameters['background']):
        namelist.write(parameters_name['background'][p] + " = {}\n".format(par))
    namelist.write("/\n\n")
    # Miscellanous:
    namelist.write("&miscellanous\n")
    for p,par in enumerate(parameters['miscellanous']):
        namelist.write(parameters_name['miscellanous'][p] + " = {}\n".format(par))
    namelist.write("/\n")
    
    namelist.close()


    # # Solver:
    # namelist.write("&solver\n")
    # namelist.write("no = {}\n".format(no))
    # namelist.write("ni = {}\n".format(ni))
    # namelist.write("lmp_mode = {}\n".format(lmp_mode))
    # namelist.write("test_ortho = {}\n".format(test_ortho))
    # namelist.write("shutoff_type = {}\n".format(shutoff_type))
    # namelist.write("shutoff_value = {}\n".format(shutoff_value))
    # namelist.write("method = {}\n".format(method))
    # namelist.write("transitive_interp = {}\n".format(transitive_interp))
    # namelist.write("projective_Bmatrix = {}\n".format(projective_Bmatrix))
    # namelist.write("/\n\n")
    # # Resolutions:
    # namelist.write("&resolutions\n")
    # namelist.write("nx = {}\n".format(nx))
    # namelist.write("ny = {}\n".format(ny))
    # namelist.write("/\n\n")
    # # Observations:
    # namelist.write("&observations\n")
    # namelist.write("nobs = {}\n".format(nobs))
    # namelist.write("sigma_obs = {}\n".format(sigma_obs))
    # namelist.write("/\n\n")
    # # Background:
    # namelist.write("&background\n")
    # namelist.write("sigmabvar = {}\n".format(sigmabvar))
    # namelist.write("Lb = {}\n".format(Lb))
    # namelist.write("spvarmin = {}\n".format(spvarmin))
    # namelist.write("/\n\n")
    # # Miscellanous:
    # namelist.write("&miscellanous\n")
    # namelist.write("new_seed = {}\n".format(new_seed))
    # namelist.write("filename = {}\n".format(filename))
    # namelist.write("/\n")
    # namelist.close()

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
