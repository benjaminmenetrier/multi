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
def namelist_write(p,directory):

    namelist = open(os.path.join(directory + "/namelist"),"w")

    # Solver:
    namelist.write("&solver\n")
    namelist.write("no = {}\n".format(p[0]))
    namelist.write("ni = {}\n".format(p[1]))
    namelist.write("lmp_mode = {}\n".format(p[2]))
    namelist.write("test_ortho = {}\n".format(p[3]))
    namelist.write("shutoff_type = {}\n".format(p[4]))
    namelist.write("shutoff_value = {}\n".format(p[5]))
    namelist.write("method = {}\n".format(p[6]))
    namelist.write("transitive_interp = {}\n".format(p[7]))
    namelist.write("projective_Bmatrix = {}\n".format(p[8]))
    namelist.write("/\n\n")
    # Resolutions:
    namelist.write("&resolutions\n")
    namelist.write("nx = {}\n".format(p[9]))
    namelist.write("ny = {}\n".format(p[10]))
    namelist.write("/\n\n")
    # Observations:
    namelist.write("&observations\n")
    namelist.write("nobs = {}\n".format(p[11]))
    namelist.write("sigma_obs = {}\n".format(p[12]))
    namelist.write("/\n\n")
    # Background:
    namelist.write("&background\n")
    namelist.write("sigmabvar = {}\n".format(p[13]))
    namelist.write("Lb = {}\n".format(p[14]))
    namelist.write("spvarmin = {}\n".format(p[15]))
    namelist.write("/\n\n")
    # Miscellanous:
    namelist.write("&miscellanous\n")
    namelist.write("new_seed = {}\n".format(p[16]))
    namelist.write("filename = {}\n".format(p[17]))
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
