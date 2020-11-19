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

# Default values for the parameters in namelist:
#---------------------------
no=2
ni=4

lmp_mode='"none"'
test_ortho=".false."
shutoff_type=0
shutoff_value=0.0
method='"theoretical"'
transitive_interp=".true."
projective_Bmatrix=".true."

nx="101,121,151"
ny="101,121,151"

nobs=3000
sigma_obs=0.1

sigmabvar=0.0
Lb=0.12
spvarmin=1.0e-5

new_seed=".false."
filename='"output"'
#---------------------------


#--------------------------------------------------------------------------------
# if True, the results already existing will be rewritten:
rewrite_results=True

# Path to the executable from multi:
exec_command='./build/bin/multi '

# Location of the raw outputs of the code:
code_output='../output.nc'
#--------------------------------------------------------------------------------

# Analysis results paths:
out_dirs=[]

# Root of the results of the analysis:
results_dir_root='./analysis_results_dev/'
out_dirs.append(results_dir_root)

# Raw results of the analysis: 
res_dir_raw=results_dir_root+'raw_results/'
out_dirs.append(res_dir_raw)

# Results for the comparision between LMP modes:
compare_methods_dir=results_dir_root+'compare_methods/'
out_dirs.append(compare_methods_dir)

# # Results for the check of second-level lmp:
# res_dir_lmp_check=results_dir_root+'check_second_level_lmp/'
# out_dirs.append(res_dir_lmp_check)

# # Results for rho:
# res_dir_lmp_compare_rho=results_dir_root+'lmp_compare_rho/'
# out_dirs.append(res_dir_lmp_compare_rho)

# # Results for beta:
# res_dir_lmp_compare_beta=results_dir_root+'lmp_compare_beta/'
# out_dirs.append(res_dir_lmp_compare_beta)

# Results for the outer vectors:
# res_dir_outer_vectors=results_dir_root+'outer_vectors/'
# out_dirs.append(res_dir_outer_vectors)

# Creates the results directories:
for dir in out_dirs:
    if not os.path.exists(dir):
        os.mkdir(dir)
    
# Get the current working directory:
cwd=os.getcwd()

################################################################################
# Run the code, and store the results:

# Store the results directories produced in the following loop:
res_dir_list=[]
# Outer iteraions for plotting:
outer_iterations_list=[]

# Loop over the parameters and run the code:
for no in [4]:
    for ni in [6]:
        for lmp_mode in ['"none"']:
            for method in ['"theoretical"','"standard"','"alternative"']:
                for nx in ['101']:
                    for nobs in [2000]:
                        for sigma_obs in [0.01]:
                            for sigmabvar in [0.1]:
                                for Lb in [0.12]:
                                    # square grid:
                                    ny=nx
                                    #ny='1'
                                    # Outer iteraions for plotting:
                                    outer_iterations=[]
                                    for io in range(no+1):
                                        outer_iterations.append((ni)*io)
                                        outer_iterations_list.append(outer_iterations)
                                        
                                    # Create the results directory:
                                    name_string1='res_no{}_ni{}_lmp-{}_met-{}'
                                    name_string1=name_string1.format(no,ni,lmp_mode.replace('"',''),method.replace('"',''))
                                    name_string2='_nx{}_ny{}_n-obs{}_sigmaobs{}'
                                    name_string2=name_string2.format(nx,ny,nobs,sigma_obs)
                                    name_string3='_sigbvar{}_Lb{}'
                                    name_string3=name_string3.format(sigmabvar,Lb)
                                    
                                    res_dir=res_dir_raw+name_string1+name_string2+name_string3
                                    res_dir_list.append(res_dir)
                                    # Rewrite the results or not:
                                    if not rewrite_results and os.path.exists(res_dir):
                                        continue
                                    if not os.path.exists(res_dir):
                                        os.mkdir(res_dir)
                                        
                                    # Write the parameters in namelist file:
                                    namelist=open("../namelist","w")
                                    # Solver:
                                    namelist.write("&solver\n")
                                    namelist.write("no = {}\n".format(no))
                                    namelist.write("ni = {}\n".format(ni))
                                    namelist.write("lmp_mode = {}\n".format(lmp_mode))
                                    namelist.write("test_ortho = {}\n".format(test_ortho))
                                    namelist.write("shutoff_type = {}\n".format(shutoff_type))
                                    namelist.write("shutoff_value = {}\n".format(shutoff_value))
                                    namelist.write("method = {}\n".format(method))
                                    namelist.write("transitive_interp = {}\n".format(transitive_interp))
                                    namelist.write("projective_Bmatrix = {}\n".format(projective_Bmatrix))
                                    namelist.write("/\n\n")
                                    # Resolutions:
                                    namelist.write("&resolutions\n")
                                    namelist.write("nx = {}\n".format(nx))
                                    namelist.write("ny = {}\n".format(ny))
                                    namelist.write("/\n\n")
                                    # Observations:
                                    namelist.write("&observations\n")
                                    namelist.write("nobs = {}\n".format(nobs))
                                    namelist.write("sigma_obs = {}\n".format(sigma_obs))
                                    namelist.write("/\n\n")
                                    # Background:
                                    namelist.write("&background\n")
                                    namelist.write("sigmabvar = {}\n".format(sigmabvar))
                                    namelist.write("Lb = {}\n".format(Lb))
                                    namelist.write("spvarmin = {}\n".format(spvarmin))
                                    namelist.write("/\n\n")
                                    # Miscellanous:
                                    namelist.write("&miscellanous\n")
                                    namelist.write("new_seed = {}\n".format(new_seed))
                                    namelist.write("filename = {}\n".format(filename))
                                    namelist.write("/\n")
                                    namelist.close()
                                    
                                    # Run the code and save the results:
                                    # Compile the code:
                                    os.chdir("../build")
                                    os.system("ecbuild ..")
                                    os.system("make")
                                    os.chdir(cwd)
                                    os.chdir('..')
                                    os.system(exec_command)
                                    os.chdir(cwd)
                                    copyfile(code_output,res_dir+"/output.nc")
                                    
################################################################################
# Plot comparision between lanczos and Planczos for cost function:
# for r,res_dir in enumerate(res_dir_list):

#     # Plots the observations
#     obs_plot(res_dir)
#     hxg_plot(res_dir)
#     innovation_plot(res_dir)

#     # Plots the comparision between lanczos and planczosif for the cost function:
#     lanczos_vs_planczosif_plot(res_dir,outer_iterations_list[r])

#     # Plots the model fields as matrices and scatterplots (what is the best?):
#     ds=netcdf_extract(res_dir)
#     for io in ds.groups:
#         x_coord=np.array(ds[io]['x_coord'])
#         y_coord=np.array(ds[io]['y_coord'])
#         for field in ['sigmab','dirac_cov','dirac_cor','xb']:
#             matrix=np.array(ds[io][field][:])
#             out_name=res_dir+'/'+field+'_'+io
#             field_plot(matrix,out_name)
#         for algo in ds[io].groups:
#             for field in ['xg']:
#                 matrix=np.array(ds[io][algo][field][:])
#                 out_name=res_dir+'/'+algo+'_'+field+'_'+io
#                 field_plot(matrix,out_name)
#                 # Plots the increment:
#                 for ii in range(ds[io][algo].dimensions['nimax'].size):
#                     dx=np.array(ds[io][algo]['dx'][ii][:])
#                     out_name=res_dir+'/'+algo+'_dx_'+io+'_'+'inner_'+str(ii)
#                     field_plot(dx,out_name)
#--------------------------------------------------------------------------------    
# Comparision between the methods:
methods_list=['theoretical','standard','alternative']
for r,res_dir in enumerate(res_dir_list):
    try:    
        if methods_list[0] in res_dir:
            res_tmp=res_dir.split('theoretical')
            res_tmp1,res_tmp2=res_tmp[0],res_tmp[1]
            # Store the output files names
            compare_methods_out=res_tmp1+'compare'+res_tmp2
            compare_methods_out=compare_methods_out.replace(res_dir_raw,compare_methods_dir)
            if not os.path.exists(compare_methods_out):
                os.mkdir(compare_methods_out)
            compare_methods=[]
            for method in methods_list:
                compare_methods.append(res_tmp1+method+res_tmp2)
            compare_methods_data=[]
            for res in compare_methods:
                ds=netcdf_extract(res)
                compare_methods_data.append(ds)
            compare_methods_plot(compare_methods_data,methods_list,outer_iterations_list[r],compare_methods_out)    
    except:
        print("Cannot compare methods: the following file does not exist:\n",res)
################################################################################




# # ##############################################################################
# # Plots the outer vectors for each outer loop:

# # Plots the innovation:
# for res_dir in res_dir_list:
#    yo_vs_hxg_plot(res_dir)

# res_file_names=['lanczos_control_space_outer_vectors.dat','PlanczosIF_model_space_outer_vectors.dat']
# for f,res_file_name in enumerate(res_file_names):
#     # Plots the guess:
#     out_name='_guess.png'
#     label=r'$x^g'
#     column_of_interest=-4
#     for res_dir in res_dir_list:
#         results_file=res_dir+res_file_name
#         out_file_name=results_file[:-4]+out_name
#         print('out_file',out_file_name)
#         vec_plot(results_file,column_of_interest,label,out_file_name)

#     # Plots the background:
#     out_name='_background.png'
#     label=r'$x^b'
#     column_of_interest=5
#     for res_dir in res_dir_list:
#         results_file=res_dir+'/'+res_file_name
#         out_file_name=results_file[:-4]+out_name
#         vec_plot(results_file,column_of_interest,label,out_file_name)

#     # Plots the increment:
#     out_name='_increment.png'
#     label=r'$\delta x^b'
#     column_of_interest=4
#     for res_dir in res_dir_list:
#         results_file=res_dir+'/'+res_file_name
#         out_file_name=results_file[:-4]+out_name
#         vec_plot(results_file,column_of_interest,label,out_file_name)


#     # Plots the obs:
#     out_name='_obs.png'
#     label=r'$y^o'
#     column_of_interest=-2
#     for res_dir in res_dir_list:
#         results_file=res_dir+res_file_name
#         out_file_name=results_file[:-4]+out_name
#         vec_plot(results_file,column_of_interest,label,out_file_name)


#     # Plots Hxg:
#     out_name='_hxg.png'
#     label=r'$H x^g'
#     column_of_interest=-3
#     for res_dir in res_dir_list:
#         results_file=res_dir+'/'+res_file_name
#         out_file_name=results_file[:-4]+out_name
#         vec_plot(results_file,column_of_interest,label,out_file_name)


#     # Plots the innovation:
#     out_name='_innovation.png'
#     #label=r'$y^o_{}-H x^g'
#     label=r'$d'
#     column_of_interest=-1
#     for res_dir in res_dir_list:
#         results_file=res_dir+'/'+res_file_name
#         out_file_name=results_file[:-4]+out_name
#         vec_plot(results_file,column_of_interest,label,out_file_name)

#     # Plots the preconditionned vectors:
#     # Lanczos:
#     column_of_interest=2
#     if f==0:
#         out_name='_dva_interp.png'
#         label=r'$\Pi \delta v^{a}'
#     elif f==1:
#         out_name='_dxabar_interp.png'
#         label=r'$\Pi \delta \bar{x^{a}}'
#     for res_dir in res_dir_list:
#         results_file=res_dir+'/'+res_file_name
#         out_file_name=results_file[:-4]+out_name
#         vec_plot(results_file,column_of_interest,label,out_file_name)

#     # PlanczosIF:
#     column_of_interest=3
#     if f==0:
#         out_name='_dvb.png'
#         label=r'$\delta v^b'
#     elif f==1:
#         out_name='_dxbbar.png'
#         label=r'$\delta \bar{x^b}'
#     for res_dir in res_dir_list:
#         results_file=res_dir+'/'+res_file_name
#         out_file_name=results_file[:-4]+out_name
#         vec_plot(results_file,column_of_interest,label,out_file_name)


# ################################################################################
# # Comparision of LMP methods:

# # Output filenames:
# out_names_J=[]
# out_names_check=[]
# out_names_rho=[]
# out_names_beta=[]
# # Store the outer_itertaions:
# outer_iterations_list_tmp=[]
# # Store the results files to compare:
# lmp_to_compare=[]
# check_second_level_lmp_dirs=[]

# for r,res_dir in enumerate(res_dir_list):
#     if 'ritz' in res_dir:
#         res_tmp=res_dir.split('ritz')
#         res_tmp1,res_tmp2=res_tmp[0],res_tmp[1]
#         # Store the output files names
#         out_name=res_tmp1+'compare'+res_tmp2
#         out_name=out_name.split(res_dir_raw)[1]

#         out_name_J=res_dir_lmp_compare_J+out_name[:-1]+'.png'
#         out_name_check=res_dir_lmp_check+out_name[:-1]+'.png'
#         out_name_rho=res_dir_lmp_compare_rho+out_name[:-1]+'.png'
#         out_name_beta=res_dir_lmp_compare_beta+out_name[:-1]+'.png'
#         out_names_J.append(out_name_J)
#         out_names_check.append(out_name_check)
#         out_names_rho.append(out_name_rho)
#         out_names_beta.append(out_name_beta)
        
#         # Store the outer iterations:
#         outer_iterations_list_tmp.append(outer_iterations_list[r])

#         # Store the results files to compare:
#         lmp_to_compare_tmp=[]
#         check_second_level_lmp_tmp=[]
#         for lmp in ['ritz','spectral','none']:
#             lmp_to_compare_tmp.append(res_tmp1+lmp+res_tmp2)
#             if not lmp=='none':
#                 check_second_level_lmp_tmp.append(res_tmp1+lmp+res_tmp2)
#         lmp_to_compare.append(lmp_to_compare_tmp)
#         check_second_level_lmp_dirs.append(check_second_level_lmp_tmp)
        
# # Plots the comparision of LMP methods according to J:
# column_of_interest=3
# ylabel1=r'$J=J_o+J_b$'
# ylabel2=r'$J_{B^{1/2}}-J_{B}$'
# lmp_compare(out_names_J,lmp_to_compare,column_of_interest,ylabel1,ylabel2,outer_iterations_list)
# # Plots the comparision of LMP methods according to rho:
# column_of_interest=6
# ylabel1=r'$\rho$'
# ylabel2=r'$\rho_{B^{1/2}}-\rho_{B}$'
# lmp_compare(out_names_rho,lmp_to_compare,column_of_interest,ylabel1,ylabel2,outer_iterations_list)
# # Plots the comparision of LMP methods according to the beta:
# column_of_interest=7
# ylabel1=r'$\beta$'
# ylabel2=r'$\beta_{B^{1/2}}-\rho_{B}$'
# lmp_compare(out_names_beta,lmp_to_compare,column_of_interest,ylabel1,ylabel2,outer_iterations_list)
# ################################################################################

# # Check that the difference between the second level preconditionners:
# #check_second_level_lmp(out_names_check,check_second_level_lmp_dirs,outer_iterations_list)




# See later the diff_plots:
################################################################################
# Plot the differences as a function of the resolution, see later.
# # Output filenames:
# out_names=[]
# # Store the outer_itertaions:
# outer_iterations_list_tmp=[]
# # Store the results files to compare:
# diff_vs_nres=[]

# for r,res_dir in enumerate(res_dir_list):
#         res_tmp=res_dir.split('nres')
#         res_tmp1,res_tmp2=res_tmp[0],res_tmp[1]
#         # Store the output files names
#         out_name=res_tmp1+'compare'+res_tmp2
#         out_name=out_name.split(res_dir_raw)[1]
#         out_name=res_dir_diff_vs_nres+out_name[:-1]+'.png'
#         out_names.append(out_name)
#         # Store the outer iterations:
#         outer_iterations_list_tmp.append(outer_iterations_list[r])
#         # Store the results files to compare:
#         diff_vs_nres_tmp=[]
#         for nres in [128,2048]:
#             diff_vs_nres_tmp.append(res_tmp1+str(nres)+res_tmp2)
#         diff_vs_nres.append(diff_vs_nres_tmp)
#         #print(diff_vs_nres)
# # Plots the comparision of LMP methods:
# #diff_plot(out_names,diff_vs_nres,outer_iterations_list)        
# ################################################################################


