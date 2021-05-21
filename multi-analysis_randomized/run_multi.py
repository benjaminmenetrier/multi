#!/usr/bin/env python

################################################################################
# run_multi.py

# purpose: Contains the function that runs the code multi.

# Author: Nicolas Baillot d'Etivaux

################################################################################
from distutils.dir_util import copy_tree
from shutil import copyfile, copytree, rmtree
import os
import numpy as np
import concurrent.futures
# Personal packages:
from analysis_tools import *
################################################################################

################################################################################
################################################################################
def run_multi(nm, method, na, algo, no, ni, lmp_mode, test_ortho, shutoff_type,
              shutoff_value, interp_method, projective_Bmatrix,
              nx, ny, nobs, sigma_obs,
              Hnl_coeff, sigmabvar, Lb, spvarmin, new_seed, filename,
              rand_seed, directories, res_dir):
   
    parameters={}
    parameters['solver'] = [nm, method, na, algo, no, ni, lmp_mode, test_ortho,
                            shutoff_type, shutoff_value,
                            interp_method, projective_Bmatrix]
    parameters['resolution'] = [nx, ny]
    parameters['obs'] = [nobs, sigma_obs, Hnl_coeff]
    parameters['background'] = [sigmabvar, Lb, spvarmin]
    parameters['miscellanous'] = [new_seed, filename]

    namelist_write(parameters,directories['multi'],rand_seed)

    # Path to the executable from multi:
    exec_command = os.path.join(directories['multi'] + '/build/bin/multi ')

    # Location of the raw output of the code:
    code_output = os.path.join(directories['multi'] + f'/output_{rand_seed}.nc')
    
    # Run the code:
    os.chdir(directories['multi'])
    # Copy the source code:
    if os.path.exists(os.path.join(res_dir + "/src")):
        rmtree(os.path.join(res_dir + "/src"))
    copytree("src", os.path.join(res_dir + "/src"))
    os.system(exec_command + f" namelist_{rand_seed} > multi.log")
    # Copy the namelist in results directory:
    copyfile(f"namelist_{rand_seed}", os.path.join(res_dir + f"/namelist_{rand_seed}"))
    # Copy the results in the results directory:
    copyfile(code_output, os.path.join(res_dir + "/output.nc"))
    # Remove the results file and namelist from multi:
    os.remove(f"namelist_{rand_seed}")
    os.remove(f"output_{rand_seed}.nc")
    os.chdir(directories['run_analysis'])
        
################################################################################
################################################################################
def run_multi_loops(rand_seed, iter_params, directories, res_dir_list):
    """Run the code multi with the parametrizations given in iter_params
    """
    #---------------------------
    # Default values for the parameters in the namelist:
    # Solver:
    nm = 3
    method = '"theoretical", "standard", "alternative"'
    na = 2
    algo = '"lanczos", "planczosif"'
    no = 2
    ni = 4
    lmp_mode = '"none"'
    test_ortho = "F"
    shutoff_type = 0
    shutoff_value = 1e-2
    interp_method = '"spectral"'
    projective_Bmatrix = "T"
    # Resolution:
    nx = "11,31,51,101"
    ny = "11,31,51,101"
    # Observations:
    nobs = 100
    sigma_obs = 0.1
    Hnl_coeff = 1
    # Background:
    sigmabvar = 0.0
    Lb = 0.12
    spvarmin = 1.0e-5
    # Miscellanous:
    new_seed = "T"
    filename = f'"output_{rand_seed}.nc"'
    #---------------------------

    res_dir_counter = 0
    #with concurrent.futures.ProcessPoolExecutor() as executor:
    #processes.append(executor.submit(run_multi, *args))
    if True:
        for no_ni, nx, nobs, sigma_obs, Hnl_coeff, sigmabvar, Lb, interp_method, projective_Bmatrix, test_ortho in iter_params:
            # args = (nm, method, na, algo, no, ni, lmp_mode, test_ortho, shutoff_type,
            #         shutoff_value, interp_method, projective_Bmatrix, nx, ny, nos, sigma_obs,
            #         Hnl_coeff, sigmabvar, Lb, spvarmin, new_seed, filename,
            #         rand_seed, directories, res_dir_list)

            res_dir = res_dir_list[res_dir_counter]
            
            no = no_ni[0]
            ni = no_ni[1]
            
            # square grid:
            ny=nx
 
            run_multi(nm, method, na, algo, no, ni, lmp_mode, test_ortho, shutoff_type,
                         shutoff_value, interp_method, projective_Bmatrix, nx, ny, nobs, sigma_obs,
                         Hnl_coeff, sigmabvar, Lb, spvarmin, new_seed, filename,
                         rand_seed, directories, res_dir)
            res_dir_counter += 1
            
################################################################################
################################################################################
def create_res_dirs(rand_seed, iter_params, results_dir_root):
    """Creates the results directories and the corresponding outer iterations lists
    """
    #---------------------------
    # Default values for the parameters in the namelist:
    # Solver:
    nm = 3
    method = '"theoretical", "standard", "alternative"'
    na = 2
    algo = '"lanczos", "planczosif"'
    no = 2
    ni = 4
    lmp_mode = '"none"'
    test_ortho = "F"
    shutoff_type = 0
    shutoff_value = 1e-2
    interp_method = '"spectral"'
    projective_Bmatrix = "T"
    # Resolution:
    nx = "11,31,51,101"
    ny = "11,31,51,101"
    # Observations:
    nobs = 100
    sigma_obs = 0.1
    Hnl_coeff = 1
    # Background:
    sigmabvar = 0.0
    Lb = 0.12
    spvarmin = 1.0e-5
    # Miscellanous:
    new_seed = "T"
    filename = f'"output_{rand_seed}.nc"'
    #---------------------------

    # Store the results directories:
    res_dir_list=[]
    # Store the outer iteraions for plotting:
    outer_iterations_list = []

    for no_ni, nx, nobs, sigma_obs, Hnl_coeff, sigmabvar, Lb, interp_method, projective_Bmatrix, test_ortho in iter_params:

        no = no_ni[0]
        ni = no_ni[1]
        
        # square grid:
        ny=nx
        
        # outer iteraions for plotting:
        outer_iterations = []
        for io in range(no+1):
            outer_iterations.append((ni)*io)
            
        res_file_name = ""

        dir_name = f'Hnl{Hnl_coeff}'
        res_file_name = dir_name
        res_dir = os.path.join(results_dir_root + '/' +  dir_name)
        
        dir_name = f'no{no}_ni{ni}'
        res_file_name += dir_name
        res_dir = os.path.join(res_dir + '_' +  dir_name)
        #if not os.path.exists(res_dir):
        #    os.mkdir(res_dir)

        dir_name = 'lmp-{}'.format(lmp_mode.replace('"',''))
        res_file_name += "_"+dir_name
        res_dir = os.path.join(res_dir + '_' +  dir_name)
        #if not os.path.exists(res_dir):
        #    os.mkdir(res_dir)

        dir_name = 'nx{}_ny{}'.format(nx.replace(',','-'), ny.replace(',','-'))
        res_file_name += "_"+dir_name
        res_dir = os.path.join(res_dir + '_' +  dir_name)
        #if not os.path.exists(res_dir):
        #    os.mkdir(res_dir)

        dir_name = f'nobs{nobs}_sigmaobs{sigma_obs}'
        res_file_name += "_"+dir_name
        res_dir = os.path.join(res_dir + '_' +  dir_name)
        #if not os.path.exists(res_dir):
        #    os.mkdir(res_dir)

        dir_name = f'sigbvar{sigmabvar}_Lb{Lb}'
        res_file_name += "_"+dir_name
        res_dir = os.path.join(res_dir + '_' +  dir_name)
        if not os.path.exists(res_dir):
            os.mkdir(res_dir)

        name_string = 'interp-{}'.format(interp_method.replace('"',''))
        res_file_name += "_"+dir_name
        res_dir = os.path.join(res_dir + '/' +  name_string)
        if not os.path.exists(res_dir):
            os.mkdir(res_dir)

        name_string = f'proj{projective_Bmatrix}'
        res_file_name += "_"+dir_name
        res_dir = os.path.join(res_dir + '/' +  name_string)
        if not os.path.exists(res_dir):
            os.mkdir(res_dir)


        name_string = f'seed_{rand_seed}'
        res_file_name += "_"+dir_name
        res_dir = os.path.join(res_dir + '/' +  name_string)
        if not os.path.exists(res_dir):
            os.mkdir(res_dir)

        res_dir_list.append(res_dir)
        outer_iterations_list.append(outer_iterations)

    return res_dir_list, outer_iterations_list
################################################################################
################################################################################        
