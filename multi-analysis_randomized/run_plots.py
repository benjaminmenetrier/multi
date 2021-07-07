#!/usr/bin/env python

################################################################################
# run_plots.py

# purpose: Calls all the plot functions.

# Author: Nicolas Baillot d'Etivaux

################################################################################
# Imported packages:
import numpy as np
import os
import concurrent.futures
# Personal packages:
from analysis_1D_plot import *
from matrix_plot import *
from analysis_tools import netcdf_extract
import netCDF4 as nc
################################################################################

################################################################################
################################################################################
def run_plots_loops(res_dir_list, outer_iterations_list, extra_monitoring):
    """Run all the plots:
    """
    processes = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
    #if True:
        # Comparision between the different methods:
        for r, res_dir in enumerate(res_dir_list):
            outer_iterations = outer_iterations_list[r]
            args = (res_dir, outer_iterations, extra_monitoring)
            processes.append(executor.submit(run_all_plots, *args))
            #run_all_plots(res_dir, outer_iterations, extra_monitoring)
################################################################################
################################################################################
def run_all_plots(res_dir, outer_iterations, extra_monitoring):
    """Rune the comparision plots and the extra monitoring plots for one given 
       directory
    """
    run_compare_plots(res_dir, outer_iterations)
    run_extra_plots(res_dir, outer_iterations, extra_monitoring)
    
################################################################################
################################################################################
def run_compare_plots(res_dir, outer_iterations):
    """Runs all the comparision plots between the different methods
    """
    try:
    #if True:
        #ds = netcdf_extract(res_dir)
        with nc.Dataset(res_dir+"/output.nc","r") as ds:
            # Plots in observation space:
            obs_plot(ds, res_dir)
            hxg_plot(ds, res_dir)
            innovation_plot(ds, res_dir)

            # Plots the truth:
            x_true = np.array(ds['xt'][:])
            out_name = os.path.join(res_dir + '/xt.png')
            field_plot(x_true, out_name)

            # Comparision for 1D variables (cost functions, rho, beta ...):
            compare_methods_plot(ds, outer_iterations, res_dir)
            compare_methods_plot2(ds, outer_iterations, res_dir)

            # Comparision for 2D variables at outer loop level:
            compare_methods_2D_outer(ds, res_dir)
    except:
        print("Cannot compare methods with file \n", res_dir)

################################################################################
################################################################################
def run_extra_plots(res_dir, outer_iterations, extra_monitoring):
    """Run the plots for extra monitoring (any outer and inner loops variables,
       matrices and fields.
    """
    # Get the data:
    #ds = netcdf_extract(res_dir)
    with nc.Dataset(res_dir+"/output.nc","r") as ds:
        
        # Plots comparision between lanczos and planczosif:
        #if True:
        try:
            lanczos_vs_planczosif_plot(ds, res_dir, outer_iterations)
        except:
            print('Error, cannot plot lanczos vs planczosif for file \n', res_dir)

        # Extra monitoring at outer and inner loop level (all the variables):
        if extra_monitoring:        
            try:
                lanczos_vs_planczosif_2D_outer(ds, res_dir)
            except:
                print("Cannot plot lanczos vs planczosif 2D for file \n", res_dir)
            # Plots in model space:
            # At outer loop level:
            for met in ds.groups:
                met_dir = os.path.join(res_dir + f'/{met}')
                if not os.path.exists(met_dir):
                    os.mkdir(met_dir)

                for io in ds[met].groups:

                    x_coord = np.array(ds[met][io]['x_coord'])
                    y_coord = np.array(ds[met][io]['y_coord'])

                    background_fields = ['sigmab', 'dirac_cov', 'dirac_cor', 'dirac_cov_bis',
                                         'dirac_cor_bis', 'xb']

                    background_dir = os.path.join(met_dir + '/background')
                    if not os.path.exists(background_dir):
                        os.mkdir(background_dir)

                    for field in background_fields:
                        try:
                            matrix = np.array(ds[met][io][field][:])
                            out_name = os.path.join(background_dir + f'/{field}_{io}')
                            field_plot(matrix, out_name)
                        except:
                            print('Cannot plot matrix for ', met, io, field)

                        # At algorithms level:
                        for algo in ds[met][io].groups:
                            algo_dir = os.path.join(met_dir + f'/{algo}')
                            if not os.path.exists(algo_dir):
                                os.mkdir(algo_dir)

                            algo_fields = ['xg']
                            for field in algo_fields :
                                try:
                                    matrix = np.array(ds[met][io][algo][field][:])
                                    out_name = os.path.join(algo_dir + f'/{algo}_{field}_{io}')
                                    field_plot(matrix, out_name)
                                except:
                                    print('Cannot plot matrix for ', met, io, field)

                                # At inner loop level:
                                inner_fields = ['dx']

                                inner_dir = os.path.join(algo_dir + f'/inner_loops')
                                if not os.path.exists(inner_dir):
                                    os.mkdir(inner_dir)

                                for field in inner_fields:
                                    inner_field_dir = os.path.join(inner_dir + f'/{field}')
                                    if not os.path.exists(inner_field_dir):
                                        os.mkdir(inner_field_dir)

                                    for ii in range(ds[met][io][algo].dimensions['nimax'].size):
                                        try:
                                            dx = np.array(ds[met][io][algo][field][ii][:])
                                            out_name = os.path.join(inner_field_dir + f'/{algo}_{field}_{io}_inner_{ii}')
                                            field_plot(dx, out_name)
                                        except:
                                            print('Cannot plot matrix for ', met, io, field)
    return    
################################################################################
################################################################################
