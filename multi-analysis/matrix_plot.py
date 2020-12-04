#!/usr/bin/env python

################################################################################
# matrix_plot.py

# purpose: Contains function plotting 2D fields.
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
from analysis_tools import netcdf_extract

# Allow the use of tex format for the labels:
plt.rcParams.update({"text.usetex": True, "font.size" : 15})

################################################################################
################################################################################
def obs_plot(ds, res_dir):
    """Plots the observations:
    """
    out_name = os.path.join(res_dir + '/obs.png')
    x_obs = np.array(ds['x_obs'][:])
    y_obs = np.array(ds['y_obs'][:])
    obs_val = np.array(ds['obs_val'][:])
    
    print("plotting:", out_name)
    fig = plt.figure()
    plt.scatter(x_obs, y_obs, c=obs_val, cmap=plt.get_cmap('copper'))
    plt.colorbar()
    plt.savefig(out_name)
    plt.clf()
    plt.close()
################################################################################
################################################################################
def coord_plot(ds, res_dir):
    """Check the model grid (should plot a straight line x=y):
    """
    for io in ds.groups:
        out_name = os.path.join(res_dir + f'/grid_coord_{io}.png')
        x_coord = np.array(ds[io]['x_coord'][:])
        y_coord = np.array(ds[io]['y_coord'][:])
        fig = plt.figure()
        plt.scatter(x_coord, y_coord)
        plt.savefig(out_name)
        plt.clf()
        plt.close()
################################################################################
################################################################################
def hxg_plot(ds, res_dir):
    """Plots hxg
    """
    x_obs = np.array(ds['x_obs'][:])
    y_obs = np.array(ds['y_obs'][:])
    hxg_dir = os.path.join(res_dir + '/hxg')
    if not os.path.exists(hxg_dir):
        os.mkdir(hxg_dir)
    for io in ds.groups:
        for algo in ds[io].groups:
            hxg = np.array(ds[io][algo]['hxg'][:])
            out_name = os.path.join(hxg_dir + f'/{algo}_hxg_{io}.png')
            print("plotting:", out_name)
            fig = plt.figure()
            plt.scatter(x_obs, y_obs, c=hxg, cmap=plt.get_cmap('copper'))
            plt.colorbar()
            plt.savefig(out_name)
            plt.clf()
            plt.close()
################################################################################
################################################################################
def innovation_plot(ds, res_dir):
    """Plots the innovation
    """
    x_obs = np.array(ds['x_obs'][:])
    y_obs = np.array(ds['y_obs'][:])
    innov_dir = os.path.join(res_dir + '/innovation')
    if not os.path.exists(innov_dir):
        os.mkdir(innov_dir)
    for io in ds.groups:
        for algo in ds[io].groups:
            d = np.array(ds[io][algo]['d'][:])
            out_name = os.path.join(innov_dir + f'/{algo}_innovation_{io}.png')
            print("plotting:", out_name)
            fig = plt.figure()
            plt.scatter(x_obs, y_obs, c=d, cmap=plt.get_cmap('copper'))
            plt.colorbar()
            plt.savefig(out_name)
            plt.clf()
            plt.close()
################################################################################
################################################################################
def field_plot(matrix, out_name):
    """Represents matrix using matshow.
    """
    print('plotting:', out_name)
    fig = plt.figure()
    plt.matshow(np.array(matrix), cmap=plt.get_cmap('copper'))
    plt.colorbar()
    plt.savefig(out_name)
    plt.clf()
    plt.close()
################################################################################
################################################################################
def field_plot2(matrix, x_coord, y_coord, out_name):
    """Represents matrix using scatterplot.
    """
    print('plotting:',out_name)
    fig = plt.figure()
    val = []
    for i in range(len(x_coord)):
        val.append(matrix[i][j])        
    plt.scatter(x_coord, y_coord, c=val, cmap=plt.get_cmap('copper'))
    plt.colorbar()
    plt.savefig(out_name)
    plt.clf()
    plt.close()
################################################################################
################################################################################
def compare_methods_2D_outer(compare_methods_data, methods_list,
                             compare_methods_dir):
    """Plots the comparision between the different methods:
    """
    diff_dict = {}

    keys = ['xg']
    labels = [r'$x^g$']

    ds_th = compare_methods_data[0]
    for io in ds_th.groups:
        for algo in ds_th[io].groups:
            diff_dict[algo] = {}
            for key in keys:
                diff_dict[algo][key] = []
                for m, met in enumerate(methods_list):
                    ds = compare_methods_data[m]
                    data_to_compare = np.array(ds[io][algo][key][:])
                    diff_dict[algo][key].append(np.array(data_to_compare))
                # Compute the difference as a 2D matrix -- diff = (method - theoretical)
                for m,met in enumerate(methods_list):
                    diff_matrix = np.zeros(np.shape(diff_dict[algo][key][m]))
                    for line in range(len(diff_dict[algo][key][m])):
                        for col in range(len(diff_dict[algo][key][m][line])):
                            diff = diff_dict[algo][key][m][line][col] - diff_dict[algo][key][0][line][col]
                            diff_matrix[line][col] = diff
                    # Plot the difference:
                    out_name = f'/theoretical_vs_{met}/{algo}_{key}_{io}.png'
                    out_name = os.path.join(compare_methods_dir + out_name)
                    print('plotting:', out_name)
                    field_plot(diff_matrix, out_name)
################################################################################
################################################################################
