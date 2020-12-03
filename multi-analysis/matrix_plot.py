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
import netCDF4 as nc
from analysis_tools import netcdf_extract

# Allow the use of tex format for the labels:
plt.rcParams.update({"text.usetex": True, "font.size" : 15})

################################################################################
################################################################################
def obs_plot(ds,res_dir):
    """Plots the observations:
    """
    out_name=res_dir+'/obs.png'
    x_obs=np.array(ds['x_obs'][:])
    y_obs=np.array(ds['y_obs'][:])
    obs_val=np.array(ds['obs_val'][:])
    
    print("plotting:", out_name)
    fig = plt.figure()
    plt.scatter(x_obs,y_obs,c=obs_val,cmap=plt.get_cmap('copper'))
    plt.colorbar()
    plt.savefig(out_name)
    plt.clf()
    plt.close()
################################################################################
################################################################################
def coord_plot(ds,res_dir):
    """Check the model grid (should plot a straight line x=y):
    """
    for io in ds.groups:
        out_name=res_dir+'/grid_coord_{}.png'.format(io)
        x_coord=np.array(ds[io]['x_coord'][:])
        y_coord=np.array(ds[io]['y_coord'][:])
        fig = plt.figure()
        plt.scatter(x_coord,y_coord)
        plt.savefig(out_name)
        plt.clf()
        plt.close()
################################################################################
################################################################################
def hxg_plot(ds,res_dir):
    """Plots hxg
    """
    x_obs=np.array(ds['x_obs'][:])
    y_obs=np.array(ds['y_obs'][:])
    for io in ds.groups:
        for algo in ds[io].groups:
            hxg=np.array(ds[io][algo]['hxg'][:])
            out_name=res_dir+'/'+algo+'_hxg_'+io+'.png'
            print("plotting:", out_name)
            fig=plt.figure()
            plt.scatter(x_obs,y_obs,c=hxg,cmap=plt.get_cmap('copper'))
            plt.colorbar()
            plt.savefig(out_name)
            plt.clf()
            plt.close()
################################################################################
################################################################################
def innovation_plot(ds,res_dir):
    """Plots the innovation
    """
    x_obs=np.array(ds['x_obs'][:])
    y_obs=np.array(ds['y_obs'][:])
    for io in ds.groups:
        for algo in ds[io].groups:
            d=np.array(ds[io][algo]['d'][:])
            out_name=res_dir+'/'+algo+'_d_'+io+'.png'
            print("plotting:", out_name)
            fig=plt.figure()
            plt.scatter(x_obs,y_obs,c=d,cmap=plt.get_cmap('copper'))
            plt.colorbar()
            plt.savefig(out_name)
            plt.clf()
            plt.close()
################################################################################
################################################################################
def field_plot(matrix,out_name):
    """Represents matrix using matshow.
    """
    print('plotting:',out_name)
    fig=plt.figure()
    plt.matshow(np.array(matrix),cmap=plt.get_cmap('copper'))
    plt.colorbar()
    plt.savefig(out_name)
    plt.clf()
    plt.close()
################################################################################
################################################################################
def field_plot2(matrix,x_coord,y_coord,out_name):
    """Represents matrix using scatterplot.
    """
    print('plotting:',out_name)
    fig=plt.figure()
    val=[]
    for i in range(len(x_coord)):
        val.append(matrix[i][j])        
    plt.scatter(x_coord,y_coord,c=val,cmap=plt.get_cmap('copper'))
    plt.colorbar()
    plt.savefig(out_name)
    plt.clf()
    plt.close()
################################################################################
################################################################################
def compare_methods_2D_outer(compare_methods_data,methods_list,compare_methods_dir):
    """Plots the comparision between the different methods:
    """
    diff_dict={}
    keys=['xg']
    labels=[r'$x^g$']

    ds_th=compare_methods_data[0]
    for io in ds_th.groups:
        for algo in ds_th[io].groups:
            diff_dict[algo]={}
            for key in keys:
                diff_dict[algo][key]=[]
                for m,met in enumerate(methods_list):
                    ds=compare_methods_data[m]
                    data_to_compare=np.array(ds[io][algo][key][:])
                    diff_dict[algo][key].append(np.array(data_to_compare))
                # Compute the difference as a 2D matrix -- diff = (method - theoretical)
                for m,met in enumerate(methods_list):
                    diff_matrix=np.zeros(np.shape(diff_dict[algo][key][m]))
                    for line in range(len(diff_dict[algo][key][m])):
                        for col in range(len(diff_dict[algo][key][m][line])):
                            diff=diff_dict[algo][key][m][line][col]-diff_dict[algo][key][0][line][col]
                            diff_matrix[line][col]=diff
                    # Plot the difference:
                    out_name = compare_methods_dir+'/theoretical_vs_{}/{}_{}_{}.png'.format(met,algo,key,io)
                    print('plotting:',out_name)
                    field_plot(diff_matrix,out_name)
################################################################################
################################################################################
