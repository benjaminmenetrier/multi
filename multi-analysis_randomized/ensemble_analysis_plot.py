#!/usr/bin/env python

################################################################################
# ensemble_analysis_plot.py

# purpose: Contains functions plotting 1D variables averaged on the elements of
# the ensemble.

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
import itertools
# Personal packages:
from analysis_tools import netcdf_extract
################################################################################

# Allow the use of tex format for the labels:
plt.rcParams.update({"text.usetex": True, "font.size" : 20})

################################################################################
################################################################################
def ensemble_compare_methods_plot(res_dir_dict, outer_iterations_dict):
    """Plots the comparision between the different methods:
    """

    # Dirty trick here: j_nl is doubled because of a problem of figure size ...
    # ... cannot find how to fix it...
    keys = ['j_nl', 'j_nl', 'jo_nl', 'jb_nl', 'j', 'jo', 'jb', 'rho_sqrt', 'beta']
    
    labels = [r'$J^{nl}$', r'$J^{nl}$', r'$J_o^{nl}$', r'$J_b^{nl}$', r'$J$',
              r'$J_o$', r'$J_b$', r'$\sqrt{\rho}$', r'$\beta$']

    ens_lanczos, ens_planczosif, ens_obj_list, ens_diff_list =  {}, {}, {}, {}
    ens_res_dir_lists = []
    
    for rs in res_dir_dict:
        res_dir_list = res_dir_dict[rs]
        outer_iterations_list = outer_iterations_dict[rs]
        
        ens_res_dir_lists.append(res_dir_list)

        ens_lanczos[rs] = {}
        ens_planczosif[rs] = {}
        ens_obj_list[rs] = {}
        ens_diff_list[rs] =  {}

        legend = []
        for r, res_dir in enumerate(res_dir_list):
            # Loop over the different seeds:
            rand_seed = res_dir.split('seed_')[1]
            
            ens_lanczos[rs][r] = {}
            ens_planczosif[rs][r] = {}
            ens_obj_list[rs][r] = {}
            ens_diff_list[rs][r] = {}
            
            ds = netcdf_extract(res_dir)

            # Save the methods name and number of configurations for plotting:
            if (rs == 0) and (r == 0):
                methods = ds.groups
                ndirs = len(res_dir_list)
                
            # Build the legend:
            legend_met = []
            for met in ds.groups:
                for algo in ds[met]['outer_1'].groups:
                    legend_met.append(f'{met}-{algo}-{rand_seed}')
                    legend.append(legend_met)
                    
            for k, key in enumerate(keys):
                ens_lanczos[rs][r][key] = {}
                ens_planczosif[rs][r][key] = {}
                ens_obj_list[rs][r][key] = {}
                ens_diff_list[rs][r][key] = {}

                for m, met in enumerate(ds.groups):
                    ens_lanczos[rs][r][key][met] = []
                    ens_planczosif[rs][r][key][met] = []
                    ens_obj_list[rs][r][key][met] = []
                    ens_diff_list[rs][r][key][met] = []

                    for j, io in enumerate(ds[met].groups):
                        if j == 0:
                            ens_lanczos[rs][r][key][met].append(np.array(ds[met][io]['lanczos'][key][:]))
                            ens_planczosif[rs][r][key][met].append(np.array(ds[met][io]['planczosif'][key][:]))
                        else:    
                            ens_lanczos[rs][r][key][met].append(np.array(ds[met][io]['lanczos'][key][1:]))
                            ens_planczosif[rs][r][key][met].append(np.array(ds[met][io]['planczosif'][key][1:]))
                    # Reshape the lists along iteration axis:        
                    #ens_lanczos[rs][r][key][met] = list(itertools.chain(*ens_lanczos[rs][r][key][met]))
                    #ens_planczosif[rs][r][key][met] = list(itertools.chain(*ens_planczosif[rs][r][key][met]))
                    
                    ens_obj_list[rs][r][key][met] = [ens_lanczos[rs][r][key][met], ens_planczosif[rs][r][key][met]]

                    # Compute the difference between lanczos and planczosif:
                    diff_tmp = []
                    for j in range(len(ens_lanczos[rs][r][key][met])):
                        diff_tmp.append(ens_lanczos[rs][r][key][met][j] - ens_planczosif[rs][r][key][met][j])
                    ens_diff_list[rs][r][key][met] = diff_tmp


    new_obj_list = {}
    new_diff_list = {}
    
    for m, met in enumerate(methods):
        new_obj_list[met] = {}
        new_diff_list[met] = {}
        for k, key in enumerate(keys):
            new_obj_list[met][key] = {}
            new_diff_list[met][key] = {}
            for r in range(ndirs):
                new_obj_list[met][key][r] = []
                new_diff_list[met][key][r] = []
                for rs in res_dir_dict:
                    new_obj_list[met][key][r].append(ens_obj_list[rs][r][key][met])
                    new_diff_list[met][key][r].append(ens_diff_list[rs][r][key][met])
                    outer_iterations_list = outer_iterations_dict[rs][r]

                    # Results directory for the analysis averaged over the seeds:
                    ensemble_res_dir = res_dir_dict[rs][r].split("seed")[0]
                    ensemble_res_dir = os.path.join(ensemble_res_dir, "averaged")
                    if not os.path.exists(ensemble_res_dir):
                        os.mkdir(ensemble_res_dir)



                    
                ylabel1 = labels[k]
                ylabel2 = r"Lanczos - PlanczosIF"
                xmax = 0
                for obj in new_obj_list[met][key][r]:
                    xmax = max(len(obj[0]), len(obj[1]))
                x = list(range(xmax))
                xlabel = "iterations"
                out_name = os.path.join(ensemble_res_dir + f'/compare_{key}.png')
                print('plotting:', out_name)
                compare_plots_2N(new_obj_list[met][key][r], ylabel1, new_diff_list[met][key][r], ylabel2, x,
                             xlabel, outer_iterations_list, legend, out_name)

################################################################################
################################################################################
def compare_methods_plot2(ds, outer_iterations, res_dir):
    """Plots the comparision between the different methods:
    """
    diff_dict = {}

    keys=['j_nl', 'jo_nl', 'jb_nl', 'j','jo', 'jb', 'rho_sqrt', 'beta']

    labels=[r'$J^{nl}$', r'$J_o^{nl}$', r'$J_b^{nl}$', r'$J$', r'$J_o$',
            r'$J_b$', r'$\sqrt{\rho}$', r'$\beta$']

    ds_th = ds['theoretical']
    for k, key in enumerate(keys):
        diff_dict[key] = {}
        for algo in ds_th['outer_1'].groups:
            diff_dict[key][algo] = {}
            for m, met in enumerate(ds.groups):
                data_to_compare = []
                outer_iterations = []
                for j, io in enumerate(ds[met].groups):
                    data = np.array(ds[met][io][algo][key][:])
                    if j == 0:
                        data_to_compare.append(list(data[:]))
                    else:
                        data_to_compare.append(list(data[1:]))
                # Concatenate the outer iterations:
                data_to_compare=list(itertools.chain(*data_to_compare))
                diff_dict[key][algo][met] = data_to_compare
        # Plotting:
        color_map = cm.get_cmap('copper', len(ds_th['outer_1'].groups))
        colors = color_map(range(len(ds_th['outer_1'].groups)))
        for m, met in enumerate(ds.groups):
            fig = plt.figure()

            compare_dir = os.path.join(res_dir + f'/theoretical_vs_{met}')
            if not os.path.exists(compare_dir):
                os.mkdir(compare_dir)
            out_name = os.path.join(compare_dir + f'/{key}.png')
            print('plotting :', out_name)
            
            label = []
            for a, algo in enumerate(ds_th['outer_1'].groups):
                label = algo
                diff_methods_algo = []
                iterations = []
                for i in range(len(diff_dict[key][algo][met])):
                    # Compute the difference between met and "theoretical":
                    diff = diff_dict[key][algo][met][i] - diff_dict[key][algo]['theoretical'][i]
                    diff_methods_algo.append(diff)
                    iterations.append(i)
                
                plt.plot(iterations[:],diff_methods_algo[:],label=label, color=colors[a])    
            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                       ncol=2, mode="expand", borderaxespad=0.1)
            plt.vlines(outer_iterations, min(diff_methods_algo), max(diff_methods_algo),
                       colors='black', linestyles='dashed')
            # Bottom plot:
            plt.xlabel('iterations')
            plt.ylabel(labels[k] + f'({met} - theoretical)')
            plt.gcf().subplots_adjust(left=0.2, bottom=0.15)
            plt.savefig(out_name)
            plt.clf()
            plt.close()
            
################################################################################
################################################################################
def compare_plots_2N(obj_list, ylabel1, diff_list, ylabel2, x, xlabel,
                     outer_iterations, legend, out_name):
    """Produces usefull comparision with 2N plots:
    Args:
        out_name (string): Name of the output png file.
        obj_list (list)      : List containing lists to be compared together.
        yabel (string)   : Label of the y-axis.
        diff_list (list)      : List of numbers regarding which one wants to compare.
        ylabel (string)  : Label of diff_list.
        x (list)         : List of numbers for x-axis.
        xlabel (string)  : Label of x-axis
    Return:
        (void): Creates the plot and save it in the png file named by out_name arg.
    """
    print(np.shape(obj_list))
    color_map = cm.get_cmap('copper', len(obj_list))
    colors = color_map(range(len(obj_list)))
    # Create figure window to plot data:
    fig = plt.figure(1, figsize=(9,9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

    # Top plot: [obj_list1,obj_list2,..obj_listN]:
    ax1 = fig.add_subplot(gs[0])
    
    ymax = 0.
    for o, obj in enumerate(obj_list):
        ymax_tmp = max(max(obj[0]), max(obj[1]))
        if (ymax_tmp > ymax):
            ymax = ymax_tmp
    if ymax > 100:
        plt.yscale("log")

    for o, obj in enumerate(obj_list):
        ax1.plot(x[:len(obj[0])], obj[0], color=colors[o], label=legend[o][0])
        ax1.plot(x[:len(obj[1])], obj[1], color=colors[o], linestyle='dashed', label=legend[o][1])
        plt.legend(bbox_to_anchor=(0., 1.102 , 1., .102), loc='lower left',
                   ncol=3, mode="expand", borderaxespad=0.1)
        plt.subplots_adjust(top=0.8)
    ax1.set_ylabel(ylabel1)
    ax1.vlines(outer_iterations, 0, ymax, colors='blue', linestyles='dashed')
    # Bottom plot: diff_list
    ax2 = fig.add_subplot(gs[1])
    for d, diff in enumerate(diff_list):
        ax2.plot(x[:len(diff)], diff,color=colors[d])
        
    ax2.axhline(color="gray", zorder=-1)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel2)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.4f'))
    plt.savefig(out_name)
    plt.clf()
    plt.close()
    
################################################################################w
################################################################################
