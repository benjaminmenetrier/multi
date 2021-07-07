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
import sys
import itertools
# Personal packages:
from analysis_tools import netcdf_extract
################################################################################

# Allow the use of tex format for the labels:
plt.rcParams.update({"text.usetex": True, "font.size" : 25})


################################################################################
################################################################################
def build_results_object(res_dir_dict, outer_iterations_dict):
    """Build usefull object which concatenates the results for each seeds,
       variable and method in order to produce statistical analysis and plots.
    """
    #keys = ['j_nl', 'jo_nl', 'jb_nl', 'j', 'jo', 'jb', 'rho_sqrt', 'beta']
    keys = ['j_full', 'jo_full', 'jb_full', 'j_quad', 'jo_quad', 'jb_quad', 'rho_sqrt', 'beta']
    
    
    ens_lanczos, ens_planczosif, ens_obj_list, ens_diff_list =  {}, {}, {}, {}

    legend = []

    
    # Save the methods name and number of configurations for plotting:
    ds0 = netcdf_extract(res_dir_dict[0][0])
    methods = ds0.groups
    
    # Recover the results for each seeds:
    # ---------------------------------------------------------------------------
    for rs in res_dir_dict:
        # Loop over the seeds:
        res_dir_list = res_dir_dict[rs]
        outer_iterations_list = outer_iterations_dict[rs]

        ens_lanczos[rs] = {}
        ens_planczosif[rs] = {}
        ens_obj_list[rs] = {}
        ens_diff_list[rs] =  {}

        if rs == 0:
            ndirs = 0

        for r, res_dir in enumerate(res_dir_list):
            rand_seed = res_dir.split('seed_')[1]

            try:
                ndirs += 1
                ens_lanczos[rs][r] = {}
                ens_planczosif[rs][r] = {}
                ens_obj_list[rs][r] = {}
                ens_diff_list[rs][r] = {}

                ds = netcdf_extract(res_dir)

                # Build the legend for plotting:
                if (r == 0):
                    avg_legend = []
                    for met in ds.groups:
                        if met == 'theoretical':
                            met  = 'th'
                        elif met == 'standard':
                            met = 'std'
                        elif met == 'alternative':
                            met = 'alt'
                        else:
                            print('problem in build_results_obj with method name: {met}')
                            return 1
                        legend.append([f'({rand_seed}){met}' + '-$B^{1/2}$', f'({rand_seed}){met}' + '-$B$'])
                        avg_legend.append([f'{met}' + '-$B^{1/2}$', f'{met}-$B$'])

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
                        ens_lanczos[rs][r][key][met] = list(itertools.chain(*ens_lanczos[rs][r][key][met]))
                        ens_planczosif[rs][r][key][met] = list(itertools.chain(*ens_planczosif[rs][r][key][met]))

                        ens_obj_list[rs][r][key][met] = [ens_lanczos[rs][r][key][met], ens_planczosif[rs][r][key][met]]

                        # Compute the difference between lanczos and planczosif:
                        diff_tmp = []
                        for j in range(len(ens_lanczos[rs][r][key][met])):
                            diff_tmp.append(ens_lanczos[rs][r][key][met][j] - ens_planczosif[rs][r][key][met][j])
                        ens_diff_list[rs][r][key][met] = diff_tmp
            except:
                print("Error while building ensemble results object with file \n", res_dir)
    
    results_obj = {}
    results_obj['objs'] = ens_obj_list, ens_diff_list
    results_obj['legends'] = legend, avg_legend
    results_obj['miscellaneous'] = methods, ndirs
    
    return results_obj 
################################################################################
################################################################################
def ensemble_compare_methods_plot(res_dir_dict, outer_iterations_dict, results_obj):
    """Plots the comparision between the different methods:
    """

    # Use the produced objects to compute the average over the seeds:

    ens_obj_list, ens_diff_list = results_obj['objs']
    legend, avg_legend = results_obj['legends']
    methods, ndirs = results_obj['miscellaneous']
    
    # keys = ['j_nl', 'jo_nl', 'jb_nl', 'j', 'jo', 'jb', 'rho_sqrt', 'beta']
    
    # labels = [r'$J^{nl}$', r'$J_o^{nl}$', r'$J_b^{nl}$', r'$J$',
    #           r'$J_o$', r'$J_b$', r'$\sqrt{\rho}$', r'$\beta$']

    keys = ['j_full', 'jo_full', 'jb_full', 'j_quad', 'jo_quad', 'jb_quad', 'rho_sqrt', 'beta']
    
    labels = [r'$J^{full}$', r'$J_o^{full}$', r'$J_b^{full}$', r'$J^{quad}$', r'$J_o^{quad}$',
              r'$J_b^{quad}$', r'$\sqrt{\rho}$', r'$\beta$']
    
    final_obj_list = {}
    final_diff_list = {}

    avg_obj_list = {}
    avg_diff_list = {}

    for r in range(ndirs):

        try:
            final_obj_list[r] = {}
            final_diff_list[r] = {}

            avg_obj_list[r] = {}
            avg_diff_list[r] = {}

            # Outer iterations list and results directory:
            outer_iterations_list = outer_iterations_dict[0][r]
            n_iterations = max(outer_iterations_list)+1

            ensemble_res_dir = res_dir_dict[0][r].split("seed")[0]
            ensemble_res_dir = os.path.join(ensemble_res_dir, "averaged")
            if not os.path.exists(ensemble_res_dir):
                os.mkdir(ensemble_res_dir)

            for k, key in enumerate(keys):
                final_obj_list[r][key] = []
                final_diff_list[r][key] = []

                avg_obj_list[r][key] = []
                avg_diff_list[r][key] = []

                for m, met in enumerate(methods):

                    # Store the results for all the seeds in order to plot each curve:
                    nseeds = 0
                    for rs in res_dir_dict:
                        final_obj_list[r][key].append(ens_obj_list[rs][r][key][met])
                        final_diff_list[r][key].append(ens_diff_list[rs][r][key][met])

                        nseeds += 1

                    # Compute the averaged values over the seeds:
                    avg_lanczos = []
                    avg_planczosif = []
                    avg_diff = []
                    for j in range(n_iterations):
                        # Loop over the iterations:
                        avg_j_lanczos = 0
                        avg_j_planczosif = 0
                        avg_j_diff = 0
                        for rs in range(nseeds):
                            avg_j_lanczos += ens_obj_list[rs][r][key][met][0][j]
                            avg_j_planczosif += ens_obj_list[rs][r][key][met][1][j]
                            avg_j_diff += ens_diff_list[rs][r][key][met][j]

                        avg_j_lanczos = avg_j_lanczos / (nseeds*1.)
                        avg_j_planczosif = avg_j_planczosif / (nseeds*1.)
                        avg_j_diff = avg_j_diff / (nseeds*1.)

                        avg_lanczos.append(avg_j_lanczos)
                        avg_planczosif.append(avg_j_planczosif)
                        avg_diff.append(avg_j_diff)

                    avg_obj_list[r][key].append([avg_lanczos , avg_planczosif])
                    avg_diff_list[r][key].append(avg_diff)

                # Plots the averaged curves:            
                ylabel1 = labels[k]
                ylabel2 = r"$B^{1/2}$ - $B$"
                xmax = 0
                for obj in avg_obj_list[r][key]:
                    xmax = max(len(obj[0]), len(obj[1]))
                x = list(range(xmax))
                xlabel = "iterations"
                avg_out_name = os.path.join(ensemble_res_dir + f'/compare_{key}.png')
                print('plotting:', avg_out_name)
                compare_plots_2N(avg_obj_list[r][key], ylabel1, avg_diff_list[r][key], ylabel2, x,
                                 xlabel, outer_iterations_list, avg_legend, avg_out_name)

                # Plots all the curves corresponding to each seed:
                ylabel1 = labels[k]
                ylabel2 = r"$B^{1/2}$ - $B$"
                xmax = 0
                for obj in final_obj_list[r][key]:
                    xmax = max(len(obj[0]), len(obj[1]))
                x = list(range(xmax))
                xlabel = "iterations"
                out_name = os.path.join(ensemble_res_dir + f'/compare_{key}_all.png')
                print('plotting:', out_name)
                compare_plots_2N(final_obj_list[r][key], ylabel1, final_diff_list[r][key], ylabel2, x,
                                 xlabel, outer_iterations_list, legend, out_name)
        except:
            print('Error: cannot plot compare_plots_2N for results_obj index: ',  r)
################################################################################
################################################################################
def linearization_check(res_dir_dict, outer_iterations_dict, results_obj):
    """Check if the linearization of the cost functions is correct:
    """
    
    ens_obj_list = results_obj['objs'][0]
    legend = results_obj['legends'][1]
    methods, ndirs = results_obj['miscellaneous']
    
    nkeys = 3
    #nl_keys = ['j_nl', 'jo_nl', 'jb_nl']
    #l_keys = ['j', 'jo', 'jb']
    
    nl_keys = ['j_full', 'jo_full', 'jb_full']
    l_keys = ['j_quad', 'jo_quad', 'jb_quad']
    
    labels = [r'$R(J)$', r'$R(J_o)$', r'$R(J_b)$']
    
    nl_avg_obj_list = {}
    l_avg_obj_list = {}
    linearization_test = {}

    for r in range(ndirs):

        nl_avg_obj_list[r] = {}
        l_avg_obj_list[r] = {}
        linearization_test[r] = {}
        
        # Outer iterations list and results directory:
        outer_iterations_list = outer_iterations_dict[0][r]
        n_iterations = max(outer_iterations_list)+1

        no = len(outer_iterations_list)-1
        ni = outer_iterations_list[1]

        try:
            ensemble_res_dir = res_dir_dict[0][r].split("seed")[0]
            ensemble_res_dir = os.path.join(ensemble_res_dir, "averaged")
            if not os.path.exists(ensemble_res_dir):
                os.mkdir(ensemble_res_dir)

            for k, key in enumerate(l_keys):

                linearization_test[r][key] = []

                for m, met in enumerate(methods):
                    nseeds = len(res_dir_dict)

                    linearization_test_lanczos = []
                    linearization_test_planczosif = []

                    linearization_test_met = []

                    for io in range(no):
                        # Loop over the outer iterations:
                        linearization_test_io_lanczos = []
                        linearization_test_io_planczosif = []

                        linearization_test_io = []

                        for ii in range(ni):
                            # Loop over the inner iterations:
                            linearization_test_ii_lanczos = 0
                            linearization_test_ii_planczosif = 0

                            # Compute the average linearization test over the seeds:
                            for rs in range(nseeds):
                                # # Lanczos:
                                # nl_diff_lanczos = ens_obj_list[rs][r][nl_keys[k]][met][0][(io+1)*ii+1] - ens_obj_list[rs][r][nl_keys[k]][met][0][(io+1)*ii]
                                # l_diff_lanczos = ens_obj_list[rs][r][key][met][0][(io+1)*ii+1] - ens_obj_list[rs][r][key][met][0][(io+1)*ii]
                                # linearization_test_ii_lanczos += nl_diff_lanczos / (l_diff_lanczos*1.)

                                # # PlanczosIF:
                                # nl_diff_planczosif = ens_obj_list[rs][r][nl_keys[k]][met][1][(io+1)*ii+1] - ens_obj_list[rs][r][nl_keys[k]][met][1][(io+1)*ii]
                                # l_diff_planczosif = ens_obj_list[rs][r][key][met][1][(io+1)*ii+1] - ens_obj_list[rs][r][key][met][1][(io+1)*ii]
                                # linearization_test_ii_planczosif += nl_diff_planczosif / (l_diff_planczosif*1.)

                                #if (ii != 0) :
                                #    print(io*ni, io*ni+ii)
                                 
                                # Lanczos:
                                if (ii != 0):
                                    nl_diff_lanczos = ens_obj_list[rs][r][nl_keys[k]][met][0][io*ni] - ens_obj_list[rs][r][nl_keys[k]][met][0][io*ni+ii]
                                    l_diff_lanczos = ens_obj_list[rs][r][key][met][0][io*ni] - ens_obj_list[rs][r][key][met][0][io*-ni+ii]
                                    linearization_test_ii_lanczos += nl_diff_lanczos / (l_diff_lanczos*1.)

                                    # PlanczosIF:
                                    nl_diff_planczosif = ens_obj_list[rs][r][nl_keys[k]][met][1][io*ni] - ens_obj_list[rs][r][nl_keys[k]][met][1][io*ni+ii]
                                    l_diff_planczosif = ens_obj_list[rs][r][key][met][1][io*ni] - ens_obj_list[rs][r][key][met][1][io*ni+ii]
                                    linearization_test_ii_planczosif += nl_diff_planczosif / (l_diff_planczosif*1.)

                            linearization_test_ii_lanczos = linearization_test_ii_lanczos / (nseeds*1.)
                            linearization_test_ii_planczosif = linearization_test_ii_planczosif / (nseeds*1.)

                            linearization_test_io_lanczos.append(linearization_test_ii_lanczos)
                            linearization_test_io_planczosif.append(linearization_test_ii_planczosif)

                            linearization_test_io.append([linearization_test_io_lanczos, linearization_test_io_planczosif])

                        linearization_test_lanczos.append(linearization_test_io_lanczos)
                        linearization_test_planczosif.append(linearization_test_io_planczosif)

                        linearization_test_met.append(linearization_test_io)

                    #linearization_test[r][key].append(linearization_test_met)
                    linearization_test[r][key].append([linearization_test_lanczos, linearization_test_planczosif])

                linearization_test[r][key] = np.array(linearization_test[r][key]).transpose(0,2,1,3)
            
                # Plots the averaged curves:            
                ylabel1 = labels[k]
                ylabel2 = r"$B^{1/2}$ - $B$"
                xlabel = "iterations"
                out_name = os.path.join(ensemble_res_dir + f'/linearization_test_{key}.png')
                print('plotting:', out_name)
                linearization_plot(linearization_test[r][key], ylabel1, linearization_test[r][key], ylabel2,
                                 xlabel, outer_iterations_list, legend, out_name)

        except:
            print("Error with linearization check")
            
################################################################################
################################################################################
def linearization_plot(obj_list, ylabel1, diff_list, ylabel2, xlabel,
                     outer_iterations_list, legend, out_name):
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
    
    ymax = obj_list[0][0][0][0]
    ymin = obj_list[0][0][0][0]
    for m, met_obj in enumerate(obj_list):
        for io, io_obj in enumerate(met_obj):
            ymax_tmp = max(max(io_obj[0]), max(io_obj[1]))
            ymin_tmp = max(min(io_obj[0]), min(io_obj[1]))
            if (ymax_tmp > ymax):    
                ymax = ymax_tmp
            if (ymin_tmp < ymin):    
                ymin = ymin_tmp
    
    ni = outer_iterations_list[1]

    color_map = cm.get_cmap('copper', len(obj_list))
    colors = color_map(range(len(obj_list)))
    # Create figure window to plot data:
    fig = plt.figure(1, figsize=(10,10))
    
    for m, met_obj in enumerate(obj_list):
        xi = 0
        xf = ni
        for io, io_obj in enumerate(met_obj):
            x = list(range(xi,xf))
            if io == 0:
                plt.plot(x, io_obj[0], color=colors[m], label=legend[m][0])
                plt.plot(x, io_obj[1], color=colors[m], linestyle='dashed', label=legend[m][1])
            else:
                plt.plot(x, io_obj[0], color=colors[m])
                plt.plot(x, io_obj[1], color=colors[m], linestyle='dashed')
            xi += ni
            xf += ni
    plt.legend(bbox_to_anchor=(0., 1.1 , 1., .102), loc='lower left',
               ncol=3, mode='expand',  borderaxespad=0)
                
    plt.ylabel(ylabel1)
    plt.vlines(outer_iterations_list, ymin, ymax, colors='blue', linestyles='dashed')

    #plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.savefig(out_name, bbox_inches='tight')
    plt.clf()
    plt.close()
################################################################################w
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
    
    color_map = cm.get_cmap('copper', len(obj_list))
    colors = color_map(range(len(obj_list)))
    # Create figure window to plot data:
    fig = plt.figure(1, figsize=(10,10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 2])

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
        #plt.legend(bbox_to_anchor=(0., 1.102 , 1., .102), loc='best',
        #           ncol=3, mode='expand',  borderaxespad=0.1)
        #handles, labels = ax1.get_legend_handles_labels()
        ax1.legend(bbox_to_anchor=(0., 1.01 , 1., .102), loc='lower left',
                   ncol=3, mode='expand',  borderaxespad=0)
        plt.subplots_adjust(top=1)
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
    plt.savefig(out_name, bbox_inches='tight')
    #fig.savefig(out_name, bbox_extra_artists=(lgd), bbox_inches='tight')
    plt.clf()
    plt.close()
################################################################################w
################################################################################
