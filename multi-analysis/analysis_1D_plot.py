#!/usr/bin/env python

################################################################################
# analysis_1D_plot.py

# purpose: Contains functions plotting 1D variables.
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
def compare_plots(obj1, obj2, ylabel12, obj3, ylabel3, x, xlabel,
                  outer_iterations, legend,out_name):
    """Produces usefull comparision plots bteween obj1 and obj2:
    Args:
        out_name (string): Name of the output png file.
        obj1, obj2 (list): The two lists of numbers to compare.
        yabel12 (string) : Label of the two objects obj1 and obj2.
        obj3 (list)      : List of numbers regarding which one wants to compare
                           obj1 and obj2 (for example obj1-obj2).
        ylabel3 (string) : Label of obj3.
        x (list)         : List of numbers for x-axis.
        xlabel (string)  : Label of x-axis.
        legend (list(string): Labels for obj1 and obj2.
    Return:
        (void): Creates the plot and save it in the png file named by out_name arg.
    """
    print("plotting:", out_name)
    # Create figure window to plot data
    fig = plt.figure(1, figsize=(9,9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[6, 3])

    # Top plot: obj1 vs obj2
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(x[:len(obj1)], obj1, color='black', label=legend[0])
    ax1.plot(x[:len(obj2)], obj2, color='black', linestyle='dashed', label=legend[1])
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
               ncol=2, mode="expand", borderaxespad=0.)

    plt.subplots_adjust(hspace=0.5)
    
    ymax = max(max(obj1), max(obj2))
    if ymax > 100:
        plt.yscale("log")
    ax1.set_ylabel(ylabel12)
    start, end = ax1.get_xlim()
    ax1.vlines(outer_iterations, 0, ymax, colors='blue', linestyles='dashed')

    # Bottom plot: obj3
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(x[:len(obj3)], obj3, color='black')
    ax2.axhline(color="gray", zorder=-1)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel3)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.savefig(out_name)
    plt.clf()
    plt.close()
    
################################################################################
################################################################################
def lanczos_vs_planczosif_plot(ds, res_dir, outer_iterations):
    """ Plot the J, Jb and Jo:
        Args:
            ds (netcdf object): results file.  
            res_dir (string): result directory in which is stored the netcdf 
                              output file.
            outer_iterations (list): position of the outer iterations for
                                     plotting.
        Return:
            (void): Plots the cost functions.
    """
    lanczos, planczosif={}, {}
    keys = ['j_nl', 'jo_nl', 'jb_nl', 'j', 'jo', 'jb', 'rho_sqrt', 'beta']

    labels = [r'$J^{nl}$', r'$J_o^{nl}$', r'$J_b^{nl}$', r'$J$', r'$J_o$',
              r'$J_b$', r'$\sqrt{\rho}$', r'$\beta$']

    for k in keys:
        lanczos[k] = []
        planczosif[k] = []
    
    for outer in ds.groups:
        for k in keys:
            lanczos[k].append(np.array(ds[outer]['lanczos'][k][:]))
            planczosif[k].append(np.array(ds[outer]['planczosif'][k][:]))
            
    diff = {}    
    for k,key in enumerate(lanczos.keys()):
        diff[key] = []
        lanczos[key] = np.reshape(lanczos[key], -1)
        planczosif[key] = np.reshape(planczosif[key], -1)
        for i in range(len(lanczos[key])):
            diff[key].append(lanczos[key][i] - planczosif[key][i])
    
        # Plot the results:
        out_name = os.path.join(res_dir + f'/compare_{key}.png')
        print("plotting:", out_name)
        ylabel12 = labels[k]
        ylabel3 = 'diff'
        xmax = max(len(lanczos[key]), len(planczosif[key]))
        x = range(xmax)
        xlabel = 'iterations'
        legend = ['Lanczos', 'PlanczosIF']
        compare_plots(lanczos[key], planczosif[key], ylabel12, diff[key], ylabel3, x,
                      xlabel, outer_iterations, legend,out_name)
        
################################################################################
################################################################################
def compare_methods_plot(compare_methods_data, methods_list, outer_iterations,
                         compare_methods_dir):
    """Plots the comparision between the different methods:
    """
    lanczos, planczosif, obj_list, diff_list={}, {}, {}, {}

    keys=['j_nl', 'jo_nl', 'jb_nl', 'j', 'jo', 'jb', 'rho_sqrt', 'beta']

    labels=[r'$J^{nl}$', r'$J_o^{nl}$', r'$J_b^{nl}$', r'$J$', r'$J_o$', r'$J_b$',
            r'$\sqrt{\rho}$', r'$\beta$']

    legend = []
    for met in methods_list:
        legend_met = []
        for space in ['control', 'model']:
            legend_met.append(met + '-' + space)
        legend.append(legend_met)
    
    for k,key in enumerate(keys):
        lanczos[key] = []
        planczosif[key] = []
        obj_list[key] = []
        diff_list[key] = []
        for i,ds in enumerate(compare_methods_data):
            lanczos[key].append([])
            planczosif[key].append([])
            for io in ds.groups:
                lanczos[key][i].append(np.array(ds[io]['lanczos'][key][:]))
                planczosif[key][i].append(np.array(ds[io]['planczosif'][key][:]))
            lanczos[key][i]=np.reshape(lanczos[key][i], -1)
            planczosif[key][i]=np.reshape(planczosif[key][i], -1)
            obj_list[key].append([lanczos[key][i],planczosif[key][i]])
            # Compute the difference between lanczos and planczosif:
            diff_tmp = []
            for j in range(len(lanczos[key][i])):
                diff_tmp.append(lanczos[key][i][j] - planczosif[key][i][j])
            diff_list[key].append(diff_tmp)
        ylabel1 = labels[k]
        ylabel2 = r"lanczos - PlanczosIF"
        xmax = 0
        for obj in obj_list[key]:
            xmax = max(len(obj[0]), len(obj[1]))
        x = list(range(xmax))
        xlabel = "iterations"
        out_name = compare_methods_dir + f'/compare_{key}.png'
        print('plotting:', out_name)
        compare_plots_2N(obj_list[key], ylabel1, diff_list[key], ylabel2, x,
                         xlabel, outer_iterations, legend,out_name)
        
################################################################################
################################################################################
def compare_methods_plot2(compare_methods_data, methods_list,
                          outer_iterations, compare_methods_dir):
    """Plots the comparision between the different methods:
    """
    diff_dict = {}

    keys=['j_nl', 'jo_nl', 'jb_nl', 'j','jo', 'jb', 'rho_sqrt', 'beta']

    labels=[r'$J^{nl}$', r'$J_o^{nl}$', r'$J_b^{nl}$', r'$J$', r'$J_o$',
            r'$J_b$', r'$\sqrt{\rho}$', r'$\beta$']

    ds_th = compare_methods_data[0]
    for k, key in enumerate(keys):
        diff_dict[key] = {}
        for algo in ds_th['outer_1'].groups:
            diff_dict[key][algo] = []
            for m, met in enumerate(methods_list):
                ds = compare_methods_data[m]
                data_to_compare = []
                outer_iterations = []
                for j, io in enumerate(ds.groups):
                    data_to_compare.append(ds[io][algo][key][:])
                    # Get the outer iterations for plotting:
                    #inner_iterations = len(ds[io][algo][key][:])
                    #outer_iterations.append((inner_iterations-1) * (j+1))
                # Concatenate the outer iterations:
                data_to_compare = np.reshape(data_to_compare,-1)
                diff_dict[key][algo].append(data_to_compare)
                
        # Plotting:
        color_map = cm.get_cmap('copper', len(ds_th['outer_1'].groups))
        colors = color_map(range(len(ds_th['outer_1'].groups)))
        for m, met in enumerate(methods_list):
            fig = plt.figure()
            #fig = plt.figure(1, figsize=(9,9))
            #gs = gridspec.GridSpec(2, 1, height_ratios=[6, 3])
            # Top plot:
            #ax1 = fig.add_subplot(gs[0])
            out_name = os.path.join(compare_methods_dir + f'/theoretical_vs_{met}/{key}.png')
            label = []
            for a, algo in enumerate(ds_th['outer_1'].groups):
                label = algo
                diff_methods_algo = []
                iterations = []
                for i in range(len(diff_dict[key][algo][m])):
                    # Compute the difference between met and "theoretical":
                    diff = diff_dict[key][algo][m][i] - diff_dict[key][algo][0][i]
                    diff_methods_algo.append(diff)
                    iterations.append(i)
                
                plt.plot(iterations[:],diff_methods_algo[:],label=label, color=colors[a])
                
            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                       ncol=2, mode="expand", borderaxespad=0.1)
            
            plt.vlines(outer_iterations, min(diff_methods_algo), max(diff_methods_algo),
                       colors='black', linestyles='dashed')
            # Bottom plot:
            # ax2 = fig.add_subplot(gs[1])
            # ax2.plot(iterations[:],diff,color='black')
            # ax2.axhline(color="gray", zorder=-1)
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
    color_map = cm.get_cmap('copper', len(obj_list))
    colors = color_map(range(len(obj_list)))
    # Create figure window to plot data:
    fig = plt.figure(1, figsize=(9,9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

    # Top plot: [obj_list1,obj_list2,..obj_listN]:
    ax1 = fig.add_subplot(gs[0])
    plt.yscale("log")
    ymax = 0.
    for o, obj in enumerate(obj_list):
        ymax_tmp = max(max(obj[0]), max(obj[1]))
        if (ymax_tmp > ymax):
            ymax = ymax_tmp
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



# Obsolete functions (but might be used later with LMPs):
################################################################################
# /!\ Obsolete: (tmp)
def lmp_compare(out_names, lmp_to_compare, column_of_interest, ylabel1, ylabel2,
                outer_iterations_list):
    """Compare spectral and ritz lmp modes:
    """
    try:
        legend=[]
        for lmp_mode in ['ritz','spectral','none']:
            legend.append([lmp_mode+'-model',lmp_mode+'-control'])
            
        for r,res_dirs in enumerate(lmp_to_compare):
            diff_list=[]
            obj_list=[]
            for res_dir in res_dirs:
                res1=np.genfromtxt(res_dir+'PlanczosIF_model_space.dat', comments='#')    
                res2=np.genfromtxt(res_dir+'lanczos_control_space.dat', comments='#')
                diff=[]
                l_iter=min(len(res1[0]),len(res2[0]))
                for c in range(len(res1)):
                    diff_col=[]
                    for l in range(l_iter):
                        diff_col.append(res1[c,l]-res2[c,l]) 
                    diff.append(diff_col)
                diff=np.array(diff)
                diff_list.append(diff[:,column_of_interest])
                obj_list.append([res1[:,column_of_interest],res2[:,column_of_interest]])

            # maybe dirtyish but...
            itot=list(range(len(diff)))
            print("plotting lmp comparision for:\n",out_names[r],"\n for column:", column_of_interest)
            x=itot
            xlabel='iterations'
            out_name=out_names[r]
            outer_iterations=outer_iterations_list[r]
            compare_plots_2N(out_name, obj_list, ylabel1, diff_list, ylabel2, x,
                             xlabel, outer_iterations, legend)
    except:
        print("Error with lmp comparision of:\n",res_dirs,"\n")
##################################################################################
################################################################################
# /!\ Obsolete (tmp):
def vec_plot(results_file,column_of_interest,label,out_name):
    """Plots 1D vectors:
    """
    print("plotting:",out_name)
    # Get the data:
    outer_vectors=np.genfromtxt(results_file, comments='#')
    io=1
    vec=[]
    indices=[]
    vec_io=[]
    indices_io=[]
    for i  in range(len(outer_vectors)):
        # Add the elemnts of the outer vector for iteration io:
        if int(outer_vectors[i,0])==io:
            vec_io.append(outer_vectors[i,column_of_interest])
            indices_io.append(outer_vectors[i,1])
        # Add the outer vector for outer iteration io to vec:    
        else:
            io=io+1
            vec.append(vec_io)
            indices.append(indices_io)
            vec_io=[]
            indices_io=[]
            vec_io.append(outer_vectors[i,column_of_interest])
            indices_io.append(outer_vectors[i,1])
    # Add the last outer vector:
    vec.append(vec_io)
    indices.append(indices_io)
    
    if not len(vec)==1:
        fig, subplots = plt.subplots(len(vec),1)
        for i, ax in enumerate(subplots):
            plt.ticklabel_format(axis="y", style="sci", scilimits=(-3,3))
            ax.plot(indices[i][:],vec[i][:],color='blue')
            ax.set_ylabel(label+'_{}$'.format(i))
        plt.tight_layout()
        fig.align_ylabels(subplots[:])    
        plt.subplots_adjust(hspace=0.5)
    else:
        fig=plt.figure()
        plt.plot(indices[0][:],vec[0][:],color='blue')
        plt.ylabel(label+r'_{}$'.format(0))
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.savefig(out_name)
    plt.clf()
    plt.close()
################################################################################
################################################################################
