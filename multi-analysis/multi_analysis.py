#!/usr/bin/env python

################################################################################
# multi-analysis.py

# purpose: This code plots the results of the code multi.
# Author: Nicolas Baillot d'Etivaux

################################################################################
# Imported packages:
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec 
import matplotlib.ticker as ticker
from matplotlib import cm
import numpy as np


# Allow the use of tex format for the labels:
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "sans-serif",
#     "font.sans-serif": ["Helvetica"]})
plt.rcParams.update({"text.usetex": True, "font.size" : 15})

# To use this module in command line:
#results_directory=str(sys.argv[1])

################################################################################
def compare_plots(out_name,obj1,obj2,ylabel12,obj3,ylabel3,x,xlabel,io):
    """Produces usefull comparision plots bteween obj1 and obj2:
    Args:
        out_name (string): Name of the output png file.
        obj1, obj2 (list): The two lists of numbers to compare.
        yabel12 (string) : Label of the two objects obj1 and obj2.
        obj3 (list)      : List of numbers regarding which one wants to compare
                           obj1 and obj2 (for example obj1-obj2).
        ylabel3 (string) : Label of obj3.
        x (list)         : List of numbers for x-axis.
        xlabel (string)  : Label of x-axis
    Return:
        (void): Creates the plot and save it in the png file named by out_name arg.
    """
    
    # Create figure window to plot data
    fig = plt.figure(1, figsize=(9,9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

    # Top plot: obj1 vs obj2
    ax1 = fig.add_subplot(gs[0])
    plt.yscale("log")
    ax1.plot(x[:],obj1,color='blue')
    ax1.plot(x[:],obj2,color='blue',linestyle='dashed')
    ymax=max(max(obj1),max(obj2))
    ax1.set_ylabel(ylabel12)
    start, end = ax1.get_xlim()
    #ax1.xaxis.set_ticks(np.arange(start, end, 1.))
    ax1.vlines(io, 0, ymax, colors='black', linestyles='dashed')

    # Bottom plot: obj3
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(x[:],obj3,color='black')
    ax2.axhline(color="gray", zorder=-1)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel3)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.4f'))
    plt.savefig(out_name)
    plt.clf()
################################################################################



################################################################################
def compare_plots_2N(out_name,obj_list,ylabel1,diff_list,ylabel2,x,xlabel,io,legend):
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

    color_map=cm.get_cmap('copper', len(obj_list))
    colors=color_map(range(len(obj_list)))
    # Create figure window to plot data:
    fig = plt.figure(1, figsize=(9,9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])
    
    # Top plot: [obj_list1,obj_list2,..obj_listN]:
    ax1 = fig.add_subplot(gs[0])
    plt.yscale("log")
    ymax=0.
    for o,obj in enumerate(obj_list):
        ymax_tmp=max(max(obj[0]),max(obj[1]))
        if (ymax_tmp>ymax):
            ymax=ymax_tmp
        ax1.plot(x[:],obj[0],color=colors[o],label=legend[o][0])
        ax1.plot(x[:],obj[1],color=colors[o],linestyle='dashed',label=legend[o][1])
        #ax1.plot(x[:],obj[2],color=c[o],linestyle='dotted',label=legend[o][2])
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=3, mode="expand", borderaxespad=0.)
    ax1.set_ylabel(ylabel1)
    #start, end = ax1.get_xlim()
    #ax1.xaxis.set_ticks(np.arange(start, end, 1.))
    ax1.vlines(io, 0, ymax, colors='black', linestyles='dashed')

    # Bottom plot: diff_list
    ax2 = fig.add_subplot(gs[1])
    for d,diff in enumerate(diff_list):
        ax2.plot(x[:],diff,color=colors[d])
    ax2.axhline(color="gray", zorder=-1)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel2)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.4f'))
    plt.savefig(out_name)
    plt.clf()
################################################################################



################################################################################
def evolution_plot(results_directory,io):
    """Plots J, Jo and Jb for lanczos and PlanczosIf, as well as their difference.
    """
    # Get the results files from the results directory and store them in lists:
    lanczos=np.genfromtxt(results_directory+'lanczos_control_space.dat', comments='#')
    PlanczosIF=np.genfromtxt(results_directory+'PlanczosIF_model_space.dat', comments='#')
    diff=np.genfromtxt(results_directory+'lanczos_control_vs_PlanczosIF_model.dat', comments='#')

    # lanczos=np.atleast_2d(np.genfromtxt(results_directory+'/lanczos_control_space.dat', comments='#'))
    # PlanczosIF=np.atleast_2d(np.genfromtxt(results_directory+'/PlanczosIF_model_space.dat', comments='#'))
    # diff=np.atleast_2d(np.genfromtxt(results_directory+'/lanczos_control_vs_PlanczosIF_model.dat', comments='#'))

    # Total iterations (outer x inner):
    itot=list(range(len(lanczos)))

    #--------------- Plots for J -----------------------------
    out_name=results_directory+"/lanczos_vs_PlanczosIF_J.png"
    compare_plots(out_name,lanczos[:,3],PlanczosIF[:,3],r'$J=J_b+J_o$',
                 diff[:,3],r'$\Delta J$',
                 itot,r'iterations',io)

    #--------------- Plots for Jb ----------------------------
    out_name=results_directory+"/lanczos_vs_PlanczosIF_Jb.png"
    compare_plots(out_name, lanczos[:,4],PlanczosIF[:,4],r'$J_b$',
                 diff[:,4],r'$\Delta J_b$',
                 itot,'iterations',io)

    #--------------- Plots for Jo ----------------------------
    out_name=results_directory+"/lanczos_vs_PlanczosIF_Jo.png"
    compare_plots(out_name,lanczos[:,5],PlanczosIF[:,5],r'$J_o$',
                 diff[:,5],r'$\Delta J_o$',
                 itot,'iterations',io)
################################################################################


################################################################################
def multi_plot(res_dir_list,outer_iterations):
    """Run the analysis over the results:
    """
    for r,res_dir in enumerate(res_dir_list):
        print("plotting for: ",res_dir)
        try:
            evolution_plot(res_dir,outer_iterations[r])
        except:
            print("Error with directory:",res_dir)
            pass
################################################################################


################################################################################
def lmp_compare(out_names,lmp_to_compare,outer_iterations_list):
    """Compare spectral and ritz lmp modes:
    """

    labels=[]
    legend=[]
    for lmp_mode in ['ritz','spectral','none']:
        labels.append(lmp_mode)
        legend.append([lmp_mode+'-model',lmp_mode+'-control'])
        

    for r,res_dirs in enumerate(lmp_to_compare):
        diff_list=[]
        obj_list=[]
        for res_dir in res_dirs:
            res1=np.genfromtxt(res_dir+'lanczos_control_vs_PlanczosIF_model.dat', comments='#')
            res2=np.genfromtxt(res_dir+'PlanczosIF_model_space.dat', comments='#')    
            res3=np.genfromtxt(res_dir+'lanczos_control_space.dat', comments='#')
            diff_list.append(res1[:,3])
            obj_list.append([res2[:,3],res3[:,3]])
        
        # maybe dirtyish but...
        itot=list(range(len(res1)))
        print("plotting lmp comparision for:\n",out_names[r])
        ylabel1=r'$J=J_o+J_b$'
        ylabel2=r'$\Delta J$'
        x=itot
        xlabel='iterations'
        out_name=out_names[r]
        outer_iterations=outer_iterations_list[r]
        try:
            compare_plots_2N(out_name,obj_list,ylabel1,diff_list,ylabel2,x,xlabel,outer_iterations,legend)
        except:
            print("Error with lmp comparision of:",res_dirs)
################################################################################

################################################################################
def diff_plot(out_names,param_list,res_dir_list):

    labels=[]
    legend=[]
    for lmp_mode in ['ritz','spectral','none']:
        labels.append(lmp_mode)
        legend.append([lmp_mode+'-model',lmp_mode+'-control'])    
    ylabel=r'\delta J'
    xlabel=r'\N_obs'
    color_map=cm.get_cmap('copper', 3)
    colors=color_map(range(3))

    for r,res_dir in enumerate(res_dir_list):
        diff_list=[]
        for res_dir in res_dirs:
            res=np.genfromtxt(res_dir+'lanczos_control_vs_PlanczosIF_model.dat', comments='#')
            diff_list.append(res[-1,3])
            
    # Create figure window to plot data:
    fig = plt.figure(1, figsize=(9,9))
    
    # Top plot: [obj_list1,obj_list2,..obj_listN]:
    ax1.plot(param_list[:],diff_list[:],color=colors[o],label=legend[o][0])
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                               ncol=3, mode="expand", borderaxespad=0.)
    ax1.set_ylabel(ylabel)
    ax1.set_xlabel(xlabel)
