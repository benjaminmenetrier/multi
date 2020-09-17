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
def compare_plot(out_name,obj1,obj2,ylabel12,obj3,ylabel3,x,xlabel,io):
    """Produces usefull comparision plots:
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
    ax1.plot(x[:],obj2,color='green')
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
    ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.4f'))
    plt.savefig(out_name)
    plt.clf()
################################################################################


################################################################################
def multi_plot(results_directory,io):

    # Get the results files from the results directory and store them in lists:
    lanczos=np.atleast_2d(np.genfromtxt(results_directory+'/lanczos_control_space.dat', comments='#'))
    PlanczosIF=np.atleast_2d(np.genfromtxt(results_directory+'/PlanczosIF_model_space.dat', comments='#'))
    diff=np.atleast_2d(np.genfromtxt(results_directory+'/lanczos_control_vs_PlanczosIF_model.dat', comments='#'))
    
    # Total iterations (outer x inner):
    itot=list(range(len(lanczos)))

    #--------------- Plots for J -----------------------------
    out_name=results_directory+"/lanczos_vs_PlanczosIF_J.png"
    compare_plot(out_name,lanczos[:,3],PlanczosIF[:,3],r'$J=J_b+J_o$',
                 diff[:,3],r'$\Delta J$',
                 itot,r'iterations',io)

    # # Create figure window to plot data
    # fig = plt.figure(1, figsize=(9,9))
    # gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

    # # Top plot:
    # ax1 = fig.add_subplot(gs[0])
    # plt.yscale("log")
    # ax1.plot(itot[:],lanczos[:,3],color='black')
    # ax1.plot(itot[:],PlanczosIF[:,3],color='green')
    # ymax=max(max(lanczos[:,3]),max(PlanczosIF[:,3]))
    # ax1.vlines(io, 0, ymax, colors='black', linestyles='dashed')
    # ax1.set_ylabel('J=Jb+Jo')

    # # Bottom plot: residuals
    # ax2 = fig.add_subplot(gs[1])
    # ax2.plot(itot[:],diff[:,3],color='black')
    # ax2.axhline(color="gray", zorder=-1)
    # ax2.set_xlabel('iterations')
    # ax2.set_ylabel('deltaJ')
    # plt.savefig(out_name)
    # plt.clf()

    #--------------- Plots for Jb ----------------------------
    out_name=results_directory+"/lanczos_vs_PlanczosIF_Jb.png"
    compare_plot(out_name, lanczos[:,4],PlanczosIF[:,4],r'$J_b$',
                 diff[:,4],r'$\Delta J_b$',
                 itot,'iterations',io)

    # # Create figure window to plot data
    # fig = plt.figure(1, figsize=(9,9))
    # gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

    # # Top plot:
    # ax1 = fig.add_subplot(gs[0])
    # plt.yscale("log")
    # ax1.plot(itot[:],lanczos[:,4],color='black')
    # ax1.plot(itot[:],PlanczosIF[:,4],color='green')
    # # ax1.set_xlabel('iterations')
    # ymax=max(max(lanczos[:,4]),max(PlanczosIF[:,4]))
    # ax1.vlines(io, 0, ymax, colors='black', linestyles='dashed')
    # ax1.set_ylabel('Jb')

    # # Bottom plot: residuals
    # ax2 = fig.add_subplot(gs[1])
    # ax2.plot(itot[:],diff[:,4],color='black')
    # ax2.axhline(color="gray", zorder=-1)
    # ax2.set_xlabel('iterations')
    # ax2.set_ylabel('deltaJb')
    # plt.savefig(out_name)
    # plt.clf()

    #--------------- Plots for Jo ----------------------------
    out_name=results_directory+"/lanczos_vs_PlanczosIF_Jo.png"
    compare_plot(out_name,lanczos[:,5],PlanczosIF[:,5],r'$J_o$',
                 diff[:,5],r'$\Delta J_o$',
                 itot,'iterations',io)

    # # Create figure window to plot data
    # fig = plt.figure(1, figsize=(9,9))
    # gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

    # # Top plot:
    # ax1 = fig.add_subplot(gs[0])
    # plt.yscale("log")
    # ax1.plot(itot[:],lanczos[:,5],color='black')
    # ax1.plot(itot[:],PlanczosIF[:,5],color='green')
    # # ax1.set_xlabel('iteratiosns')
    # ymax=max(max(lanczos[:,5]),max(PlanczosIF[:,5]))
    # ax1.vlines(io, 0, ymax, colors='black', linestyles='dashed')
    # ax1.set_ylabel('Jo')

    # # Bottom plot: residuals
    # ax2 = fig.add_subplot(gs[1])
    # ax2.plot(itot[:],diff[:,5],color='black')
    # ax2.axhline(color="gray", zorder=-1)
    # ax2.set_xlabel('iterations')
    # ax2.set_ylabel('deltaJo')
    # plt.savefig(out_name)
    # plt.clf()
################################################################################
