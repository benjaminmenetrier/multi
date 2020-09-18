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
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.4f'))
    plt.savefig(out_name)
    plt.clf()
################################################################################



################################################################################
def compare_plot_N(out_name,obj1,ylabel1,obj2,ylabel2,x,xlabel,io):
    """Produces usefull comparision with N plots:
    Args:
        out_name (string): Name of the output png file.
        obj1 (list)      : List containing lists to be compared together.
        yabel (string)   : Label of the y-axis.
        obj2 (list)      : List of numbers regarding which one wants to compare.
        ylabel (string)  : Label of obj2.
        x (list)         : List of numbers for x-axis.
        xlabel (string)  : Label of x-axis
    Return:
        (void): Creates the plot and save it in the png file named by out_name arg.
    """
    
    # Create figure window to plot data
    fig = plt.figure(1, figsize=(9,9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

    # Top plot: [obj11,obj12,..obj1N]
    ax1 = fig.add_subplot(gs[0])
    plt.yscale("log")
    ymax=0.
    for obj in obj1:
        ymax_tmp=max(obj)
        if (ymax_tmp>ymax):
            ymax=ymax_tmp
        c=np.random.random(3)
        ax1.plot(x[:],obj,color=c)
    ax1.set_ylabel(ylabel1)
    start, end = ax1.get_xlim()
    #ax1.xaxis.set_ticks(np.arange(start, end, 1.))
    ax1.vlines(io, 0, ymax, colors='black', linestyles='dashed')

    # Bottom plot: obj2
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(x[:],obj2,color='black')
    ax2.axhline(color="gray", zorder=-1)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel2)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.4f'))
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

    #--------------- Plots for Jb ----------------------------
    out_name=results_directory+"/lanczos_vs_PlanczosIF_Jb.png"
    compare_plot(out_name, lanczos[:,4],PlanczosIF[:,4],r'$J_b$',
                 diff[:,4],r'$\Delta J_b$',
                 itot,'iterations',io)

    #--------------- Plots for Jo ----------------------------
    out_name=results_directory+"/lanczos_vs_PlanczosIF_Jo.png"
    compare_plot(out_name,lanczos[:,5],PlanczosIF[:,5],r'$J_o$',
                 diff[:,5],r'$\Delta J_o$',
                 itot,'iterations',io)

    # test plot:
    obj1=[lanczos[:,5],PlanczosIF[:,3],lanczos[:,3]]
    ylabel1='lala'
    ylabel2='toto'
    obj2=diff[:,5]
    x=itot
    xlabel='iterations'
    io=io
    out_name=results_directory+"/test.png"

    compare_plot_N(out_name,obj1,ylabel1,obj2,ylabel2,x,xlabel,io)


################################################################################

