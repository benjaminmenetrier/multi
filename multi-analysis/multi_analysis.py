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
from matplotlib.offsetbox import AnchoredText

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
    ax1.plot(x[:len(obj1)],obj1,color='blue')
    ax1.plot(x[:len(obj2)],obj2,color='blue',linestyle='dashed')
    ymax=max(max(obj1),max(obj2))
    if ymax > 100:
        plt.yscale("log")
    ax1.set_ylabel(ylabel12)
    start, end = ax1.get_xlim()
    #ax1.xaxis.set_ticks(np.arange(start, end, 1.))
    ax1.vlines(io, 0, ymax, colors='blue', linestyles='dashed')

    # Bottom plot: obj3
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(x[:len(obj3)],obj3,color='black')
    ax2.axhline(color="gray", zorder=-1)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel3)
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
    #diff=np.genfromtxt(results_directory+'lanczos_control_vs_PlanczosIF_model.dat', comments='#')
    diff=[]
    l_iter=min(len(lanczos[0]),len(PlanczosIF[0]))
    for c in range(len(lanczos)):
        diff_col=[]
        for l in range(l_iter):
            diff_col.append(lanczos[c,l]-PlanczosIF[c,l]) 
        diff.append(diff_col)
    diff=np.array(diff)
    
    # lanczos=np.atleast_2d(np.genfromtxt(results_directory+'/lanczos_control_space.dat', comments='#'))
    # PlanczosIF=np.atleast_2d(np.genfromtxt(results_directory+'/PlanczosIF_model_space.dat', comments='#'))
    # diff=np.atleast_2d(np.genfromtxt(results_directory+'/lanczos_control_vs_PlanczosIF_model.dat', comments='#'))

    # Total iterations (outer x inner):
    itot=list(range(len(lanczos)))

    #--------------- Plots for J -----------------------------
    out_name=results_directory+"/lanczos_vs_PlanczosIF_J.png"
    compare_plots(out_name,lanczos[:,3],PlanczosIF[:,3],r'$J=J_b+J_o$',
                 diff[:,3],r'$J_{B^{1/2}}-J_{B}$',
                 itot,r'iterations',io)

    #--------------- Plots for Jb ----------------------------
    out_name=results_directory+"/lanczos_vs_PlanczosIF_Jb.png"
    compare_plots(out_name, lanczos[:,4],PlanczosIF[:,4],r'$J_b$',
                 diff[:,4],r'$J_{b,B^{1/2}}-J_{b,B}$',
                 itot,'iterations',io)

    #--------------- Plots for Jo ----------------------------
    out_name=results_directory+"/lanczos_vs_PlanczosIF_Jo.png"
    compare_plots(out_name,lanczos[:,5],PlanczosIF[:,5],r'$J_o$',
                 diff[:,5],r'$J_{o,B^{1/2}}-J_{o,B}$',
                 itot,'iterations',io)
################################################################################


################################################################################
def multi_plot(res_dir_list,outer_iterations):
    """Run the analysis over the results:
    """
    for r,res_dir in enumerate(res_dir_list):
        print("plotting for raw results for: \n",res_dir)
        try:
            evolution_plot(res_dir,outer_iterations[r])
        except:
            evolution_plot(res_dir,outer_iterations[r])
            #print("Error with directory:",res_dir)
################################################################################


################################################################################
def yo_vs_hxg_plot(res_dir):
    """Plot the 1D observations versus the result of applying the observation operator to the guess:
    """
    try:
        for name in ['lanczos_control_space_outer_vectors.dat','PlanczosIF_model_space_outer_vectors.dat']:
            results_file=res_dir+name
            out_name=results_file[:-4]+'_obs_vs_Hxg.png'
            outer_vectors=np.genfromtxt(results_file, comments='#')
            
            io=1
            hxg=[]
            yo=[]
            d=[]
            indices=[]

            hxg_io=[]
            yo_io=[]
            d_io=[]
            indices_io=[]
            
            for i  in range(len(outer_vectors)):
                if outer_vectors[i,0]==io:
                    hxg_io.append(outer_vectors[i,-3])
                    yo_io.append(outer_vectors[i,-2])
                    d_io.append(outer_vectors[i,-1])
                    indices_io.append(outer_vectors[i,1])
                else:
                    io=io+1
                    hxg.append(hxg_io)
                    yo.append(yo_io)
                    d.append(d_io)
                    indices.append(indices_io)
                    hxg_io=[]
                    yo_io=[]
                    d_io=[]
                    indices_io=[]
                    hxg_io.append(outer_vectors[i,-3])
                    yo_io.append(outer_vectors[i,-2])
                    d_io.append(outer_vectors[i,-1])
                    indices_io.append(outer_vectors[i,1])
            print('Plotting the innovation for :', results_file)   
            # Create figure window to plot data
            fig = plt.figure(1, figsize=(9,9))
            gs = gridspec.GridSpec(2, 1, height_ratios=[6, 2])

            for r in range(len(yo)):
                # Top plot: yo and hxg
                ax1 = fig.add_subplot(gs[0])
                ax1.plot(indices[r][:],hxg[r][:],color='blue',label=r'$H x_g$')
                ax1.plot(indices[r][:],yo[r][:],color='red',label=r'$y^o$')
                plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
                   ncol=2, mode="expand", borderaxespad=0.)

                # Bottom plot: innovation
                ax2 = fig.add_subplot(gs[1])
                ax2.plot(indices[r][:],d[r][:],color='black')
                ax2.set_ylabel(r'Innovation')
                plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
                plt.savefig(out_name)
                plt.clf()
    except:
        print('Error in yo_vs_hxg for: ',res_dir)
################################################################################


################################################################################
def vec_plot(results_file,column_of_interest,label,out_file_name):
    """Plot the outer vectors:
    """
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
    # Add the last outer vector.
    vec.append(vec_io)
    indices.append(indices_io)
        
    if not len(vec)==1:
        fig, subplots = plt.subplots(len(vec),1)
        for i, ax in enumerate(subplots):
            plt.ticklabel_format(axis="y", style="sci", scilimits=(-3,3))
            ax.plot(indices[i][:],vec[i][:],color='blue')
            ax.set_ylabel(label+'_{}$'.format(i))
            # at = AnchoredText(r"io={}".format(i),
            #           prop=dict(size=15), frameon=True,loc='upper left',)
            # at.patch.set_boxstyle("round,pad=0.,rounding_size=0.1")
            # ax.add_artist(at)
            
        plt.tight_layout()
        fig.align_ylabels(subplots[:])    
        plt.subplots_adjust(hspace=0.5)
        #fig.text(0., 0.5, r'$x_g$', va='center', rotation='vertical')
    else:
        fig=plt.figure()
        plt.plot(indices[0][:],vec[0][:],color='blue')
        plt.ylabel(label+r'_{}$'.format(0))
        plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

    plt.savefig(out_file_name)
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
    # at = AnchoredText(r"Figure 1",
    #               prop=dict(size=15), frameon=True,loc='upper left',)
    # at.patch.set_boxstyle("round,pad=0.,rounding_size=0.1")
    # ax1.add_artist(at)
    plt.yscale("log")
    ymax=0.
    for o,obj in enumerate(obj_list):
        ymax_tmp=max(max(obj[0]),max(obj[1]))
        if (ymax_tmp>ymax):
            ymax=ymax_tmp
        ax1.plot(x[:len(obj[0])],obj[0],color=colors[o],label=legend[o][0])
        ax1.plot(x[:len(obj[1])],obj[1],color=colors[o],linestyle='dashed',label=legend[o][1])
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=3, mode="expand", borderaxespad=0.)
    ax1.set_ylabel(ylabel1)
    #start, end = ax1.get_xlim()
    #ax1.xaxis.set_ticks(np.arange(start, end, 1.))
    ax1.vlines(io, 0, ymax, colors='blue', linestyles='dashed')

    # Bottom plot: diff_list
    ax2 = fig.add_subplot(gs[1])
    for d,diff in enumerate(diff_list):
        ax2.plot(x[:len(diff)],diff,color=colors[d])
    ax2.axhline(color="gray", zorder=-1)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel(ylabel2)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.4f'))
    plt.savefig(out_name)
    plt.clf()
################################################################################


################################################################################
def lmp_compare(out_names,lmp_to_compare,column_of_interest,ylabel1,ylabel2,outer_iterations_list):
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
                #diff=np.genfromtxt(res_dir+'lanczos_control_vs_PlanczosIF_model.dat', comments='#')
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
            compare_plots_2N(out_name,obj_list,ylabel1,diff_list,ylabel2,x,xlabel,outer_iterations,legend)
    except:
        print("Error with lmp comparision of:\n",res_dirs,"\n")
################################################################################


# def check_second_level_lmp(out_names,lmp_to_compare,outer_iterations_list):
#     """Compare spectral and ritz lmp modes:
#     """

#     legend=[]
#     for space in ['model','control']:
#         legend.append([space+'-ritz',space+'-spectral'])
        
#     for r,res_dirs in enumerate(lmp_to_compare):
#         print(res_dirs)
        
#         res_ritz1=np.genfromtxt(res_dirs[0]+'PlanczosIF_model_space.dat', comments='#')    
#         res_ritz2=np.genfromtxt(res_dirs[0]+'lanczos_control_space.dat', comments='#')

#         res_spec1=np.genfromtxt(res_dirs[1]+'PlanczosIF_model_space.dat', comments='#')    
#         res_spec2=np.genfromtxt(res_dirs[1]+'lanczos_control_space.dat', comments='#')

#         diff1=res_ritz1[:,3]-res_spec1[:,3]
#         diff2=res_ritz2[:,3]-res_spec2[:,3]
#         diff_list=[diff1,diff2]
#         print(diff_list)
        
#         obj1=[res_ritz1[:,3],res_spec1[:,3]]
#         obj2=[res_ritz2[:,3],res_spec2[:,3]]
#         obj_list=[obj1,obj2]

#         for o, obj in obj_list:
#             # maybe dirtyish but...
#             itot=list(range(len(res_ritz1)))
#             print("checkin second level lmp for:\n",out_names[r],"\n")
#             ylabel1=r'$J=J_o+J_b$'
#             ylabel2=r'$J_{ritz}-J_{spec}$'
#             x=itot
#             xlabel='iterations'
#             out_name=out_names[r]
#             outer_iterations=outer_iterations_list[r]

#             print(np.shape(obj_list))
#             try:
#                 compare_plots_2N(out_name,obj_list,ylabel1,diff_list,ylabel2,x,xlabel,outer_iterations,legend)
#             except:
#                 compare_plots_2N(out_name,obj_list,ylabel1,diff_list,ylabel2,x,xlabel,outer_iterations,legend)
#                 #print("Error with lmp comparision of:\n",res_dirs,"\n")



# ################################################################################
# def diff_plot(out_names,param_list,res_dir_list):

#     legend=[]
#     for lmp_mode in ['ritz','spectral','none']:
#             for par in param_list:
#                 legend.append([lmp_mode+'-model',lmp_mode+'-control'])

#     ylabel=r'$J_{B^{1/2}}-J_{B}$'
#     xlabel=r'\N_obs'
#     color_map=cm.get_cmap('copper', 3)
#     colors=color_map(range(len(param_list)))

#     for r,res_dirs in enumerate(res_dir_list):
#         diff_list=[]
#         for res_dir in res_dirs:
#             res=np.genfromtxt(res_dir+'lanczos_control_vs_PlanczosIF_model.dat', comments='#')
#             diff_list.append(res[-1,3])
            
#     # Create figure window to plot data:
#     fig = plt.figure(1, figsize=(9,9))

#     for par in param_list: 
#         ax1.plot(param_list[:],diff_list[:],color=colors[o],label=legend[o][0])
#     plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
#                                    ncol=3, mode="expand", borderaxespad=0.)
#     ax1.set_ylabel(ylabel)
#     ax1.set_xlabel(xlabel)
#     plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
#     plt.savefig(out_name)
#     plt.clf()
# ################################################################################

################################################################################
# matrix monitoring
def matrix_monitoring(res_dir,results_file_name,out_file_name):
    """Draw the B matrix in color code
    """
    res_file=res_dir+results_file_name
    bmatrix=np.genfromtxt(res_file,comments='#')

    b_vec=[]
    b_mat=[]
    bmatrices=[]
    io_old=1
    ib_old=1

    for b in bmatrix:
        io=b[0]
        ib=b[1]
        elem=b[2]
        #print(b)
        if io==io_old:
            if ib==ib_old:
                b_vec.append(elem)
                #print('vec1',b_vec,io,io_old,ib,ib_old)
            else:
                #print('vec21',b_vec,io,io_old,ib,ib_old)
                ib_old=ib
                b_mat.append(np.array(b_vec))
                b_vec=[]
                b_vec.append(elem)
                #print('vec22',b_vec,io,io_old,ib,ib_old)
                #print('mat1',b_mat)
        else:
            #print('matrix',b_mat,io,io_old,ib,ib_old)
            io_old=io
            ib_old=1
            b_mat.append(np.array(b_vec))
            #print('mat',b_mat,io,io_old,ib,ib_old)
            bmatrices.append(np.array(b_mat))
            b_mat=[]
            b_vec=[]
            b_vec.append(elem)
    b_mat.append(b_vec)        
    bmatrices.append(b_mat)  
    bmatrices=np.array(bmatrices)
    #print(bmatrices)
    #print(np.shape(bmatrices))
    
    for b,bmat in enumerate(bmatrices):
        #print('final',np.array(b_mat))
        fig=plt.figure()
        outname=res_dir+out_file_name+'_io{}.png'.format(b+1)
        plt.matshow(np.array(bmat),cmap=plt.get_cmap('copper'))
        plt.colorbar()
        plt.savefig(outname)
################################################################################
